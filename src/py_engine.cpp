// py_engine.cpp — модуль расширения CPython `synchrotronsim._engine` (nanobind) над ядром кольца.
//
// Шаги 0 и 2 плана docs/LIBRARY_ARCHITECTURE.md: то же ядро, что у сервера MatrixTracker --ring-serve
// (mc::RingModel из matrix_ring_serve.h), но в процессе Python — без канала stdin/stdout. Решётка — либо
// из таблиц ws (как сервер), либо из массивов Python-модели Machine (spec + cavities: диск не читается ни при
// сборке, ни при RESET). Снимок/откат (snapshot/restore/reset_to) и клон (clone) — копии mc::RingModel.
//
// Правила границы:
//   * каждый вызов трекинга — пакет частиц; на время счёта GIL отпущен, частицы параллельно (OpenMP);
//   * вход копируется в std::vector<State> ДО отпуска GIL, numpy-выход создаётся ПОСЛЕ — без копии
//     (буфер передаётся массиву вместе с владельцем-капсулой);
//   * у кольца свой мьютекс: одно кольцо из двух потоков Python — по очереди; разные кольца — параллельно.
// Результаты бит-в-бит совпадают с сервером при той же последовательности запросов (номер запроса входит в
// зерно ISR — RingModel::particleSeed).
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/map.h>
#include <nanobind/stl/tuple.h>

#include <mutex>
#include <tuple>

#include "matrix_ring_serve.h"

namespace nb = nanobind;
using namespace nb::literals;

namespace {

using StatesIn = nb::ndarray<const double, nb::shape<-1, 6>, nb::c_contig, nb::device::cpu>;
using FramesIn = nb::ndarray<const double, nb::shape<12>, nb::c_contig, nb::device::cpu>;
using NpD = nb::ndarray<nb::numpy, double>;

// std::vector → numpy без копии: буфер живёт, пока жив массив (владелец — капсула).
NpD to_numpy(std::vector<double>&& v, std::initializer_list<size_t> shape) {
    auto* p = new std::vector<double>(std::move(v));
    nb::capsule owner(p, [](void* q) noexcept { delete static_cast<std::vector<double>*>(q); });
    return NpD(p->data(), shape, owner);
}

std::vector<mc::State> to_states(const StatesIn& a) {
    const size_t n = a.shape(0);
    const double* d = a.data();
    std::vector<mc::State> st(n);
    for (size_t p = 0; p < n; ++p) {
        const double* b = d + p * 6;
        st[p] = { b[0], b[1], b[2], b[3], b[4], b[5] };
    }
    return st;
}

mc::Frame to_frame(const FramesIn& a) {
    const double* f = a.data();
    return { {f[0], f[1], f[2]}, {f[3], f[4], f[5]}, {f[6], f[7], f[8]}, {f[9], f[10], f[11]} };
}

using Col = nb::ndarray<const double, nb::ndim<1>, nb::c_contig, nb::device::cpu>;

std::vector<double> col_of(const nb::dict& d, const char* key, size_t n) {
    if (!d.contains(key)) throw std::runtime_error(std::string("spec: нет столбца ") + key);
    Col a = nb::cast<Col>(d[key]);
    if (a.shape(0) != n) throw std::runtime_error(std::string("spec: столбец ") + key + " другой длины");
    return std::vector<double>(a.data(), a.data() + n);
}

// Описание линии из словаря столбцов (synchrotronsim.engine.line_spec(Machine)): обмотки уже сложены.
std::shared_ptr<const mc::LineSpec> spec_of(const nb::dict& d) {
    auto L = std::make_shared<mc::LineSpec>();
    L->dev = nb::cast<std::vector<std::string>>(d["dev"]);
    L->kind = nb::cast<std::vector<std::string>>(d["kind"]);
    const size_t n = L->dev.size();
    L->x = col_of(d, "x", n); L->y = col_of(d, "y", n); L->z = col_of(d, "z", n);
    L->yaw = col_of(d, "yaw", n); L->pitch = col_of(d, "pitch", n); L->roll = col_of(d, "roll", n);
    L->leff = col_of(d, "leff", n);
    L->B = col_of(d, "B", n); L->G = col_of(d, "G", n); L->Gv = col_of(d, "Gv", n);
    L->S = col_of(d, "S", n); L->O = col_of(d, "O", n);
    L->e1 = col_of(d, "e1", n); L->e2 = col_of(d, "e2", n); L->lf1 = col_of(d, "lf1", n); L->lf2 = col_of(d, "lf2", n);
    L->nsl_box = nb::cast<std::vector<int>>(d["nsl_box"]);
    L->check();
    return L;
}

// Резонаторы: [(id, x, y, z, nx, ny, nz, volt_share)] — как CavitySpec из rf_cavities.csv.
using CavTuple = std::tuple<std::string, double, double, double, double, double, double, double>;
std::shared_ptr<const std::vector<mc::CavitySpec>> cavs_of(const std::vector<CavTuple>& cs) {
    auto out = std::make_shared<std::vector<mc::CavitySpec>>();
    for (const auto& c : cs) {
        mc::CavitySpec s; s.id = std::get<0>(c);
        s.pl.c = { std::get<1>(c), std::get<2>(c), std::get<3>(c) };
        s.pl.n = { std::get<4>(c), std::get<5>(c), std::get<6>(c) };
        if (!(mc::norm(s.pl.n) > 0)) throw std::runtime_error("резонатор " + s.id + ": нулевая нормаль");
        s.share = std::get<7>(c);
        out->push_back(s);
    }
    return out;
}

// Снимок живого состояния кольца (магниты с возмущениями, срезы, ВЧ) — копия модели; номер запроса не входит.
struct Snapshot { mc::RingModel R; };

struct Ring {
    mc::RingModel R;
    std::mutex mu;

    void configure(int lattice, double lead, double design_mev, bool sr, bool isr, int nsl, int order,
                   bool rf, double volt_mv, double harmonic, double phase_deg, uint64_t seed,
                   std::map<std::string, int> nsl_dev, bool nsl_dev_explicit) {
        R.cfg.lattice = lattice; R.cfg.lead = lead; R.cfg.design_mev = design_mev;
        R.cfg.sr = sr; R.cfg.isr = isr; R.cfg.nsl = nsl; R.cfg.order = order;
        R.cfg.rf = rf; R.cfg.volt_mv = volt_mv; R.cfg.harmonic = harmonic; R.cfg.phase_deg = phase_deg;
        R.cfg.seed = seed;
        R.nsl_dev = std::move(nsl_dev); R.nsl_dev_explicit = nsl_dev_explicit;
    }

    // из таблиц рабочей папки ws (как сервер)
    Ring(int lattice, std::string ws, double lead, double design_mev, bool sr, bool isr, int nsl, int order,
         bool rf, double volt_mv, double harmonic, double phase_deg, uint64_t seed,
         std::map<std::string, int> nsl_dev, bool nsl_dev_explicit) {
        if (!ws.empty() && ws.back() != '/' && ws.back() != '\\') ws += '/';
        R.ws = ws;
        configure(lattice, lead, design_mev, sr, isr, nsl, order, rf, volt_mv, harmonic, phase_deg, seed,
                  std::move(nsl_dev), nsl_dev_explicit);
        nb::gil_scoped_release nogil;
        R.build();
    }

    // из массивов Python-модели (шаг 2): диск не читается ни при сборке, ни при RESET
    Ring(nb::dict spec, std::vector<CavTuple> cavities, int lattice, double lead, double design_mev, bool sr,
         bool isr, int nsl, int order, bool rf, double volt_mv, double harmonic, double phase_deg, uint64_t seed,
         std::map<std::string, int> nsl_dev, bool nsl_dev_explicit) {
        R.spec = spec_of(spec);
        R.cav_spec = cavs_of(cavities);
        configure(lattice, lead, design_mev, sr, isr, nsl, order, rf, volt_mv, harmonic, phase_deg, seed,
                  std::move(nsl_dev), nsl_dev_explicit);
        nb::gil_scoped_release nogil;
        R.build();
    }

    explicit Ring(const mc::RingModel& m) : R(m) {}

    Ring* clone() {                                  // независимая копия: магниты, возмущения, ВЧ, номер запроса
        std::lock_guard<std::mutex> g(mu);
        return new Ring(R);
    }
    Snapshot* snapshot() {
        std::lock_guard<std::mutex> g(mu);
        return new Snapshot{ R };
    }
    void restore(const Snapshot& s) {                // номер запроса (зерно ISR) не откатывается
        std::lock_guard<std::mutex> g(mu);
        uint64_t rq = R.requests;
        R = s.R;
        R.requests = rq;
    }
    // RESET без пересборки: магниты и срезы — из снимка свежего кольца, а ВЧ (cfg) — текущий, как у
    // команды RESET сервера (параметры ВЧ живут в конфиге и переживают RESET).
    void reset_to(const Snapshot& s) {
        std::lock_guard<std::mutex> g(mu);
        uint64_t rq = R.requests;
        mc::RingConfig cfg = R.cfg;
        R = s.R;
        R.requests = rq; R.cfg = cfg;
        R.updateRf();
    }

    // ---- пакетный трекинг (GIL отпущен на время счёта) -------------------------------------------
    NpD track(const StatesIn& a, int turns) {
        auto st = to_states(a); std::vector<double> out;
        { nb::gil_scoped_release nogil; std::lock_guard<std::mutex> g(mu); out = R.track(st, turns); }
        return to_numpy(std::move(out), { st.size(), 6 });
    }
    NpD tbt(const StatesIn& a, int turns, double limit) {
        auto st = to_states(a); std::vector<double> out;
        { nb::gil_scoped_release nogil; std::lock_guard<std::mutex> g(mu); out = R.tbt(st, turns, limit); }
        return to_numpy(std::move(out), { st.size(), (size_t)turns + 1, 6 });
    }
    NpD elements(const StatesIn& a) {
        auto st = to_states(a); std::vector<double> out; size_t n = 0;
        { nb::gil_scoped_release nogil; std::lock_guard<std::mutex> g(mu); out = R.elements(st); n = R.mags.size(); }
        return to_numpy(std::move(out), { st.size(), n + 1, 6 });
    }
    NpD pieces(const StatesIn& a) {
        auto st = to_states(a); std::vector<double> out; size_t n = 0;
        { nb::gil_scoped_release nogil; std::lock_guard<std::mutex> g(mu);
          out = R.piecesTrack(st); n = R.currentPieces().size(); }
        return to_numpy(std::move(out), { st.size(), n + 1, 6 });
    }
    NpD junction(const StatesIn& a, const FramesIn& fa, const FramesIn& fb) {   // как команда JUNCTION
        auto st = to_states(a); mc::Frame A = to_frame(fa), B = to_frame(fb);
        std::vector<double> out(st.size() * 6, std::numeric_limits<double>::quiet_NaN());
        {
            nb::gil_scoped_release nogil; std::lock_guard<std::mutex> g(mu);
            for (size_t p = 0; p < st.size(); ++p) {
                try {
                    mc::State q = mc::modified_drift(st[p], A, B, false);
                    double* o = &out[p * 6]; o[0] = q.x; o[1] = q.px; o[2] = q.y; o[3] = q.py; o[4] = q.t; o[5] = q.delta;
                } catch (const mc::ParticleLost&) {}
            }
        }
        return to_numpy(std::move(out), { st.size(), 6 });
    }

    // ---- полное описание узлов (шаг 3): ТЕКУЩЕЕ состояние ядра — рамки, силы, флаги тела ------------
    // Плоский массив RingModel::describeFlat (раскладка — в matrix_ring_serve.h, разбор — ring_cpp.parse_describe);
    // та же функция отвечает на команду DESCRIBE сервера.
    NpD describe_flat() {
        std::vector<double> out;
        { std::lock_guard<std::mutex> g(mu); out = R.describeFlat(); }
        size_t n = out.size();
        return to_numpy(std::move(out), { n });
    }

    // ---- описание (как INFO / PIECEMETA сервера, но числами, без текста) ----------------------------
    nb::dict info() {
        std::lock_guard<std::mutex> g(mu);
        size_t n = R.mags.size();
        std::vector<double> s_in(R.s_in), s_out(R.s_out), leff(n);
        std::vector<std::string> kinds(n), names(n);
        for (size_t j = 0; j < n; ++j) {
            leff[j] = R.mags[j].phys.Leff; kinds[j] = R.mags[j].kind; names[j] = R.mags[j].dev;
        }
        nb::dict d;
        d["s_in"] = to_numpy(std::move(s_in), { n });
        d["s_out"] = to_numpy(std::move(s_out), { n });
        d["leff"] = to_numpy(std::move(leff), { n });
        d["kinds"] = kinds; d["names"] = names;
        return d;
    }
    std::vector<std::tuple<int, bool, bool, double, double, double, double, double>> piece_meta() {
        std::lock_guard<std::mutex> g(mu);
        std::vector<std::tuple<int, bool, bool, double, double, double, double, double>> out;
        for (const auto& pc : R.currentPieces())
            out.emplace_back(pc.mag, pc.junction, pc.el.sector, pc.el.Leff, pc.el.theta0, pc.el.k1, pc.el.k2, pc.el.k3);
        return out;
    }

    // ---- живые возмущения (как команды KICK..BSCALE, RF, RESET) -----------------------------------
    void kick(int j, double v)   { std::lock_guard<std::mutex> g(mu); R.kick(j, v); }
    void vkick(int j, double v)  { std::lock_guard<std::mutex> g(mu); R.vkick(j, v); }
    void yaw(int j, double v)    { std::lock_guard<std::mutex> g(mu); R.yaw(j, v); }
    void dk(int j, double v)     { std::lock_guard<std::mutex> g(mu); R.dk(j, v); }
    void dkv(int j, double v)    { std::lock_guard<std::mutex> g(mu); R.dkv(j, v); }
    void dsk(int j, double v)    { std::lock_guard<std::mutex> g(mu); R.dsk(j, v); }
    void dk2(int j, double v)    { std::lock_guard<std::mutex> g(mu); R.dk2(j, v); }
    void dk3(int j, double v)    { std::lock_guard<std::mutex> g(mu); R.dk3(j, v); }
    void bscale(int j, double v) { std::lock_guard<std::mutex> g(mu); R.bscale(j, v); }
    void set_rf(bool on, double volt_mv, double harmonic, double phase_deg, double freq_hz) {
        std::lock_guard<std::mutex> g(mu); R.setRf(on, volt_mv, harmonic, phase_deg, freq_hz);
    }
    void reset() { nb::gil_scoped_release nogil; std::lock_guard<std::mutex> g(mu); R.build(); }
};

int omp_threads() {
#ifdef _OPENMP
    return omp_get_max_threads();
#else
    return 1;
#endif
}

} // namespace

NB_MODULE(_engine, m) {
    m.doc() = "SynchrotronSim: C++-ядро кольца (mc::RingModel) как модуль расширения; шаг 0 плана библиотеки";
    m.def("omp_threads", &omp_threads, "Число потоков OpenMP, доступных ядру");
    m.attr("DESCRIBE_COLS") = nb::str(mc::DESCRIBE_COLS);    // раскладка Ring.describe_flat (= команда DESCRIBE)
    m.attr("build_info") = nb::str(
#if defined(__clang__)
        "clang " __clang_version__
#elif defined(__GNUC__)
        "gcc " __VERSION__
#elif defined(_MSC_VER)
        "msvc"
#else
        "unknown"
#endif
#ifdef _OPENMP
        ", openmp"
#endif
    );

    nb::class_<Snapshot>(m, "RingSnapshot", "Снимок живого состояния кольца (Ring.snapshot / Ring.restore)");

    nb::class_<Ring>(m, "Ring", "Кольцо (порт ring_analysis.Ring) — то же ядро, что MatrixTracker --ring-serve")
        .def(nb::init<int, std::string, double, double, bool, bool, int, int, bool, double, double, double,
                      uint64_t, std::map<std::string, int>, bool>(),
             "lattice"_a, "ws"_a, "lead"_a, "design_mev"_a, "sr"_a, "isr"_a, "nsl"_a, "order"_a,
             "rf"_a, "volt_mv"_a, "harmonic"_a, "phase_deg"_a, "seed"_a, "nsl_dev"_a, "nsl_dev_explicit"_a)
        .def(nb::init<nb::dict, std::vector<CavTuple>, int, double, double, bool, bool, int, int, bool, double,
                      double, double, uint64_t, std::map<std::string, int>, bool>(),
             "spec"_a, "cavities"_a, "lattice"_a, "lead"_a, "design_mev"_a, "sr"_a, "isr"_a, "nsl"_a, "order"_a,
             "rf"_a, "volt_mv"_a, "harmonic"_a, "phase_deg"_a, "seed"_a, "nsl_dev"_a, "nsl_dev_explicit"_a)
        .def("clone", &Ring::clone, nb::rv_policy::take_ownership)
        .def("snapshot", &Ring::snapshot, nb::rv_policy::take_ownership)
        .def("restore", &Ring::restore, "snapshot"_a)
        .def("reset_to", &Ring::reset_to, "snapshot"_a)
        .def_prop_ro("from_spec", [](Ring& r) { return (bool)r.R.spec; })
        .def("track", &Ring::track, "states"_a, "turns"_a)
        .def("track_tbt", &Ring::tbt, "states"_a, "turns"_a, "limit"_a)
        .def("track_elements", &Ring::elements, "states"_a)
        .def("track_pieces", &Ring::pieces, "states"_a)
        .def("junction", &Ring::junction, "states"_a, "A"_a, "B"_a)
        .def("info", &Ring::info)
        .def("describe_flat", &Ring::describe_flat)
        .def("piece_meta", &Ring::piece_meta)
        .def("kick", &Ring::kick, "j"_a, "theta"_a)
        .def("vkick", &Ring::vkick, "j"_a, "theta"_a)
        .def("yaw", &Ring::yaw, "j"_a, "alpha"_a)
        .def("dk", &Ring::dk, "j"_a, "dk"_a)
        .def("dkv", &Ring::dkv, "j"_a, "dk"_a)
        .def("dsk", &Ring::dsk, "j"_a, "dk"_a)
        .def("dk2", &Ring::dk2, "j"_a, "dk"_a)
        .def("dk3", &Ring::dk3, "j"_a, "dk"_a)
        .def("bscale", &Ring::bscale, "j"_a, "f"_a)
        .def("set_rf", &Ring::set_rf, "on"_a, "volt_mv"_a, "harmonic"_a, "phase_deg"_a, "freq_hz"_a)
        .def("reset", &Ring::reset)
        .def_prop_ro("n_mag", [](Ring& r) { return r.R.mags.size(); })
        .def_prop_ro("n_pieces", [](Ring& r) { return r.R.pieces.size(); })
        .def_prop_ro("circumference", [](Ring& r) { return r.R.circumference; })
        .def_prop_ro("c_rf", [](Ring& r) { return r.R.c_rf; })
        .def_prop_ro("requests", [](Ring& r) { return r.R.requests; });
}
