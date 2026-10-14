// matrix_ring_serve.h — СЕРВЕР КОЛЬЦА для анализа (ring_analysis / сверка с elegant) на C++-ядре.
//
// Дословный порт ring_analysis.Ring (Python): кольцо из решётки БД в её собственном порядке, плоскость
// наблюдения перед первым магнитом (start_before_first_m), оборот = плоскость → все магниты → плоскость,
// время t после оборота — запаздывание относительно света, ВЧ как RFCA elegant с FIDUCIAL="light".
// Поверх — пакетные запросы: трек N оборотов, пооборотные координаты у плоскости, координаты на выходе
// каждого магнита, координаты после каждого среза изгиба (для диффузии квантового излучения). Частицы
// независимы и считаются параллельно (OpenMP); у каждой свой генератор квантового возбуждения.
//
// Протокол (stdin/stdout, двоичный): строка-команда, затем полезная нагрузка — float64 little-endian.
//   TRACK n turns        + n·6 чисел     → OK m + m чисел (NaN у потерянных)
//   TBT n turns limit    + n·6           → OK m + n·(turns+1)·6 (NaN после потери; limit ≤ 0 — без порога)
//   ELEM n               + n·6           → OK m + n·(n_mag+1)·6 (выход каждого магнита + плоскость)
//   PIECES n             + n·6           → OK m + n·(n_pieces+1)·6 (старт + после каждого куска)
//   JUNCTION n           + 24 + n·6      → OK m + n·6 (стык плоскостей A→B: p, axis, u, v каждой)
//   PIECEMETA / INFO                     → OK k + k текстовых строк
//   KICK j theta | VKICK j theta | YAW j alpha | DK j dk | RESET   → OK 0
#pragma once
#include "matrix_ring.h"
#include <cstdio>
#include <cstdint>
#include <cstring>
#include <limits>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace mc {

struct RingConfig {
    int lattice = 1;
    double lead = std::numeric_limits<double>::quiet_NaN();   // NaN ⇒ половина замыкающего промежутка
    double design_mev = 6000.0;
    bool sr = false, isr = false;
    int nsl = SECTOR_NSLICES, order = 4;
    bool rf = false; double volt_mv = 0, harmonic = 0, phase_deg = 0;
    uint64_t seed = 12345;
};

struct Piece { int mag; bool junction; Element el; };

struct RingModel {
    RingConfig cfg;
    std::string ws;
    std::vector<Node> mags;
    Node plane;
    std::vector<double> s_in, s_out;
    double circumference = 0.0;
    std::vector<Piece> pieces;
    uint64_t requests = 0;

    void build() {
        Settings stg;
        mags = buildLine(ws, cfg.lattice, cfg.design_mev, cfg.sr || cfg.isr, stg, cfg.nsl, cfg.order, cfg.isr);
        if (mags.empty()) throw std::runtime_error("решётка пуста");
        const Node& m0 = mags.front(); const Node& m1 = mags.back();
        double gap = dot(m0.fr_in.p - m1.fr_out.p, m1.fr_out.axis);
        double d0 = std::isnan(cfg.lead) ? gap / 2.0 : cfg.lead;
        plane = makePlane(m0.fr_in.p - m0.fr_in.axis * d0, m0.fr_in.axis);
        s_in.clear(); s_out.clear();
        double s = d0;
        for (size_t j = 0; j < mags.size(); ++j) {
            if (j) s += dot(mags[j].fr_in.p - mags[j - 1].fr_out.p, mags[j - 1].fr_out.axis);
            s_in.push_back(s); s += mags[j].phys.Leff; s_out.push_back(s);
        }
        circumference = s + dot(plane.fr_in.p - m1.fr_out.p, m1.fr_out.axis);
        buildPieces();
    }

    // Изгиб → его срезы (копии с L/N, θ/N, одним срезом и краями только на крайних) — ring_analysis._bend_slices.
    void buildPieces() {
        pieces.clear();
        for (size_t j = 0; j < mags.size(); ++j) {
            const Element& m = mags[j].phys;
            if (m.sector) {
                int N = m.nsl;
                for (int i = 0; i < N; ++i) {
                    Element c = m;
                    c.Leff = m.Leff / N; c.theta0 = m.theta0 / N; c.nsl = 1;
                    c.kEdge1 = (i == 0) ? m.kEdge1 : 0.0;
                    c.kEdge2 = (i == N - 1) ? m.kEdge2 : 0.0;
                    c.precache();
                    pieces.push_back({ (int)j, i == 0, c });
                }
            } else {
                pieces.push_back({ (int)j, true, m });
            }
        }
    }

    double beta(double d) const {
        double pc = cfg.design_mev * (1.0 + d);
        return pc / std::sqrt(pc * pc + MC2_MEV * MC2_MEV);
    }

    // RFCA: ΔE = e·V·sin(φ + 2πh·t/C); канонические поперечные импульсы не меняются — Ring.rf_kick.
    double rfKick(double t, double d) const {
        const double PI = 3.14159265358979323846;
        double dE = cfg.volt_mv * std::sin(cfg.phase_deg * PI / 180.0 + 2.0 * PI * cfg.harmonic * t / circumference);
        double pc = cfg.design_mev * (1.0 + d);
        double E = std::sqrt(pc * pc + MC2_MEV * MC2_MEV) + dE;
        return std::sqrt(E * E - MC2_MEV * MC2_MEV) / cfg.design_mev - 1.0;
    }

    static bool finite(const State& s) {
        return std::isfinite(s.x) && std::isfinite(s.px) && std::isfinite(s.y) && std::isfinite(s.py)
            && std::isfinite(s.t) && std::isfinite(s.delta);
    }

    // Один оборот — Ring.turn. false ⇒ частица потеряна.
    bool turn(State& st, IsrRng* rng) const {
        double t0 = st.t;
        try {
            State s = modified_drift(st, plane.fr_out, mags[0].fr_in, mags[0].same_plane_as_prev);
            s = mags[0].body(s, rng);
            for (size_t j = 1; j < mags.size(); ++j) {
                s = modified_drift(s, mags[j - 1].fr_out, mags[j].fr_in, mags[j].same_plane_as_prev);
                s = mags[j].body(s, rng);
            }
            s = modified_drift(s, mags.back().fr_out, plane.fr_in, false);
            s.t = t0 + (s.t - t0) / beta(s.delta) - circumference;
            if (cfg.rf) s.delta = rfKick(s.t, s.delta);
            if (!finite(s)) return false;
            st = s;
            return true;
        } catch (const ParticleLost&) {
            return false;
        }
    }

    uint64_t particleSeed(size_t p) const {        // splitmix64 от (seed, номер запроса, номер частицы)
        uint64_t z = cfg.seed + 0x9E3779B97F4A7C15ULL * (requests * 1000003ULL + p + 1);
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
        return z ^ (z >> 31);
    }

    // ---- пакетные запросы ------------------------------------------------------------------------
    std::vector<double> track(const std::vector<State>& in, int turns) {
        ++requests;
        const double NaN = std::numeric_limits<double>::quiet_NaN();
        std::vector<double> out(in.size() * 6, NaN);
        #pragma omp parallel for schedule(dynamic)
        for (long p = 0; p < (long)in.size(); ++p) {
            IsrRng rng(particleSeed(p)); IsrRng* r = cfg.isr ? &rng : nullptr;
            State s = in[p]; bool ok = true;
            for (int k = 0; k < turns && ok; ++k) ok = turn(s, r);
            if (ok) { double* o = &out[p * 6]; o[0] = s.x; o[1] = s.px; o[2] = s.y; o[3] = s.py; o[4] = s.t; o[5] = s.delta; }
        }
        return out;
    }

    std::vector<double> tbt(const std::vector<State>& in, int turns, double limit) {
        ++requests;
        const double NaN = std::numeric_limits<double>::quiet_NaN();
        size_t row = (size_t)(turns + 1) * 6;
        std::vector<double> out(in.size() * row, NaN);
        #pragma omp parallel for schedule(dynamic)
        for (long p = 0; p < (long)in.size(); ++p) {
            IsrRng rng(particleSeed(p)); IsrRng* r = cfg.isr ? &rng : nullptr;
            State s = in[p];
            auto put = [&](int k, const State& q) {
                double* o = &out[p * row + (size_t)k * 6];
                o[0] = q.x; o[1] = q.px; o[2] = q.y; o[3] = q.py; o[4] = q.t; o[5] = q.delta;
            };
            put(0, s);
            for (int k = 0; k < turns; ++k) {
                if (!turn(s, r)) break;
                if (limit > 0 && (std::fabs(s.x) > limit || std::fabs(s.y) > limit)) break;
                put(k + 1, s);
            }
        }
        return out;
    }

    std::vector<double> elements(const std::vector<State>& in) {
        ++requests;
        const double NaN = std::numeric_limits<double>::quiet_NaN();
        size_t n = mags.size(), row = (n + 1) * 6;
        std::vector<double> out(in.size() * row, NaN);
        #pragma omp parallel for schedule(dynamic)
        for (long p = 0; p < (long)in.size(); ++p) {
            auto put = [&](size_t k, const State& q) {
                double* o = &out[p * row + k * 6];
                o[0] = q.x; o[1] = q.px; o[2] = q.y; o[3] = q.py; o[4] = q.t; o[5] = q.delta;
            };
            try {                                     // Ring.track_elements: выход каждого магнита, затем плоскость
                State s = modified_drift(in[p], plane.fr_out, mags[0].fr_in, mags[0].same_plane_as_prev);
                s = mags[0].body(s);
                if (!finite(s)) continue;
                put(0, s);
                bool ok = true;
                for (size_t j = 1; j < n && ok; ++j) {
                    s = modified_drift(s, mags[j - 1].fr_out, mags[j].fr_in, mags[j].same_plane_as_prev);
                    s = mags[j].body(s);
                    ok = finite(s);
                    if (ok) put(j, s);
                }
                if (ok) {
                    State e = modified_drift(s, mags.back().fr_out, plane.fr_in, false);
                    if (finite(e)) put(n, e);
                }
            } catch (const ParticleLost&) {}
        }
        return out;
    }

    std::vector<double> piecesTrack(const std::vector<State>& in) {
        ++requests;
        const double NaN = std::numeric_limits<double>::quiet_NaN();
        size_t n = pieces.size(), row = (n + 1) * 6;
        std::vector<double> out(in.size() * row, NaN);
        #pragma omp parallel for schedule(dynamic)
        for (long p = 0; p < (long)in.size(); ++p) {
            auto put = [&](size_t k, const State& q) {
                double* o = &out[p * row + k * 6];
                o[0] = q.x; o[1] = q.px; o[2] = q.y; o[3] = q.py; o[4] = q.t; o[5] = q.delta;
            };
            State s = in[p]; put(0, s);
            const Frame* prev = &plane.fr_out;
            try {
                for (size_t q = 0; q < n; ++q) {
                    const Piece& pc = pieces[q];
                    if (pc.junction) {
                        s = modified_drift(s, *prev, mags[pc.mag].fr_in, mags[pc.mag].same_plane_as_prev);
                        prev = &mags[pc.mag].fr_out;
                    }
                    s = pc.el.body(s);
                    if (!finite(s)) break;
                    put(q + 1, s);
                }
            } catch (const ParticleLost&) {}
        }
        return out;
    }

    // ---- возмущения (как ring_analysis: кик корректора, поворот магнита) --------------------------
    void kick(int j, double theta) {                 // доп. гор. кик; body() применяет его к любому элементу
        if (j < 0 || j >= (int)mags.size()) throw std::runtime_error("KICK: нет магнита");
        mags[j].phys.hkick = theta;                  // (включая секторные — для отклика DQ/диполей в orbit_infer)
    }
    void vkick(int j, double theta) {                // доп. верт. кик (корректор орбиты; и отклик секторных)
        if (j < 0 || j >= (int)mags.size()) throw std::runtime_error("VKICK: нет магнита");
        mags[j].phys.vkick = theta;
    }
    void dk(int j, double dq) {                       // BBA: транзитно изменить нормальный градиент k1 (обе плоскости)
        if (j < 0 || j >= (int)mags.size()) throw std::runtime_error("DK: нет магнита");
        mags[j].phys.k1 += dq; mags[j].phys.k1v += dq;
        mags[j].phys.precache();                      // у секторного обновит b1=k1·ρ (body_sector читает прекэш)
    }
    void dsk(int j, double dq) {                       // коррекция связи: транзитно изменить скью-градиент ks
        if (j < 0 || j >= (int)mags.size()) throw std::runtime_error("DSK: нет магнита");
        mags[j].phys.ks += dq;                         // киковый путь читает ks напрямую, прекэш не нужен
        mags[j].phys.precache();
    }
    void dk2(int j, double dq) {                       // нелинейный тюнер: транзитно изменить секступольную силу k2
        if (j < 0 || j >= (int)mags.size()) throw std::runtime_error("DK2: нет магнита");
        mags[j].phys.k2 += dq;                         // у секторного precache обновит b2=k2·ρ
        mags[j].phys.precache();
    }
    void dk3(int j, double dq) {                       // нелинейный тюнер: транзитно изменить октупольную силу k3
        if (j < 0 || j >= (int)mags.size()) throw std::runtime_error("DK3: нет магнита");
        mags[j].phys.k3 += dq;                         // у секторного precache обновит b3=k3·ρ
        mags[j].phys.precache();
    }
    void yaw(int j, double alpha) {                  // ring_analysis.yaw_element: поворот рамок вокруг v входа
        if (j < 0 || j >= (int)mags.size()) throw std::runtime_error("YAW: нет магнита");
        Node& m = mags[j];
        if (m.phys.sector) throw std::runtime_error("YAW: поворот секторного изгиба не поддержан");
        Vec3 P = m.fr_in.p, ax = m.fr_in.axis, u = m.fr_in.u, v = m.fr_in.v;
        Vec3 ax2 = rot(ax, v, alpha), u2 = rot(u, v, alpha);
        m.fr_in = { P, ax2, u2, v };
        m.fr_out = { P + ax2 * m.phys.Leff, ax2, u2, v };
    }
};

// ---- цикл сервера -----------------------------------------------------------------------------------
inline bool readStates(std::vector<State>& st, long n) {
    std::vector<double> buf((size_t)n * 6);
    if (n > 0 && std::fread(buf.data(), sizeof(double), buf.size(), stdin) != buf.size()) return false;
    st.resize((size_t)n);
    for (long p = 0; p < n; ++p) {
        const double* b = &buf[(size_t)p * 6];
        st[p] = { b[0], b[1], b[2], b[3], b[4], b[5] };
    }
    return true;
}

inline void sendDoubles(const std::vector<double>& v) {
    std::fprintf(stdout, "OK %zu\n", v.size());
    if (!v.empty()) std::fwrite(v.data(), sizeof(double), v.size(), stdout);
    std::fflush(stdout);
}

inline void sendText(const std::vector<std::string>& lines) {
    std::fprintf(stdout, "OK %zu\n", lines.size());
    for (auto& l : lines) std::fprintf(stdout, "%s\n", l.c_str());
    std::fflush(stdout);
}

inline int serveRing(RingModel& R) {
    R.build();
    std::fprintf(stdout, "READY %zu %.17g %zu\n", R.mags.size(), R.circumference, R.pieces.size());
    std::fflush(stdout);
    char line[512];
    std::vector<State> st;
    while (std::fgets(line, sizeof(line), stdin)) {
        std::istringstream is(line); std::string cmd; is >> cmd;
        try {
            if (cmd == "QUIT" || cmd.empty()) break;
            if (cmd == "TRACK") {
                long n; int turns; is >> n >> turns;
                if (!readStates(st, n)) break;
                sendDoubles(R.track(st, turns));
            } else if (cmd == "TBT") {
                long n; int turns; double limit; is >> n >> turns >> limit;
                if (!readStates(st, n)) break;
                sendDoubles(R.tbt(st, turns, limit));
            } else if (cmd == "ELEM") {
                long n; is >> n;
                if (!readStates(st, n)) break;
                sendDoubles(R.elements(st));
            } else if (cmd == "PIECES") {
                long n; is >> n;
                if (!readStates(st, n)) break;
                sendDoubles(R.piecesTrack(st));
            } else if (cmd == "JUNCTION") {              // стык двух произвольных плоскостей (проверка матриц стыков)
                long n; is >> n;
                double fb[24];
                if (std::fread(fb, sizeof(double), 24, stdin) != 24) break;
                if (!readStates(st, n)) break;
                Frame A{ {fb[0], fb[1], fb[2]}, {fb[3], fb[4], fb[5]}, {fb[6], fb[7], fb[8]}, {fb[9], fb[10], fb[11]} };
                Frame B{ {fb[12], fb[13], fb[14]}, {fb[15], fb[16], fb[17]}, {fb[18], fb[19], fb[20]}, {fb[21], fb[22], fb[23]} };
                std::vector<double> out(st.size() * 6, std::numeric_limits<double>::quiet_NaN());
                for (size_t p = 0; p < st.size(); ++p) {
                    try {
                        State q = modified_drift(st[p], A, B, false);
                        double* o = &out[p * 6]; o[0] = q.x; o[1] = q.px; o[2] = q.y; o[3] = q.py; o[4] = q.t; o[5] = q.delta;
                    } catch (const ParticleLost&) {}
                }
                sendDoubles(out);
            } else if (cmd == "INFO") {
                std::vector<std::string> out;
                char b[512];
                for (size_t j = 0; j < R.mags.size(); ++j) {
                    const Element& e = R.mags[j].phys;
                    std::snprintf(b, sizeof(b), "%zu %.17g %.17g %.17g %s %s", j, R.s_in[j], R.s_out[j], e.Leff,
                                  R.mags[j].kind.empty() ? "-" : R.mags[j].kind.c_str(),
                                  R.mags[j].dev.empty() ? "-" : R.mags[j].dev.c_str());
                    out.push_back(b);
                }
                sendText(out);
            } else if (cmd == "PIECEMETA") {
                std::vector<std::string> out;
                char b[512];
                for (auto& pc : R.pieces) {
                    std::snprintf(b, sizeof(b), "%d %d %d %.17g %.17g %.17g %.17g %.17g", pc.mag, pc.junction ? 1 : 0,
                                  pc.el.sector ? 1 : 0, pc.el.Leff, pc.el.theta0, pc.el.k1, pc.el.k2, pc.el.k3);
                    out.push_back(b);
                }
                sendText(out);
            } else if (cmd == "KICK") {
                int j; double th; is >> j >> th; R.kick(j, th); sendText({});
            } else if (cmd == "VKICK") {
                int j; double th; is >> j >> th; R.vkick(j, th); sendText({});
            } else if (cmd == "YAW") {
                int j; double a; is >> j >> a; R.yaw(j, a); sendText({});
            } else if (cmd == "DK") {                        // транзитно изменить градиент квадруполя (для BBA)
                int j; double dk; is >> j >> dk; R.dk(j, dk); sendText({});
            } else if (cmd == "DSK") {                       // транзитно изменить скью-градиент (коррекция связи)
                int j; double dk; is >> j >> dk; R.dsk(j, dk); sendText({});
            } else if (cmd == "DK2") {                       // транзитно изменить секступольную силу (нелин. тюнер)
                int j; double dk; is >> j >> dk; R.dk2(j, dk); sendText({});
            } else if (cmd == "DK3") {                       // транзитно изменить октупольную силу (нелин. тюнер)
                int j; double dk; is >> j >> dk; R.dk3(j, dk); sendText({});
            } else if (cmd == "RESET") {
                R.build(); sendText({});
            } else {
                std::fprintf(stdout, "ERR неизвестная команда %s\n", cmd.c_str()); std::fflush(stdout);
            }
        } catch (const std::exception& e) {
            std::fprintf(stdout, "ERR %s\n", e.what()); std::fflush(stdout);
        }
    }
    return 0;
}

} // namespace mc
