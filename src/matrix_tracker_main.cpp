// matrix_tracker_main.cpp — CLI матричного ядра SynchrotronSim (C++).
//
// Цель: автономное матричное ядро — все расчёты в C++ (Python только БД/анализ/визуализация).
// Python-ядро matrix_tracker.py остаётся ЭТАЛОНОМ для разработки и валидации.
//
// Режимы:
//   --element <kind> <Leff> <B> <G> <Gv> <S> <O> <edge1> <edge2> <nsl> <design_mev> <sr> <x px y py t delta>
//       Прогнать ОДНО тело элемента над состоянием и напечатать результат (для регрессии против Python,
//       см. test_matrix_core_cpp.py). kind: dipole|dq|quad|sext|oct. B,G,Gv,S,O — физкоэф. mat_*.
//   <exp_id>   (в разработке) — автономный трек всего кольца из БД. Каркас/чтение БД — следующий слой;
//              физический движок (matrix_core.h) уже готов и провалидирован.
#include <iostream>
#include <string>
#include <sstream>
#include <cstdio>
#include <vector>
#include <filesystem>
#include <system_error>
#include "matrix_core.h"
#include "matrix_ring.h"
#include "matrix_npz.h"
#include "matrix_ring_serve.h"
#ifdef _WIN32
#include <io.h>
#include <fcntl.h>
#endif

#ifdef _WIN32
#include <windows.h>
#endif

// Поиск папки workspace (как src/main.cpp): вверх от .exe / argv0 / cwd до каталога с подпапкой workspace.
static std::string findWorkspace(const char* argv0) {
    namespace fs = std::filesystem; std::error_code ec; std::vector<fs::path> starts;
    #ifdef _WIN32
    char exePath[MAX_PATH]; DWORD len = GetModuleFileNameA(NULL, exePath, MAX_PATH);
    if (len > 0 && len < MAX_PATH) starts.push_back(fs::path(exePath).parent_path());
    #endif
    if (argv0 && *argv0) { fs::path p(argv0); if (p.has_parent_path()) starts.push_back(p.parent_path()); }
    starts.push_back(fs::current_path(ec));
    for (const auto& start : starts) {
        fs::path dir = fs::absolute(start, ec);
        for (int i = 0; i < 10; ++i) {
            if (fs::is_directory(dir / "workspace", ec)) return (dir / "workspace").generic_string() + "/";
            fs::path parent = dir.parent_path();
            if (parent.empty() || parent == dir) break;
            dir = parent;
        }
    }
    return "workspace/";
}

static void printState(const mc::State& s) {
    // Полная точность (%.17g) — чтобы сравнение с Python-эталоном было предельно точным.
    std::printf("%.17g %.17g %.17g %.17g %.17g %.17g\n", s.x, s.px, s.y, s.py, s.t, s.delta);
}

static int runElement(int argc, char* argv[]) {
    // argv[2..13] — параметры элемента (12), argv[14..19] — состояние (6). Итого argc >= 20.
    if (argc < 20) {
        std::cerr << "Usage: --element <kind> <Leff> <B> <G> <Gv> <S> <O> <edge1> <edge2>"
                     " <nsl> <design_mev> <sr> <x> <px> <y> <py> <t> <delta>" << std::endl;
        return 2;
    }
    std::string kind = argv[2];
    double Leff = std::stod(argv[3]);
    double B = std::stod(argv[4]), G = std::stod(argv[5]), Gv = std::stod(argv[6]);
    double S = std::stod(argv[7]), O = std::stod(argv[8]);
    double e1 = std::stod(argv[9]), e2 = std::stod(argv[10]);
    int nsl = std::stoi(argv[11]);
    double design = std::stod(argv[12]);
    bool sr = (std::stoi(argv[13]) != 0);
    mc::State st{ std::stod(argv[14]), std::stod(argv[15]), std::stod(argv[16]),
                  std::stod(argv[17]), std::stod(argv[18]), std::stod(argv[19]) };
    mc::Element el = mc::Element::make(kind, Leff, B, G, Gv, S, O, e1, e2, nsl, design, sr);
    try {
        printState(el.body(st));
    } catch (const mc::ParticleLost&) {
        std::printf("LOST\n");
    }
    return 0;
}

// --serve — ПЕРСИСТЕНТНЫЙ сервер кик-траекторий (замена Python track_matrix): читает запросы из stdin, пишет
// траектории в stdout. Без релонча на каждый фит-eval. Один запрос = одна строка:
//   <kind> <Leff> <M> <di> <do> <nsl> <design> <sr> <delta> <B> <G> <Gv> <S> <O> <e1> <e2> <nprobe> <u1 v1 …>
// Трек на пробу: дрейф di (8 подшагов) → элемент M суб-срезов длины Leff/M → дрейф do (8 подшагов), запись
// (s,x,y) в каждой точке (как track_matrix: 1 старт + 8 + M + 8). Ответ на пробу: "<npts> s x y s x y …"
// или "LOST". Строка "QUIT" завершает.
static int runServe() {
    std::ios::sync_with_stdio(false);
    std::string line;
    while (std::getline(std::cin, line)) {
        if (line.empty()) continue;
        if (line == "QUIT" || line == "quit") break;
        std::istringstream is(line);
        std::string kind; double Leff, di, doo, design, delta, B, G, Gv, S, O, e1, e2;
        int M, nsl, srI, nprobe;
        if (!(is >> kind >> Leff >> M >> di >> doo >> nsl >> design >> srI >> delta
                 >> B >> G >> Gv >> S >> O >> e1 >> e2 >> nprobe)) {
            std::printf("ERR parse\n"); std::fflush(stdout); continue;
        }
        double sL = Leff / M;                                       // длина суб-среза (как build_matrix_magnet L/M)
        mc::Element sub = mc::Element::make(kind, sL, B, G, Gv, S, O, e1, e2, nsl, design, srI != 0);
        std::vector<std::pair<double, double>> probes(nprobe);
        for (int i = 0; i < nprobe; ++i) is >> probes[i].first >> probes[i].second;
        const int nsub = 8;                                         // подшаги дрейфа (track_matrix drift_rec n=8)
        for (int i = 0; i < nprobe; ++i) {
            mc::State st{ probes[i].first, 0.0, probes[i].second, 0.0, 0.0, delta };
            std::vector<double> pts; pts.reserve((size_t)3 * (1 + 2 * nsub + M));
            auto rec = [&](double s) { pts.push_back(s); pts.push_back(st.x); pts.push_back(st.y); };
            double s = 0.0; bool lost = false;
            rec(0.0);
            try {
                for (int k = 0; k < nsub; ++k) { st = mc::drift(st, di / nsub); s += di / nsub; rec(s); }
                for (int m = 0; m < M; ++m)    { st = sub.body(st);             s += sL;          rec(s); }
                for (int k = 0; k < nsub; ++k) { st = mc::drift(st, doo / nsub); s += doo / nsub; rec(s); }
            } catch (const mc::ParticleLost&) { lost = true; }
            if (lost) { std::printf("LOST\n"); }
            else {
                std::printf("%d", (int)(pts.size() / 3));
                for (double val : pts) std::printf(" %.17g", val);
                std::printf("\n");
            }
        }
        std::fflush(stdout);
    }
    return 0;
}

// --track <exp_id> <turns> [x0_mm y0_mm] — АВТОНОМНЫЙ трек кольца из БД: инжекция + N оборотов + ВЧ.
// Пишет results/trajectories_local_exp<exp>.npz (mode=long) для traj_analyzer + печатает финальные состояния
// (для регрессии): "<idx> x px y py t delta surv" или "<idx> LOST surv". x0/y0 (мм) ⇒ ОДНА частица (скан).
static int runTrack(const std::string& ws, int exp, int turns, double x0_mm, double y0_mm) {
    mc::Table ex = mc::readTable(ws + "tables/experiments.csv");
    int ce = ex.col("exp_id"), cl = ex.col("lattice_id"), cg = ex.col("group_id");
    int clr = ex.col("loss_radius_mm"), cps = ex.col("plane_same_tol_mm"),
        cptl = ex.col("plane_parallel_tol"), cns = ex.col("n_slices");
    int lattice = -1, group = -1; mc::Settings stg;
    for (auto& r : ex.rows) if (mc::cell(r, ce) == std::to_string(exp)) {
        lattice = std::stoi(mc::cell(r, cl)); group = std::stoi(mc::cell(r, cg));
        if (!mc::cell(r, clr).empty())  stg.loss_radius_m     = std::stod(mc::cell(r, clr)) * 1e-3;
        if (!mc::cell(r, cps).empty())  stg.plane_same_tol_m  = std::stod(mc::cell(r, cps)) * 1e-3;
        if (!mc::cell(r, cptl).empty()) stg.plane_parallel_tol = std::stod(mc::cell(r, cptl));
        if (!mc::cell(r, cns).empty())  stg.n_slices          = std::stoi(mc::cell(r, cns));
        break;
    }
    if (lattice < 0) { std::cerr << "Experiment " << exp << " not found\n"; return 1; }
    // радиация: строка в radiation_params.csv → СИ вкл (энергия/ВЧ оттуда); нет строки → выкл.
    bool sr = false, rf_ideal = false; double design = 6000.0;
    mc::Table rp = mc::readTable(ws + "tables/radiation_params.csv");
    int re = rp.col("exp_id"), rd = rp.col("design_energy_mev"), rrf = rp.col("rf");
    for (auto& r : rp.rows) if (mc::cell(r, re) == std::to_string(exp)) {
        sr = true;
        if (!mc::cell(r, rd).empty()) design = std::stod(mc::cell(r, rd));
        rf_ideal = (mc::cell(r, rrf) == "ideal");
        break;
    }
    mc::Injection inj = mc::inject(ws, lattice, group, design, sr, stg);
    if (x0_mm != 0.0 || y0_mm != 0.0) {                        // ОДНА частица со смещением (скан апертуры/затухания)
        inj.states.clear();
        inj.states.push_back({ x0_mm * 1e-3, 0.0, y0_mm * 1e-3, 0.0, 0.0, 0.0 });
    }
    std::printf("# exp %d lattice %d group %d: %zu particles, %s, %zu elements, sr=%d rf=%s design=%.0f turns=%d\n",
                exp, lattice, group, inj.states.size(), inj.closed ? "closed" : "open",
                inj.ordered.size(), (int)sr, rf_ideal ? "ideal" : "none", design, turns);
    // трек + сбор turn-by-turn для .npz (mode=long). Частицы НЕЗАВИСИМЫ → OpenMP (кольцо read-only).
    std::vector<double> tx, txp, ty, typ;
    std::vector<int64_t> tt_counts, surv;
    std::vector<mc::LongResult> results(inj.states.size());
#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < (int)inj.states.size(); ++i)
        results[i] = mc::trackLong(inj, inj.states[i], turns, rf_ideal);
    for (size_t i = 0; i < results.size(); ++i) {              // сбор В ПОРЯДКЕ + вывод (серийно, детерминизм)
        mc::LongResult& lr = results[i];
        tt_counts.push_back((int64_t)lr.tt.size());
        surv.push_back(lr.surv_turns);
        for (auto& a : lr.tt) { tx.push_back(a[0]); txp.push_back(a[1]); ty.push_back(a[2]); typ.push_back(a[3]); }
        const mc::State& s = lr.final_state;
        if (lr.lost) std::printf("%zu LOST %d\n", i, lr.surv_turns);
        else std::printf("%zu %.17g %.17g %.17g %.17g %.17g %.17g %d\n",
                         i, s.x, s.px, s.y, s.py, s.t, s.delta, lr.surv_turns);
    }
    // .npz для traj_analyzer (mode=long: turn-by-turn у плоскости старта)
    mc::NpzWriter npz;
    npz.addVecI8("tt_counts", tt_counts);
    npz.addVecF8("tt_x", tx); npz.addVecF8("tt_xp", txp); npz.addVecF8("tt_y", ty); npz.addVecF8("tt_yp", typ);
    npz.addVecI8("surv_turns", surv);
    npz.addScalarStr("mode", "long");
    npz.addScalarI8("n_turns", turns);
    npz.addScalarBool("closed", inj.closed);
    npz.addScalarF8("ring_len", inj.ring_len);
    npz.addScalarI8("lattice", lattice); npz.addScalarI8("group", group);
    npz.addScalarI8("sr", sr ? 1 : 0); npz.addScalarI8("isr", 0);
    npz.addScalarStr("rf", rf_ideal ? "ideal" : "none");
    npz.addScalarF8("design_mev", design);
    std::error_code ec; std::filesystem::create_directories(ws + "results", ec);
    std::string npzPath = ws + "results/trajectories_local_exp" + std::to_string(exp) + ".npz";
    if (npz.write(npzPath)) std::printf("# npz: %s\n", npzPath.c_str());
    else std::fprintf(stderr, "npz write failed: %s\n", npzPath.c_str());
    return 0;
}

// Чтение конфига опыта (лат/группа/Settings/радиация) — общее для --track и --full.
struct ExpCfg { int lattice = -1, group = -1; mc::Settings stg; bool sr = false, rf_ideal = false; double design = 6000.0; };
static bool readExpCfg(const std::string& ws, int exp, ExpCfg& c) {
    mc::Table ex = mc::readTable(ws + "tables/experiments.csv");
    int ce = ex.col("exp_id"), cl = ex.col("lattice_id"), cg = ex.col("group_id");
    int clr = ex.col("loss_radius_mm"), cps = ex.col("plane_same_tol_mm"),
        cptl = ex.col("plane_parallel_tol"), cns = ex.col("n_slices");
    for (auto& r : ex.rows) if (mc::cell(r, ce) == std::to_string(exp)) {
        c.lattice = std::stoi(mc::cell(r, cl)); c.group = std::stoi(mc::cell(r, cg));
        if (!mc::cell(r, clr).empty())  c.stg.loss_radius_m      = std::stod(mc::cell(r, clr)) * 1e-3;
        if (!mc::cell(r, cps).empty())  c.stg.plane_same_tol_m   = std::stod(mc::cell(r, cps)) * 1e-3;
        if (!mc::cell(r, cptl).empty()) c.stg.plane_parallel_tol = std::stod(mc::cell(r, cptl));
        if (!mc::cell(r, cns).empty())  c.stg.n_slices           = std::stoi(mc::cell(r, cns));
        break;
    }
    if (c.lattice < 0) return false;
    mc::Table rp = mc::readTable(ws + "tables/radiation_params.csv");   // строка → СИ вкл (энергия/ВЧ оттуда)
    int re = rp.col("exp_id"), rd = rp.col("design_energy_mev"), rrf = rp.col("rf");
    for (auto& r : rp.rows) if (mc::cell(r, re) == std::to_string(exp)) {
        c.sr = true;
        if (!mc::cell(r, rd).empty()) c.design = std::stod(mc::cell(r, rd));
        c.rf_ideal = (mc::cell(r, rrf) == "ideal");
        break;
    }
    return true;
}

// --full <exp> <turns> [x0_mm y0_mm] — ПОЛНЫЙ трек (x(S)/спектр/тюны): запись на границах элементов за все
// обороты → results/trajectories_local_exp<exp>.npz (mode=full). БЕЗ ВЧ (как local_trajectories).
static int runFull(const std::string& ws, int exp, int turns, double x0_mm, double y0_mm) {
    ExpCfg c;
    if (!readExpCfg(ws, exp, c)) { std::cerr << "Experiment " << exp << " not found\n"; return 1; }
    const int TARGET = 30000;
    mc::Injection inj = mc::inject(ws, c.lattice, c.group, c.design, c.sr, c.stg);
    if (x0_mm != 0.0 || y0_mm != 0.0) {
        inj.states.clear(); inj.states.push_back({ x0_mm * 1e-3, 0.0, y0_mm * 1e-3, 0.0, 0.0, 0.0 });
    }
    std::vector<std::string> elem_names; elem_names.push_back("INJ_PLANE");
    for (auto& nd : inj.ordered) elem_names.push_back(nd.dev.empty() ? "?" : nd.dev);
    std::printf("# exp %d lattice %d group %d: %zu particles, %s, %zu elements, sr=%d design=%.0f turns=%d (mode=full)\n",
                exp, c.lattice, c.group, inj.states.size(), inj.closed ? "closed" : "open",
                inj.ordered.size(), (int)c.sr, c.design, turns);
    std::vector<double> S, x, px, y, py, ct, delta, X, Y, Z;
    std::vector<int64_t> elem, counts, surv;
    for (size_t i = 0; i < inj.states.size(); ++i) {
        mc::FullResult fr = mc::trackFull(inj, inj.states[i], turns, TARGET);
        counts.push_back((int64_t)fr.size()); surv.push_back(fr.surv_turns);
        S.insert(S.end(), fr.S.begin(), fr.S.end());       px.insert(px.end(), fr.px.begin(), fr.px.end());
        x.insert(x.end(), fr.x.begin(), fr.x.end());       y.insert(y.end(), fr.y.begin(), fr.y.end());
        py.insert(py.end(), fr.py.begin(), fr.py.end());   ct.insert(ct.end(), fr.ct.begin(), fr.ct.end());
        delta.insert(delta.end(), fr.delta.begin(), fr.delta.end());  elem.insert(elem.end(), fr.elem.begin(), fr.elem.end());
        X.insert(X.end(), fr.X.begin(), fr.X.end());       Y.insert(Y.end(), fr.Y.begin(), fr.Y.end());
        Z.insert(Z.end(), fr.Z.begin(), fr.Z.end());
        std::printf("%zu points=%zu surv=%d\n", i, fr.size(), fr.surv_turns);
    }
    mc::NpzWriter npz;
    npz.addVecI8("counts", counts);
    npz.addVecF8("S", S); npz.addVecF8("x", x); npz.addVecF8("xp", px); npz.addVecF8("y", y);
    npz.addVecF8("yp", py); npz.addVecF8("ct", ct); npz.addVecF8("delta", delta); npz.addVecI8("elem", elem);
    npz.addVecF8("X", X); npz.addVecF8("Y", Y); npz.addVecF8("Z", Z);
    npz.addVecStr("elem_names", elem_names);
    npz.addVecI8("surv_turns", surv);
    npz.addScalarStr("mode", "full");
    npz.addScalarI8("n_turns", turns);
    npz.addScalarBool("closed", inj.closed);
    npz.addScalarF8("ring_len", inj.ring_len);
    npz.addScalarI8("lattice", c.lattice); npz.addScalarI8("group", c.group);
    npz.addScalarI8("sr", c.sr ? 1 : 0); npz.addScalarI8("isr", 0);
    npz.addScalarStr("rf", "none"); npz.addScalarF8("design_mev", c.design);
    std::error_code ec; std::filesystem::create_directories(ws + "results", ec);
    std::string npzPath = ws + "results/trajectories_local_exp" + std::to_string(exp) + ".npz";
    if (npz.write(npzPath)) std::printf("# npz: %s\n", npzPath.c_str());
    else std::fprintf(stderr, "npz write failed: %s\n", npzPath.c_str());
    return 0;
}

// --ring-serve <lattice> <lead_m|nan> <design_mev> <sr> <isr> <nsl> <order> <rf> <volt_mv> <harmonic>
//              <phase_deg> <seed> — сервер кольца для анализа (порт ring_analysis.Ring), см. matrix_ring_serve.h.
static int runRingServe(int argc, char* argv[]) {
    if (argc < 14) {
        std::cerr << "Usage: --ring-serve <lattice> <lead_m|nan> <design_mev> <sr> <isr> <nsl> <order> <rf>"
                     " <volt_mv> <harmonic> <phase_deg> <seed>" << std::endl;
        return 2;
    }
    #ifdef _WIN32
    _setmode(_fileno(stdin), _O_BINARY); _setmode(_fileno(stdout), _O_BINARY);
    #endif
    mc::RingModel R;
    // необязательный 13-й аргумент (после seed) — путь к рабочей папке workspace/ с ВОЗМУЩЁННОЙ
    // lattice_configs (для ансамблей ошибок); по умолчанию — авто-поиск от пути exe.
    if (argc >= 15 && argv[14][0]) {
        std::string w = argv[14];
        if (!w.empty() && w.back() != '/' && w.back() != '\\') w += '/';
        R.ws = w;
    } else {
        R.ws = findWorkspace(argv[0]);
    }
    std::string lead = argv[3];
    R.cfg.lattice = std::stoi(argv[2]);
    R.cfg.lead = (lead == "nan" || lead == "none") ? std::numeric_limits<double>::quiet_NaN() : std::stod(lead);
    R.cfg.design_mev = std::stod(argv[4]);
    R.cfg.sr = std::stoi(argv[5]) != 0; R.cfg.isr = std::stoi(argv[6]) != 0;
    R.cfg.nsl = std::stoi(argv[7]); R.cfg.order = std::stoi(argv[8]);
    R.cfg.rf = std::stoi(argv[9]) != 0; R.cfg.volt_mv = std::stod(argv[10]);
    R.cfg.harmonic = std::stod(argv[11]); R.cfg.phase_deg = std::stod(argv[12]);
    R.cfg.seed = (uint64_t)std::stoull(argv[13]);
    try {
        return mc::serveRing(R);
    } catch (const std::exception& e) {
        std::fprintf(stdout, "ERR %s\n", e.what()); std::fflush(stdout);
        return 1;
    }
}

int main(int argc, char* argv[]) {
    #ifdef _WIN32
    SetConsoleOutputCP(CP_UTF8);
    #endif
    if (argc >= 2 && std::string(argv[1]) == "--ring-serve") {
        return runRingServe(argc, argv);
    }
    if (argc >= 2 && std::string(argv[1]) == "--element") {
        return runElement(argc, argv);
    }
    if (argc >= 2 && std::string(argv[1]) == "--serve") {
        return runServe();
    }
    if ((argc >= 4) && (std::string(argv[1]) == "--track" || std::string(argv[1]) == "--full")) {
        std::string ws = findWorkspace(argv[0]);
        double x0 = (argc >= 5) ? std::stod(argv[4]) : 0.0;
        double y0 = (argc >= 6) ? std::stod(argv[5]) : 0.0;
        int exp = std::stoi(argv[2]), turns = std::stoi(argv[3]);
        return (std::string(argv[1]) == "--full") ? runFull(ws, exp, turns, x0, y0)
                                                   : runTrack(ws, exp, turns, x0, y0);
    }
    std::cout << "MatrixTracker (SynchrotronSim matrix core, C++)\n"
                 "  --element <kind Leff B G Gv S O edge1 edge2 nsl design_mev sr> <x px y py t delta>\n"
                 "      прогнать тело одного элемента (валидация против Python-эталона)\n"
                 "  --track <exp_id> <turns> [x0_mm y0_mm]\n"
                 "      автономный трек кольца: turn-by-turn у плоскости + ВЧ → .npz (mode=long, затухание)\n"
                 "  --full <exp_id> <turns> [x0_mm y0_mm]\n"
                 "      полный трек: запись на границах элементов, x(S)/спектр/тюны → .npz (mode=full)\n"
                 "  --serve   персистентный сервер кик-траекторий (stdin→stdout): дрейф→элемент(M суб)→дрейф\n"
                 "      на пробу, замена Python track_matrix в element_compare (фит/RMS/демо)\n"
                 "  --ring-serve <lattice lead design_mev sr isr nsl order rf volt h phase seed>\n"
                 "      сервер кольца для анализа (ring_analysis/сверка с elegant): TRACK/TBT/ELEM/PIECES\n";
    return (argc < 2) ? 1 : 0;
}
