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
#include <cstdio>
#include <vector>
#include <filesystem>
#include <system_error>
#include "matrix_core.h"
#include "matrix_ring.h"
#include "matrix_npz.h"

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
    // трек + сбор turn-by-turn для .npz (mode=long)
    std::vector<double> tx, txp, ty, typ;
    std::vector<int64_t> tt_counts, surv;
    for (size_t i = 0; i < inj.states.size(); ++i) {
        mc::LongResult lr = mc::trackLong(inj, inj.states[i], turns, rf_ideal);
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

int main(int argc, char* argv[]) {
    #ifdef _WIN32
    SetConsoleOutputCP(CP_UTF8);
    #endif
    if (argc >= 2 && std::string(argv[1]) == "--element") {
        return runElement(argc, argv);
    }
    if (argc >= 4 && std::string(argv[1]) == "--track") {
        std::string ws = findWorkspace(argv[0]);
        double x0 = (argc >= 5) ? std::stod(argv[4]) : 0.0;
        double y0 = (argc >= 6) ? std::stod(argv[5]) : 0.0;
        return runTrack(ws, std::stoi(argv[2]), std::stoi(argv[3]), x0, y0);
    }
    std::cout << "MatrixTracker (SynchrotronSim matrix core, C++)\n"
                 "  --element <kind Leff B G Gv S O edge1 edge2 nsl design_mev sr> <x px y py t delta>\n"
                 "      прогнать тело одного элемента (валидация против Python-эталона)\n"
                 "  --track <exp_id> <turns> [x0_mm y0_mm]\n"
                 "      автономный трек кольца из БД (инжекция + N об. + ВЧ по radiation_params.csv) →\n"
                 "      results/trajectories_local_exp<exp>.npz для traj_analyzer; x0/y0 ⇒ одна частица (скан)\n";
    return (argc < 2) ? 1 : 0;
}
