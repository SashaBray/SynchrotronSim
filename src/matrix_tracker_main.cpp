// matrix_tracker_main.cpp — CLI матричного ядра SynchrotronSim (C++).
//
// Цель: автономное матричное ядро — все расчёты в C++ (Python только БД/анализ/визуализация).
// Python-ядро matrix_tracker.py остаётся ЭТАЛОНОМ для разработки и валидации.
//
// Режимы:
//   --element <kind> <Leff> <B> <G> <Gv> <S> <O> <edge1> <edge2> <nsl> <design_mev> <sr> <x px y py t delta> [lf1 lf2]
//       Прогнать ОДНО тело элемента над состоянием и напечатать результат (для регрессии против Python,
//       см. test_matrix_core_cpp.py). kind: dipole|dq|quad|sext|oct|kicker. B,G,Gv,S,O — физкоэф. mat_*
//       (абсолютные у любого вида; мультиполи допустимы у любого вида — тело выбирается по полю).
//       lf1/lf2 (необязательно) — L_f торцов, фриндж HGAP/FINT (верт. линза грани tg(E−ψ)/ρ).
//   <exp_id>   (в разработке) — автономный трек всего кольца из БД. Каркас/чтение БД — следующий слой;
//              физический движок (matrix_core.h) уже готов и провалидирован.
#include <iostream>
#include <cctype>
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
    // необязательный хвост [lf1 lf2] — длины зоны спада торцов L_f (фриндж HGAP/FINT), м
    double lf1 = (argc >= 22) ? std::stod(argv[20]) : 0.0, lf2 = (argc >= 22) ? std::stod(argv[21]) : 0.0;
    // nsl_box = nsl: «ящик» корректора режет на те же N, что box-тело Python (n_slices строки) — паритет теста
    mc::Element el = mc::Element::make(kind, Leff, B, G, Gv, S, O, e1, e2, nsl, design, sr, 4, false, nsl, lf1, lf2);
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

// Резонатор опыта из строки radiation_params.csv (rf: none/ideal/sine + rf_voltage_mv/rf_harmonic/rf_phase_deg/
// rf_freq_hz) — как matrix_tracker.load_radiation + rf_spec. Пустая фаза или auto ⇒ синхронная.
struct RfRow { std::string type = "none"; double volt = 0, h = 0, phase_deg = 0, freq = 0, ideal_energy = 0;
               bool phase_auto = true; };

static bool readRfRow(const mc::Table& rp, const std::vector<std::string>& r, RfRow& o, std::string& err) {
    auto c = [&](const char* k) { return mc::cell(r, rp.col(k)); };
    o.type = c("rf"); if (o.type.empty()) o.type = "none";
    for (auto& ch : o.type) ch = (char)std::tolower((unsigned char)ch);
    if (o.type != "none" && o.type != "ideal" && o.type != "sine") {
        err = "неизвестный тип резонатора rf=" + o.type + " (допустимо: none, ideal, sine)"; return false;
    }
    if (o.type == "ideal" && !c("rf_ideal_energy_mev").empty())      // энергия после идеального резонатора, МэВ
        o.ideal_energy = std::stod(c("rf_ideal_energy_mev"));
    if (o.type == "sine") {
        if (c("rf_voltage_mv").empty() || c("rf_harmonic").empty()) {
            err = "тип rf=sine требует rf_voltage_mv и rf_harmonic"; return false;
        }
        o.volt = std::stod(c("rf_voltage_mv")); o.h = std::stod(c("rf_harmonic"));
        std::string ph = c("rf_phase_deg");
        for (auto& ch : ph) ch = (char)std::tolower((unsigned char)ch);
        o.phase_auto = ph.empty() || ph == "auto";
        if (!o.phase_auto) o.phase_deg = std::stod(ph);
        if (!c("rf_freq_hz").empty()) o.freq = std::stod(c("rf_freq_hz"));
    }
    return true;
}

static mc::RfSpec makeRf(const RfRow& rr, const mc::Injection& inj, double design) {
    mc::RfSpec rf;
    rf.L = inj.ring_len; rf.end_share = inj.end_share;
    if (rr.type == "ideal") {                    // δ после резонатора = E/E_design − 1 (пусто — проектная, δ = 0)
        rf.type = mc::RfSpec::IDEAL; rf.design = design;
        rf.ideal_delta = (rr.ideal_energy > 0) ? rr.ideal_energy / design - 1.0 : 0.0;
    }
    if (rr.type == "sine") {
        rf.type = mc::RfSpec::SINE; rf.volt = rr.volt; rf.h = rr.h; rf.design = design;
        rf.prepare(inj, rr.freq, rr.phase_auto, rr.phase_deg);
    }
    return rf;
}

// Конфиг опыта (решётка/группа/плоскости/Settings/радиация) — общее для --track и --full. Пусто ⇒ err.
struct ExpCfg { int lattice = -1, group = -1; mc::Settings stg; mc::ExpPlanes planes;
                bool sr = false; double design = 6000.0; RfRow rr; };
static bool readExpCfg(const std::string& ws, int exp, ExpCfg& c, std::string& err) {
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
        try { c.planes = mc::readExpPlanes(ex, r); }
        catch (const std::exception& e) { err = e.what(); return false; }
        break;
    }
    if (c.lattice < 0) { err = "Experiment " + std::to_string(exp) + " not found"; return false; }
    // радиация: строка в radiation_params.csv → СИ вкл (энергия/ВЧ оттуда); нет строки → выкл.
    mc::Table rp = mc::readTable(ws + "tables/radiation_params.csv");
    int re = rp.col("exp_id"), rd = rp.col("design_energy_mev");
    for (auto& r : rp.rows) if (mc::cell(r, re) == std::to_string(exp)) {
        c.sr = true;
        if (!mc::cell(r, rd).empty()) c.design = std::stod(mc::cell(r, rd));
        if (!readRfRow(rp, r, c.rr, err)) { err = "radiation_params: " + err; return false; }
        break;
    }
    return true;
}

static void printPath(const mc::Injection& inj) {
    std::printf("# start plane: p=(%.9g, %.9g, %.9g) n=(%.9g, %.9g, %.9g); %s; line_len=%.17g m\n",
                inj.plane.center.x, inj.plane.center.y, inj.plane.center.z,
                inj.plane.axis.x, inj.plane.axis.y, inj.plane.axis.z,
                inj.closed ? "closed ring (finish = start)" : "open line (portal finish -> start)", inj.ring_len);
    if (!inj.closed)
        std::printf("# finish plane: p=(%.9g, %.9g, %.9g) n=(%.9g, %.9g, %.9g)\n",
                    inj.finish.center.x, inj.finish.center.y, inj.finish.center.z,
                    inj.finish.axis.x, inj.finish.axis.y, inj.finish.axis.z);
    for (auto& nd : inj.ordered)
        if (nd.cavity) std::printf("# rf cavity %s: s=%.9g m from start, share=%.9g\n", nd.dev.c_str(), nd.s_c, nd.cav_share);
    if (inj.end_share > 0) std::printf("# rf cavity at the end of the pass: share=%.9g\n", inj.end_share);
    for (auto& id : inj.skipped) std::printf("# rf cavity %s: off the start->finish path, skipped\n", id.c_str());
    auto pk = [&](int i, const char* where) {
        const mc::KickerPulse& k = inj.kickers[i];
        std::printf("# pulsed kicker %s%s: s=%.9g m from start, theta=(%.9g, %.9g) rad, %s D=%.9g s rise=%.9g s, "
                    "delay=%.12g s, period=%.9g s\n", k.id.c_str(), where, k.s_k, k.thx, k.thy,
                    k.shape == 0 ? "half_sine" : "trapezoid", k.D, k.rise, k.delay, k.period > 0 ? k.period : 0.0);
    };
    for (int i : inj.kick_begin) pk(i, " (start plane, beginning of the pass)");
    for (auto& nd : inj.ordered) if (nd.kicker) pk(nd.kick_idx, "");
    for (int i : inj.kick_end) pk(i, " (finish plane, end of the pass)");
    for (auto& id : inj.kick_skipped) std::printf("# pulsed kicker %s: off the start->finish path, skipped\n", id.c_str());
}

// --track <exp_id> <turns> [x0_mm y0_mm] — АВТОНОМНЫЙ трек опыта из БД: пролёты старт → финиш (у разомкнутой
// линии — портал) + ВЧ. Пишет results/trajectories_local_exp<exp>.npz (mode=long) для traj_analyzer + печатает
// финальные состояния (регрессия): "<idx> x px y py t delta surv" или "<idx> LOST surv". x0/y0 (мм) ⇒ ОДНА частица.
static int runTrack(const std::string& ws, int exp, int turns, double x0_mm, double y0_mm) {
    ExpCfg c; std::string err;
    if (!readExpCfg(ws, exp, c, err)) { std::cerr << err << "\n"; return 1; }
    mc::Injection inj; mc::RfSpec rf;
    try {
        inj = mc::inject(ws, c.lattice, c.group, c.planes, c.design, c.sr, c.stg, c.rr.type != "none", exp);
        rf = makeRf(c.rr, inj, c.design);
    } catch (const std::exception& e) { std::cerr << e.what() << "\n"; return 1; }
    if (x0_mm != 0.0 || y0_mm != 0.0) {                        // ОДНА частица со смещением (скан апертуры/затухания)
        inj.states.clear();
        inj.states.push_back({ x0_mm * 1e-3, 0.0, y0_mm * 1e-3, 0.0, 0.0, 0.0 });
    }
    std::printf("# exp %d lattice %d group %d: %zu particles, %s, %zu elements, sr=%d rf=%s design=%.0f turns=%d\n",
                exp, c.lattice, c.group, inj.states.size(), inj.closed ? "closed" : "open",
                inj.ordered.size() - inj.n_cav - inj.n_kick, (int)c.sr, rf.name(), c.design, turns);
    printPath(inj);
    if (rf.type == mc::RfSpec::IDEAL)
        std::printf("# rf ideal: delta after cavity = %.17g (E = %.9g MeV)\n", rf.ideal_delta,
                    c.design * (1.0 + rf.ideal_delta));
    if (rf.type == mc::RfSpec::SINE)
        std::printf("# rf sine: V=%.6g MV h=%.6g C_RF=%.17g m phase=%.17g rad%s U0=%.17g MeV\n",
                    rf.volt, rf.h, rf.c_rf, rf.phase, rf.auto_phase ? " (auto)" : "", rf.U0);
    // трек + сбор turn-by-turn для .npz (mode=long). Частицы НЕЗАВИСИМЫ → OpenMP (кольцо read-only).
    std::vector<double> tx, txp, ty, typ;
    std::vector<int64_t> tt_counts, surv;
    std::vector<mc::LongResult> results(inj.states.size());
#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < (int)inj.states.size(); ++i)
        results[i] = mc::trackLong(inj, inj.states[i], turns, rf);
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
    npz.addScalarI8("lattice", c.lattice); npz.addScalarI8("group", c.group);
    npz.addScalarI8("sr", c.sr ? 1 : 0); npz.addScalarI8("isr", 0);
    npz.addScalarStr("rf", rf.name());
    npz.addScalarF8("design_mev", c.design);
    std::error_code ec; std::filesystem::create_directories(ws + "results", ec);
    std::string npzPath = ws + "results/trajectories_local_exp" + std::to_string(exp) + ".npz";
    if (npz.write(npzPath)) std::printf("# npz: %s\n", npzPath.c_str());
    else std::fprintf(stderr, "npz write failed: %s\n", npzPath.c_str());
    return 0;
}

// --full <exp> <turns> [x0_mm y0_mm] — ПОЛНЫЙ трек (x(S)/спектр/тюны): запись на границах элементов за все
// обороты → results/trajectories_local_exp<exp>.npz (mode=full). БЕЗ ВЧ (как local_trajectories).
static int runFull(const std::string& ws, int exp, int turns, double x0_mm, double y0_mm) {
    ExpCfg c; std::string err;
    if (!readExpCfg(ws, exp, c, err)) { std::cerr << err << "\n"; return 1; }
    const int TARGET = 30000;
    mc::Injection inj;
    try { inj = mc::inject(ws, c.lattice, c.group, c.planes, c.design, c.sr, c.stg, false, exp); }
    catch (const std::exception& e) { std::cerr << e.what() << "\n"; return 1; }
    if (x0_mm != 0.0 || y0_mm != 0.0) {
        inj.states.clear(); inj.states.push_back({ x0_mm * 1e-3, 0.0, y0_mm * 1e-3, 0.0, 0.0, 0.0 });
    }
    std::vector<std::string> elem_names; elem_names.push_back("INJ_PLANE");
    for (auto& nd : inj.ordered) elem_names.push_back(nd.dev.empty() ? "?" : nd.dev);
    if (!inj.closed) elem_names.push_back("FINISH_PLANE");
    std::printf("# exp %d lattice %d group %d: %zu particles, %s, %zu elements, sr=%d design=%.0f turns=%d (mode=full)\n",
                exp, c.lattice, c.group, inj.states.size(), inj.closed ? "closed" : "open",
                inj.ordered.size(), (int)c.sr, c.design, turns);
    printPath(inj);
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
                 "  --element <kind Leff B G Gv S O edge1 edge2 nsl design_mev sr> <x px y py t delta> [lf1 lf2]\n"
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
