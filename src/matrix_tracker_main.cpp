// ==========================================================================
// matrix_tracker_main.cpp — CLI для C++ матричного трекера (MatrixTracker.h).
// Порт matrix_tracker.py: собирает кольцо из тех же CSV, трекает тест-пучок,
// печатает финальные состояния (для сверки с Python) и замеряет время оборота.
//
// Запуск:  MatrixTracker.exe [lattice=107] [source=physical|current] [reps=200]
// ==========================================================================
#include <iostream>
#include <vector>
#include <string>
#include <chrono>
#include <iomanip>
#ifdef _WIN32
#include <windows.h>
#endif
#ifdef _OPENMP
#include <omp.h>
#endif
#include "MatrixTracker.h"
#include "Workspace.h"

using namespace mtx;

int main(int argc, char* argv[]) {
#ifdef _WIN32
    SetConsoleOutputCP(CP_UTF8);
#endif
    int lattice     = (argc > 1) ? std::stoi(argv[1]) : 107;
    std::string src = (argc > 2) ? std::string(argv[2]) : "physical";
    int reps        = (argc > 3) ? std::stoi(argv[3]) : 200;

    std::string ws = findWorkspace(argv[0]);
    auto tc = [](std::chrono::high_resolution_clock::time_point a) {
        return std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - a).count();
    };

    auto t0 = std::chrono::high_resolution_clock::now();
    DB db = loadDB(ws);
    auto mags = buildLine(ws, db, lattice, -1, -1, src);
    double tBuild = tc(t0);

    std::cout << "Решётка " << lattice << ", источник='" << src << "', q/p0=" << std::fixed
              << std::setprecision(4) << QP << ". Магнитов: " << mags.size()
              << " (сборка " << std::setprecision(3) << tBuild * 1e3 << " мс)\n";
    if (mags.empty()) { std::cout << "Пусто — проверь lattice_id/connectivity.\n"; return 1; }

    std::cout << "Первые магниты (dev | L,мм | yaw,мрад | theta0,мрад | k1 | k2 | k3):\n";
    for (size_t i = 0; i < std::min<size_t>(8, mags.size()); ++i) {
        const auto& m = mags[i];
        std::string devs; for (size_t j = 0; j < m.devs.size(); ++j) devs += (j ? "+" : "") + m.devs[j];
        std::cout << "  " << std::left << std::setw(26) << devs << std::right << std::setprecision(2)
                  << " L=" << std::setw(6) << m.L * 1e3 << " yaw=" << std::setw(8) << m.yaw * 1e3
                  << " th0=" << std::setw(8) << m.theta0 * 1e3
                  << " k1=" << std::setw(8) << std::setprecision(3) << m.k1
                  << " k2=" << std::setw(9) << m.k2 << " k3=" << std::setw(9) << m.k3 << "\n";
    }

    // тест-пучок (детерминированный — для сверки с Python)
    std::vector<State> beam = {
        {0, 0, 0, 0, 0, 0.0},        // опорная
        {1e-3, 0, 0, 0, 0, 0.0},     // +1мм x
        {0, 0, 1e-3, 0, 0, 0.0},     // +1мм y
        {2e-3, 0, 0, 0, 0, 0.0},     // +2мм x
        {0, 0, 0, 0, 0, 0.005},      // δ=+0.5%
    };
    const char* lbl[] = {"опорная  ", "x0=1мм   ", "y0=1мм   ", "x0=2мм   ", "δ=+0.5%  "};

    std::cout << "\nФинал (x,pₓ,y,p_y мм/мрад, t м, δ):\n" << std::setprecision(6);
    mtx::BPMAcc monitors;                                    // датчики положения пучка (BPM), если есть в решётке
    for (size_t i = 0; i < beam.size(); ++i) {
        TrackResult r = track(beam[i], mags, &monitors);
        if (r.lostAt >= 0)
            std::cout << "  " << lbl[i] << ": ПОТЕРЯНА на магните #" << r.lostAt << " (" << r.reason << ")\n";
        else
            std::cout << "  " << lbl[i] << ": x=" << std::setw(11) << r.fin.x * 1e3
                      << " pₓ=" << std::setw(11) << r.fin.px * 1e3
                      << " y=" << std::setw(11) << r.fin.y * 1e3
                      << " p_y=" << std::setw(11) << r.fin.py * 1e3
                      << " t=" << std::setw(11) << r.fin.t
                      << " δ=" << r.fin.d << "\n";
    }
    if (!monitors.empty()) {                                 // буфер в ОЗУ → несовершенства → СБРОС в БД после прогона
        mtx::BPMRegistry reg;
        reg.config = mtx::loadBPMConfig(ws + "tables/bpm_config.csv", lattice);
        reg.addTurn(monitors, 0);                            // 1 оборот; id_exp = номер решётки
        reg.flush(std::to_string(lattice), ws + "tables/bpm_readings.csv");
        std::cout << "\nДатчики BPM: " << reg.records.size() << " показаний → bpm_readings.csv (id_exp=" << lattice << ")\n";
        for (const auto& r : reg.records)
            std::cout << "  датчик " << r.sid << ": x=" << r.x * 1e3 << "мм y=" << r.y * 1e3
                      << "мм σx=" << r.sigx * 1e3 << "мм n=" << r.n << "\n";
    }

    // ---- замер времени: reps оборотов всего пучка (последоват. и OpenMP) ----
    volatile double sink = 0;
    auto tS = std::chrono::high_resolution_clock::now();
    for (int rep = 0; rep < reps; ++rep)
        for (auto& b : beam) { TrackResult r = track(b, mags); sink += r.fin.x; }
    double tSeq = tc(tS) / reps;

    double tPar = tSeq;
#ifdef _OPENMP
    auto tP = std::chrono::high_resolution_clock::now();
    for (int rep = 0; rep < reps; ++rep) {
        #pragma omp parallel for reduction(+:sink) schedule(static)
        for (int i = 0; i < (int)beam.size(); ++i) { TrackResult r = track(beam[i], mags); sink += r.fin.x; }
    }
    tPar = tc(tP) / reps;
#endif
    (void)sink;
    std::cout << "\nВремя оборота (" << mags.size() << " магнитов, " << beam.size() << " частиц):\n";
    std::cout << std::setprecision(3) << "  последовательно: " << tSeq * 1e3 << " мс ("
              << tSeq / beam.size() * 1e3 << " мс/частицу)\n";
#ifdef _OPENMP
    std::cout << "  OpenMP (" << omp_get_max_threads() << " потоков): " << tPar * 1e3 << " мс\n";
#endif
    return 0;
}
