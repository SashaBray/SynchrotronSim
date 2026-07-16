// ==========================================================================
// matrix_optimizer_main.cpp — C++ оптимизатор матричной модели (интеграция
// алгоритма настройки в ядро). Читает кольцо (MatrixTracker.h) + knob-map
// (matrix_knobs.csv, генерит build_knobmap.py). Цель — ВЫЖИВАНИЕ пучка на всём
// кольце: минимизируем макс. смещение от оси + штраф за потери. Метод —
// координатный спуск (line-search) по ручкам. Ручки: общие standard + уникальные
// junction; изгибы/корректоры заморожены (нет ручек). Результат — в
// matrix_tune_result.csv (knob_id;value), применяется apply_knobmap.py.
//
// Запуск: MatrixOptimizer.exe [lattice=101] [rounds=20] [aper=0.03]
// ==========================================================================
#include <iostream>
#include <vector>
#include <string>
#include <chrono>
#include <iomanip>
#include <cmath>
#ifdef _WIN32
#include <windows.h>
#endif
#ifdef _OPENMP
#include <omp.h>
#endif
#include "MatrixTracker.h"
#include "Workspace.h"

using namespace mtx;

// ---- ручка: тип коэф. (G/S/O), значение, индексы управляемых магнитов --------
struct Knob { char coeff; double val; std::vector<int> idx; std::string sp_type; };

std::vector<Knob> readKnobs(const std::string& ws) {
    std::vector<Knob> ks;
    // метаданные ручек: knob_id;coeff;init;sp_type;n_groups
    for (auto& r : readCSV(ws + "tables/matrix_knobs.csv")) {
        if (r.size() < 4) continue;
        Knob k; k.coeff = r[1][0]; k.val = std::stod(r[2]); k.sp_type = r[3];
        ks.push_back(k);
    }
    // индексы управляемых магнитов — из ДОЧЕРНЕЙ таблицы matrix_knob_members (knob_id;magnet_index)
    for (auto& r : readCSV(ws + "tables/matrix_knob_members.csv")) {
        if (r.size() < 2) continue;
        int kid = std::stoi(r[0]);
        if (kid >= 0 && kid < (int)ks.size()) ks[kid].idx.push_back(std::stoi(r[1]));
    }
    return ks;
}

// ---- рабочий пучок ±1мм: опорная + полное кольцо радиуса 1мм в (x,y) + дисперсия = 11 частиц ----
// (>2мм считаем неустойчивыми — оптимизируем только реальную рабочую область ±1мм)
std::vector<State> makeBeam(double r) {
    std::vector<State> P{ {0,0,0,0,0,0},
        { r,0,0,0,0,0}, {-r,0,0,0,0,0}, {0,0, r,0,0,0}, {0,0,-r,0,0,0},
        { r,0, r,0,0,0}, {-r,0, r,0,0,0}, { r,0,-r,0,0,0}, {-r,0,-r,0,0,0},
        {0,0,0,0,0, 0.005}, {0,0,0,0,0,-0.005} };
    return P;
}

// ---- применить значения ручек к магнитам (k1/k2/k3 = QP·знач; theta0 не трогаем) ----
void applyKnobs(std::vector<Magnet>& mags, const std::vector<Knob>& ks, const std::vector<double>& v) {
    for (size_t i = 0; i < ks.size(); ++i) {
        double kv = QP * v[i]; char c = ks[i].coeff;
        for (int gi : ks[i].idx) {
            if (c == 'G') mags[gi].k1 = kv; else if (c == 'S') mags[gi].k2 = kv; else if (c == 'O') mags[gi].k3 = kv;
        }
    }
    for (auto& m : mags)
        m.nonzero = (std::abs(m.theta0) > 0 || std::abs(m.k1) > 0 || std::abs(m.k2) > 0 || std::abs(m.k3) > 0);
}

// ---- трек частицы с апертурой: макс. смещение выжившей / магнит потери -------
inline void trackLoss(State s, const std::vector<Magnet>& mags, double aper, double& maxoff, int& lostAt) {
    maxoff = 0.0; lostAt = -1;
    if (!body(mags[0], s)) { maxoff = aper; lostAt = 0; return; }
    for (size_t j = 1; j < mags.size(); ++j) {
        if (!modifiedDrift(s, mags[j - 1], mags[j])) { maxoff = aper; lostAt = (int)j; return; }
        double off = std::sqrt(s.x * s.x + s.y * s.y);
        if (off > maxoff) maxoff = off;
        if (off > aper) { maxoff = aper; lostAt = (int)j; return; }
        if (!body(mags[j], s)) { maxoff = aper; lostAt = (int)j; return; }
    }
}

// ---- loss: среднее по частицам (выжившие → макс. смещение; потерянные → aper·(2−доля), раньше=хуже) ----
double lossOf(const std::vector<Magnet>& mags, const std::vector<State>& beam, double aper,
              int& horizon, int& nlost) {
    int total = (int)mags.size(); int np = (int)beam.size();
    double sum = 0; horizon = total; nlost = 0;
    std::vector<double> cost(np); std::vector<int> lost(np);
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < np; ++i) {
        double mo; int la; trackLoss(beam[i], mags, aper, mo, la);
        cost[i] = (la < 0) ? mo : aper * (2.0 - (double)la / total);
        lost[i] = la;
    }
    for (int i = 0; i < np; ++i) {
        sum += cost[i];
        if (lost[i] >= 0) { nlost++; if (lost[i] < horizon) horizon = lost[i]; }
    }
    return sum / np;
}

int main(int argc, char* argv[]) {
#ifdef _WIN32
    SetConsoleOutputCP(CP_UTF8);
#endif
    int lattice = (argc > 1) ? std::stoi(argv[1]) : 101;
    int rounds  = (argc > 2) ? std::stoi(argv[2]) : 20;
    double aper = (argc > 3) ? std::stod(argv[3]) : 0.03;
    double beamR = (argc > 4) ? std::stod(argv[4]) * 1e-3 : 1e-3;   // радиус пучка, мм → м
    std::string ws = findWorkspace(argv[0]);
    auto now = [] { return std::chrono::high_resolution_clock::now(); };
    auto sec = [](auto a, auto b) { return std::chrono::duration<double>(b - a).count(); };

    DB db = loadDB(ws);
    auto mags = buildLine(ws, db, lattice, -1, -1, "physical");
    auto knobs = readKnobs(ws);
    auto beam = makeBeam(beamR);
    if (mags.empty() || knobs.empty()) { std::cout << "Нет магнитов/ручек — проверь БД/knob-map.\n"; return 1; }
    int total = (int)mags.size();

    std::vector<double> val(knobs.size()), step(knobs.size());
    for (size_t i = 0; i < knobs.size(); ++i) { val[i] = knobs[i].val; step[i] = std::abs(val[i]) * 0.25 + 3e-4; }

    applyKnobs(mags, knobs, val);
    int horizon, nlost; double best = lossOf(mags, beam, aper, horizon, nlost);
    std::cout << "Кольцо " << lattice << ": " << total << " магнитов, " << knobs.size()
              << " ручек, пучок " << beam.size() << " частиц (радиус " << beamR * 1e3
              << "мм), апертура " << aper * 1e3 << "мм.\n";
    std::cout << std::fixed << std::setprecision(2)
              << "Старт (warm-start .lte): loss=" << best * 1e3 << "мм  горизонт=" << horizon << "/" << total
              << "  потеряно=" << nlost << "\n";

    auto t0 = now();
    for (int rnd = 0; rnd < rounds; ++rnd) {
        int improved = 0;
        for (size_t i = 0; i < knobs.size(); ++i) {
            for (double s : { step[i], -step[i] }) {
                double old = val[i]; val[i] = old + s;
                applyKnobs(mags, knobs, val);
                int h, n; double l = lossOf(mags, beam, aper, h, n);
                if (l < best - 1e-12) { best = l; horizon = h; nlost = n; improved++; break; }
                val[i] = old;
            }
        }
        if (improved == 0) for (auto& s : step) s *= 0.5;
        applyKnobs(mags, knobs, val);   // зафиксировать лучшее
        std::cout << "  раунд " << std::setw(2) << rnd + 1 << ": loss=" << std::setw(7) << best * 1e3
                  << "мм  горизонт=" << std::setw(4) << horizon << "/" << total << "  потеряно=" << std::setw(2)
                  << nlost << "  улучшено=" << std::setw(3) << improved << "  [" << std::setprecision(0)
                  << sec(t0, now()) << "с]" << std::setprecision(2) << "\n";
        if (improved == 0) { bool tiny = true; for (double s : step) if (s > 1e-7) tiny = false; if (tiny) break; }
    }
    std::cout << "Финал: loss=" << best * 1e3 << "мм  горизонт=" << horizon << "/" << total
              << "  потеряно=" << nlost << "  за " << std::setprecision(0) << sec(t0, now()) << "с\n";

    // запись настроенных ручек
    std::ofstream out(ws + "tables/matrix_tune_result.csv");
    out << "knob_id;value\n";
    for (size_t i = 0; i < val.size(); ++i) out << i << ";" << std::setprecision(12) << val[i] << "\n";
    std::cout << "Записано: workspace/tables/matrix_tune_result.csv (" << val.size()
              << " ручек). Применить: py -3.14 apply_knobmap.py\n";
    return 0;
}
