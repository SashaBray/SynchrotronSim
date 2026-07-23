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
#include "matrix_core.h"

#ifdef _WIN32
#include <windows.h>
#endif

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

int main(int argc, char* argv[]) {
    #ifdef _WIN32
    SetConsoleOutputCP(CP_UTF8);
    #endif
    if (argc >= 2 && std::string(argv[1]) == "--element") {
        return runElement(argc, argv);
    }
    std::cout << "MatrixTracker (SynchrotronSim matrix core, C++)\n"
                 "  --element <kind Leff B G Gv S O edge1 edge2 nsl design_mev sr> <x px y py t delta>\n"
                 "      прогнать тело одного элемента (валидация против Python-эталона)\n"
                 "  <exp_id>   автономный трек кольца из БД — в разработке (движок matrix_core.h готов)\n";
    return (argc < 2) ? 1 : 0;
}
