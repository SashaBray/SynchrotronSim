#ifndef FINISH_PLANE_H
#define FINISH_PLANE_H

#include <cmath>
#include "FieldMap.h"   // FieldVector

// Финишная плоскость, заданная тремя точками p1, p2, p3.
//
// Семантика срабатывания:
//   В hasCrossed() сначала строится знаковое расстояние d частицы от плоскости
//   (нормировка такая, что startPos имеет d < 0).  Срабатывание -- ИЗМЕНЕНИЕ
//   ЗНАКА d между последовательными вызовами.  Это позволяет корректно
//   обрабатывать случаи, когда орбита пересекает одну и ту же плоскость
//   дважды за оборот (например, плоскость за SP20: в X-Y след плоскости
//   пересекает кольцевую орбиту и около SP1, и около SP20).
//
//   skipSteps пропускает первые N вызовов hasCrossed() ДО любой проверки --
//   это удобный механизм пройти "первое" пересечение и поймать второе.
//
//   Каждая частица должна иметь СВОЮ копию FinishPlane (callCount/lastD --
//   per-particle state; race-free под OpenMP).
struct FinishPlane {
    FieldVector p0;
    FieldVector n;
    int skipSteps   = 0;
    int callCount   = 0;
    double lastD    = 0.0;
    bool initialized = false;

    void setup(const FieldVector& a, const FieldVector& b, const FieldVector& c,
               const FieldVector& startPos, int skip = 0) {
        p0 = a;
        double ux = b.x - a.x, uy = b.y - a.y, uz = b.z - a.z;
        double vx = c.x - a.x, vy = c.y - a.y, vz = c.z - a.z;
        double nx = uy * vz - uz * vy;
        double ny = uz * vx - ux * vz;
        double nz = ux * vy - uy * vx;
        double m  = std::sqrt(nx * nx + ny * ny + nz * nz);
        if (m > 1e-30) { nx /= m; ny /= m; nz /= m; }
        else           { nx = 1.0; ny = 0.0; nz = 0.0; }
        // Нормаль выбирается так, чтобы startPos лежал на отрицательной стороне.
        double sd = (startPos.x - a.x) * nx + (startPos.y - a.y) * ny + (startPos.z - a.z) * nz;
        if (sd > 0.0) { nx = -nx; ny = -ny; nz = -nz; }
        n = {nx, ny, nz};
        skipSteps    = skip;
        callCount    = 0;
        initialized  = false;
        lastD        = 0.0;
    }

    bool hasCrossed(const FieldVector& pos) {
        if (callCount++ < skipSteps) return false;
        double d = (pos.x - p0.x) * n.x + (pos.y - p0.y) * n.y + (pos.z - p0.z) * n.z;
        if (!initialized) {
            initialized = true;
            lastD = d;
            return false;
        }
        bool crossed = (lastD < 0.0 && d >= 0.0) || (lastD >= 0.0 && d < 0.0);
        lastD = d;
        return crossed;
    }

    void reset() { callCount = 0; initialized = false; lastD = 0.0; }
};

#endif // FINISH_PLANE_H
