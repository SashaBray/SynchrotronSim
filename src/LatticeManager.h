#ifndef LATTICE_MANAGER_H
#define LATTICE_MANAGER_H

#include "FieldMap.h"
#include <vector>
#include <memory>
#include <cmath>
#include <algorithm>

enum class FieldType { MAGNETIC, ELECTRIC };

struct TransformMatrix {
    double m[3][3];
    double tx, ty, tz;

    FieldVector transformPointInv(const FieldVector& p) const {
        double dx = p.x - tx;
        double dy = p.y - ty;
        double dz = p.z - tz;
        return {
            dx * m[0][0] + dy * m[1][0] + dz * m[2][0],
            dx * m[0][1] + dy * m[1][1] + dz * m[2][1],
            dx * m[0][2] + dy * m[1][2] + dz * m[2][2]
        };
    }

    FieldVector rotateVector(const FieldVector& v) const {
        return {
            v.x * m[0][0] + v.y * m[0][1] + v.z * m[0][2],
            v.x * m[1][0] + v.y * m[1][1] + v.z * m[1][2],
            v.x * m[2][0] + v.y * m[2][1] + v.z * m[2][2]
        };
    }
};

class MagneticElement {
public:
    std::shared_ptr<FieldMap> map;
    FieldType type;
    TransformMatrix transform;
    double currentScale;
    bool useSymmetry;

    // Новые калибровочные параметры
    bool isActive;         // Из колонки "Active" (True/False)
    double fieldDirection; // Из колонки "field_direction" (1/-1)
    double correctCoeff;   // Из колонки "correct_coeff"

    // Bounding sphere в МИРОВЫХ координатах -- для быстрого прескрининга
    // в getTotalFields (избегаем transformPointInv для далёких магнитов).
    // Считается через computeBounds() после задания transform/map/useSymmetry.
    FieldVector worldCenter = {0,0,0};
    double      worldRadius = 0.0;
    double      worldRadius2 = 0.0;  // = worldRadius² (для сравнения без sqrt)

    void setOrientation(double x, double y, double z, double yaw, double pitch, double roll) {
        transform.tx = x; transform.ty = y; transform.tz = z;
        double sy = sin(yaw),   cy = cos(yaw);
        double sp = sin(pitch), cp = cos(pitch);
        double sr = sin(roll),  cr = cos(roll);

        transform.m[0][0] = cy * cp;
        transform.m[0][1] = cy * sp * sr - sy * cr;
        transform.m[0][2] = cy * sp * cr + sy * sr;
        transform.m[1][0] = sy * cp;
        transform.m[1][1] = sy * sp * sr + cy * cr;
        transform.m[1][2] = sy * sp * cr - cy * sr;
        transform.m[2][0] = -sp;
        transform.m[2][1] = cp * sr;
        transform.m[2][2] = cp * cr;
    }

    // Должно быть вызвано ПОСЛЕ установки map, useSymmetry, transform.
    // Считает мировой центр (трансформ от центра локального bbox) и радиус,
    // покрывающий ВСЕ возможные точки локального bbox после поворота.
    void computeBounds() {
        if (!map) { worldCenter = {transform.tx, transform.ty, transform.tz};
                    worldRadius = 0.0; worldRadius2 = 0.0; return; }
        const GridHeader& h = map->headerRef();
        // В локальных координатах bbox карты:
        // - При useSymmetry поле зеркалится по X: эффективный диапазон [-maxX, +maxX].
        // - Иначе используем [minX, maxX] как есть.
        double lxMin = h.minX, lxMax = h.maxX;
        if (useSymmetry) { lxMin = -h.maxX; lxMax = h.maxX; }
        double lyMin = h.minY, lyMax = h.maxY;
        double lzMin = h.minZ, lzMax = h.maxZ;
        // Локальный центр bbox
        double lcx = 0.5 * (lxMin + lxMax);
        double lcy = 0.5 * (lyMin + lyMax);
        double lcz = 0.5 * (lzMin + lzMax);
        // Полупролёты bbox
        double hx = 0.5 * (lxMax - lxMin);
        double hy = 0.5 * (lyMax - lyMin);
        double hz = 0.5 * (lzMax - lzMin);
        // Мировой центр = transform applied to (lcx, lcy, lcz)
        // (forward: world = R * local + t). Используем поэлементно.
        worldCenter.x = transform.m[0][0]*lcx + transform.m[0][1]*lcy + transform.m[0][2]*lcz + transform.tx;
        worldCenter.y = transform.m[1][0]*lcx + transform.m[1][1]*lcy + transform.m[1][2]*lcz + transform.ty;
        worldCenter.z = transform.m[2][0]*lcx + transform.m[2][1]*lcy + transform.m[2][2]*lcz + transform.tz;
        // Радиус = sqrt(hx²+hy²+hz²) -- покрывает любой поворот bbox.
        worldRadius  = std::sqrt(hx*hx + hy*hy + hz*hz);
        worldRadius2 = hx*hx + hy*hy + hz*hz;
    }

    // Быстрый прескрининг: вернуть true, если pos СОВСЕМ ВНЕ bounding sphere.
    // Возвращаем bool вместо early-return в addContribution, потому что в
    // некоторых callsite'ах удобнее вызвать отдельно.
    inline bool isOutsideBounds(const FieldVector& pos) const {
        double dx = pos.x - worldCenter.x;
        double dy = pos.y - worldCenter.y;
        double dz = pos.z - worldCenter.z;
        return (dx*dx + dy*dy + dz*dz) > worldRadius2;
    }

    void addContribution(const FieldVector& globalPos, FieldVector& totalB, FieldVector& totalE) const {
        // Прескрининг по bounding sphere -- 7 flops + 1 compare вместо
        // 12+ flops на полный transformPointInv + interp setup.
        if (isOutsideBounds(globalPos)) return;

        FieldVector localPos = transform.transformPointInv(globalPos);
        FieldVector localF = map->getInterpolatedField(localPos, useSymmetry);

        if (localF.x == 0 && localF.y == 0 && localF.z == 0) return;

        FieldVector globalF = transform.rotateVector(localF);

        // Применяем новые коэффициенты калибровки
        double finalScale = correctCoeff * fieldDirection;

        // Масштабирование по току применяется ТОЛЬКО для активных электромагнитов
        if (isActive) {
            finalScale *= currentScale;
        }

        globalF.x *= finalScale;
        globalF.y *= finalScale;
        globalF.z *= finalScale;

        if (type == FieldType::MAGNETIC) {
            totalB.x += globalF.x; totalB.y += globalF.y; totalB.z += globalF.z;
        } else {
            totalE.x += globalF.x; totalE.y += globalF.y; totalE.z += globalF.z;
        }
    }
};

// ==========================================================================
// AcceleratorConfig: контейнер элементов + 2D spatial bucket grid по (x, y).
//
// Идея ускорения: вместо линейного цикла по 2762 элементам, на каждом шаге
// частицы дёргаем только элементы из bucket'а, содержащего текущую позицию.
// Bucket size выбирается так, чтобы покрывать максимальный worldRadius
// магнита => элемент регистрируется максимум в 2x2=4 bucket'ах.
//
// Сетка строится по X-Y (Z ≈ 0 на кольце), bucket'ы -- разреженный 2D-массив.
// ==========================================================================
class AcceleratorConfig {
private:
    std::vector<MagneticElement> elements;
    double maxRadiusSq = 0;

    // Spatial grid
    bool   gridBuilt = false;
    double bucketSize = 1.0;   // 1 метр -- покрывает любой реальный магнит
    int    bxMin = 0, byMin = 0;
    int    bxN = 0,  byN = 0;
    // bucketGrid[bx * byN + by] = список индексов элементов
    std::vector<std::vector<int>> bucketGrid;

    static inline int floorDiv(double a, double b) {
        return (int)std::floor(a / b);
    }

public:
    void addElement(const MagneticElement& el) {
        elements.push_back(el);
        elements.back().computeBounds();
        double r_sq = el.transform.tx * el.transform.tx + el.transform.ty * el.transform.ty + el.transform.tz * el.transform.tz;
        if (r_sq > maxRadiusSq) maxRadiusSq = r_sq;
        gridBuilt = false;
    }
    double getMaxRadiusSq() const { return maxRadiusSq; }

    // Построить bucket grid. Вызывать ПОСЛЕ окончания всех addElement.
    // Возвращает суммарное число (elem, bucket) пар (для диагностики).
    int buildSpatialIndex() {
        if (elements.empty()) { gridBuilt = true; return 0; }
        // 1. Определим максимальный worldRadius -- задаём по нему bucket size.
        double maxR = 0.0;
        for (const auto& el : elements) maxR = std::max(maxR, el.worldRadius);
        // bucketSize = 2*maxR гарантирует, что элемент влазит в 2x2 ≤ 4 buckets.
        // Но слишком крупный bucket снижает прескрининг. Берём max(1.0, 2*maxR).
        bucketSize = std::max(1.0, 2.0 * maxR);
        // 2. Глобальные границы (по xx, yy всех bounding sphere).
        double xLo =  1e30, xHi = -1e30;
        double yLo =  1e30, yHi = -1e30;
        for (const auto& el : elements) {
            xLo = std::min(xLo, el.worldCenter.x - el.worldRadius);
            xHi = std::max(xHi, el.worldCenter.x + el.worldRadius);
            yLo = std::min(yLo, el.worldCenter.y - el.worldRadius);
            yHi = std::max(yHi, el.worldCenter.y + el.worldRadius);
        }
        bxMin = floorDiv(xLo, bucketSize);
        byMin = floorDiv(yLo, bucketSize);
        int bxMax = floorDiv(xHi, bucketSize);
        int byMax = floorDiv(yHi, bucketSize);
        bxN = bxMax - bxMin + 1;
        byN = byMax - byMin + 1;
        bucketGrid.assign((size_t)bxN * byN, {});
        // 3. Регистрируем каждый элемент во всех bucket'ах, пересекаемых
        //    его bounding sphere (по проекции на (x,y) plane).
        int totalRegs = 0;
        for (int i = 0; i < (int)elements.size(); ++i) {
            const auto& el = elements[i];
            int bxLo = floorDiv(el.worldCenter.x - el.worldRadius, bucketSize) - bxMin;
            int bxHi = floorDiv(el.worldCenter.x + el.worldRadius, bucketSize) - bxMin;
            int byLo = floorDiv(el.worldCenter.y - el.worldRadius, bucketSize) - byMin;
            int byHi = floorDiv(el.worldCenter.y + el.worldRadius, bucketSize) - byMin;
            for (int bx = bxLo; bx <= bxHi; ++bx)
                for (int by = byLo; by <= byHi; ++by) {
                    bucketGrid[(size_t)bx * byN + by].push_back(i);
                    totalRegs++;
                }
        }
        gridBuilt = true;
        return totalRegs;
    }

    // Геттер для диагностики
    void getGridStats(int& nBucketsX, int& nBucketsY, double& bucketSizeOut,
                      int& nNonEmpty, int& maxBucketLen, double& avgBucketLen) const {
        nBucketsX = bxN; nBucketsY = byN; bucketSizeOut = bucketSize;
        nNonEmpty = 0; maxBucketLen = 0;
        long long sum = 0;
        for (const auto& b : bucketGrid) {
            if (!b.empty()) { nNonEmpty++; sum += b.size();
                              if ((int)b.size() > maxBucketLen) maxBucketLen = b.size(); }
        }
        avgBucketLen = nNonEmpty > 0 ? (double)sum / nNonEmpty : 0.0;
    }

    // Возвращает true, если в окрестности pos НЕТ магнитов (drift section).
    // Можно использовать для перехода на аналитический drift `x += v*dt`
    // вместо полного Boris step. Дёшевая проверка: 2 деления + 1 ветка.
    inline bool isDriftRegion(const FieldVector& pos) const {
        if (!gridBuilt) return false;  // safe fallback
        int bx = floorDiv(pos.x, bucketSize) - bxMin;
        int by = floorDiv(pos.y, bucketSize) - byMin;
        if (bx < 0 || bx >= bxN || by < 0 || by >= byN) return true;
        return bucketGrid[(size_t)bx * byN + by].empty();
    }

    void getTotalFields(const FieldVector& pos, FieldVector& B, FieldVector& E) const {
        B = {0,0,0}; E = {0,0,0};
        if (!gridBuilt) {
            // Fallback: линейный цикл (для backward compat если индекс не построен).
            for (const auto& el : elements) el.addContribution(pos, B, E);
            return;
        }
        // Найти bucket текущей позиции.
        int bx = floorDiv(pos.x, bucketSize) - bxMin;
        int by = floorDiv(pos.y, bucketSize) - byMin;
        if (bx < 0 || bx >= bxN || by < 0 || by >= byN) return;
        const auto& bucket = bucketGrid[(size_t)bx * byN + by];
        for (int idx : bucket) {
            elements[idx].addContribution(pos, B, E);
        }
    }

    // Доступ к элементам для оптимизатора (мутирует currentScale/correctCoeff между итерациями).
    std::vector<MagneticElement>& elementsMut() { return elements; }
    const std::vector<MagneticElement>& elementsRef() const { return elements; }
    size_t size() const { return elements.size(); }
};

#endif
