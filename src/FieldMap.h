#ifndef FIELD_MAP_H
#define FIELD_MAP_H

#include "AnsysConverter.h" // Используем GridHeader и FieldVector
#include <vector>
#include <string>
#include <fstream>
#include <memory>
#include <cmath>

class FieldMap {
private:
    GridHeader header;
    std::vector<FieldVector> data;

public:
    FieldMap() = default;

    /**
     * @brief Загрузка бинарной карты в память
     */
    bool loadFromBinary(const std::string& path) {
        std::ifstream in(path, std::ios::binary);
        if (!in) return false;

        // Читаем заголовок
        in.read(reinterpret_cast<char*>(&header), sizeof(GridHeader));

        // Выделяем память и читаем весь массив данных одним блоком
        size_t totalNodes = static_cast<size_t>(header.nx) * header.ny * header.nz;
        data.resize(totalNodes);
        in.read(reinterpret_cast<char*>(data.data()), totalNodes * sizeof(FieldVector));

        return in.good();
    }

    /**
     * @brief Трилинейная интерполяция поля в локальных координатах
     * @param p Точка в локальной системе координат (метры)
     * @param useSymmetry Флаг использования симметрии YZ (относительно X=0)
     */
    FieldVector getInterpolatedField(FieldVector p, bool useSymmetry,
                                     double rotX = 0.0, double rotY = 0.0, double rotZ = 0.0) const {
        bool mirrored = false;
        const bool hasRot = (rotX != 0.0 || rotY != 0.0 || rotZ != 0.0);

        // Поворот 3-вектора вокруг локальной оси (ax: 0=x,1=y,2=z) на угол ang (рад).
        auto rotAxis = [](FieldVector v, int ax, double ang) -> FieldVector {
            double c = std::cos(ang), s = std::sin(ang);
            FieldVector o = v;
            if (ax == 0)      { o.y = v.y * c - v.z * s; o.z = v.y * s + v.z * c; }  // Rx
            else if (ax == 1) { o.z = v.z * c - v.x * s; o.x = v.z * s + v.x * c; }  // Ry
            else              { o.x = v.x * c - v.y * s; o.y = v.x * s + v.y * c; }  // Rz
            return o;
        };

        // 0. Поворот карты (переориентация установленного магнита). R = Rz·Ry·Rx (углы rotX/Y/Z, рад).
        //    Поле ПОВЁРНУТОГО магнита в точке p: опрашиваем ИСХОДНОЕ поле в Rᵀp, затем поворачиваем
        //    вектор поля B(p)=R·B_исх(Rᵀp). Поворот — собственное преобразование (det=+1), псевдовектор
        //    B при нём преобразуется как обычный вектор (в отличие от отражения). Rᵀ = Rxᵀ·Ryᵀ·Rzᵀ.
        if (hasRot) {
            p = rotAxis(p, 2, -rotZ);
            p = rotAxis(p, 1, -rotY);
            p = rotAxis(p, 0, -rotX);
        }

        // 1. Учет симметрии по плоскости YZ (нормаль - ось X)
        if (useSymmetry && p.x < 0) {
            p.x = -p.x; // Отражаем координату для поиска в положительной сетке
            mirrored = true;
        }

        // 2. Проверка выхода за границы сетки (Bounding Box)
        if (p.x < header.minX || p.x > header.maxX ||
            p.y < header.minY || p.y > header.maxY ||
            p.z < header.minZ || p.z > header.maxZ) {
            return {0.0, 0.0, 0.0};
        }

        // 3. Расчет дробных индексов
        double fx = (p.x - header.minX) / header.stepX;
        double fy = (p.y - header.minY) / header.stepY;
        double fz = (p.z - header.minZ) / header.stepZ;

        int i = static_cast<int>(fx);
        int j = static_cast<int>(fy);
        int k = static_cast<int>(fz);

        // Веса для интерполяции (0.0 ... 1.0)
        double tx = fx - i;
        double ty = fy - j;
        double tz = fz - k;

        // Защита от выхода за пределы массива при p.x == maxX
        if (i >= header.nx - 1) { i = header.nx - 2; tx = 1.0; }
        if (j >= header.ny - 1) { j = header.ny - 2; ty = 1.0; }
        if (k >= header.nz - 1) { k = header.nz - 2; tz = 1.0; }

        // 4. Извлечение 8 соседних узлов
        auto getV = [&](int _i, int _j, int _k) -> const FieldVector& {
            return data[static_cast<size_t>(_i) * (header.ny * header.nz) +
                        static_cast<size_t>(_j) * header.nz + _k];
        };

        const auto& v000 = getV(i,   j,   k);
        const auto& v100 = getV(i+1, j,   k);
        const auto& v010 = getV(i,   j+1, k);
        const auto& v110 = getV(i+1, j+1, k);
        const auto& v001 = getV(i,   j,   k+1);
        const auto& v101 = getV(i+1, j,   k+1);
        const auto& v011 = getV(i,   j+1, k+1);
        const auto& v111 = getV(i+1, j+1, k+1);

        // 5. Линейная интерполяция по трем осям
        FieldVector res;
        auto lerp = [](double a, double b, double t) { return a * (1.0 - t) + b * t; };

        res.x = lerp(lerp(lerp(v000.x, v100.x, tx), lerp(v010.x, v110.x, tx), ty),
                     lerp(lerp(v001.x, v101.x, tx), lerp(v011.x, v111.x, tx), ty), tz);

        res.y = lerp(lerp(lerp(v000.y, v100.y, tx), lerp(v010.y, v110.y, tx), ty),
                     lerp(lerp(v001.y, v101.y, tx), lerp(v011.y, v111.y, tx), ty), tz);

        res.z = lerp(lerp(lerp(v000.z, v100.z, tx), lerp(v010.z, v110.z, tx), ty),
                     lerp(lerp(v001.z, v101.z, tx), lerp(v011.z, v111.z, tx), ty), tz);

        // 6. Применяем правило симметрии: инвертируем Bx, если мы в отзеркаленной зоне
        if (mirrored) {
            res.x = -res.x;
        }

        // 7. Поворот вектора поля R·B_исх: R = Rz·Ry·Rx (сначала Rx, затем Ry, затем Rz).
        if (hasRot) {
            res = rotAxis(res, 0, rotX);
            res = rotAxis(res, 1, rotY);
            res = rotAxis(res, 2, rotZ);
        }

        return res;
    }

    double getNominalArg() const { return header.nominalArg; }

    // Доступ к границам сетки -- нужен LatticeManager'у для расчёта
    // bounding sphere магнита (быстрый прескрининг до полного transform).
    const GridHeader& headerRef() const { return header; }
};

#endif // FIELD_MAP_H
