#ifndef LATTICE_MANAGER_H
#define LATTICE_MANAGER_H

#include "FieldMap.h"
#include <vector>
#include <memory>
#include <cmath>

// Типы полей для идентификации в расчетах
enum class FieldType { MAGNETIC, ELECTRIC };

// Легкая структура для матричных операций (3x3 и смещение)
struct TransformMatrix {
    double m[3][3]; // Матрица вращения
    double tx, ty, tz; // Вектор смещения (центр устройства)

    // Применение трансформации к точке: P_loc = R^T * (P_glob - T)
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

    // Применение вращения к вектору: V_glob = R * V_loc
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
    std::shared_ptr<FieldMap> map; // Ссылка на общую карту поля
    FieldType type;
    TransformMatrix transform;    // Геометрия в пространстве
    double currentScale;          // Коэффициент масштабирования (I_real / I_map)
    bool useSymmetry;

    /**
     * @brief Инициализация матриц на основе углов Эйлера (Yaw-Pitch-Roll)
     */
    void setOrientation(double x, double y, double z, double yaw, double pitch, double roll) {
        transform.tx = x; transform.ty = y; transform.tz = z;

        // Матрицы поворота (X-вперед, Y-влево, Z-вверх)
        // Yaw (psi) вокруг Z, Pitch (theta) вокруг Y, Roll (phi) вокруг X
        double sy = sin(yaw),   cy = cos(yaw);
        double sp = sin(pitch), cp = cos(pitch);
        double sr = sin(roll),  cr = cos(roll);

        // Матрица R = Rz * Ry * Rx
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

    /**
     * @brief Получение вклада данного элемента в глобальной точке
     */
    void addContribution(const FieldVector& globalPos, FieldVector& totalB, FieldVector& totalE) const {
        // 1. Перевод в локальную систему координат
        FieldVector localPos = transform.transformPointInv(globalPos);

        // 2. Запрос интерполированного поля из карты
        FieldVector localF = map->getInterpolatedField(localPos, useSymmetry);

        // Если поле нулевое (вне границ), выходим
        if (localF.x == 0 && localF.y == 0 && localF.z == 0) return;

        // 3. Масштабирование и поворот вектора обратно в глобальную систему
        FieldVector globalF = transform.rotateVector(localF);
        globalF.x *= currentScale;
        globalF.y *= currentScale;
        globalF.z *= currentScale;

        // 4. Суммирование в нужный канал
        if (type == FieldType::MAGNETIC) {
            totalB.x += globalF.x; totalB.y += globalF.y; totalB.z += globalF.z;
        } else {
            totalE.x += globalF.x; totalE.y += globalF.y; totalE.z += globalF.z;
        }
    }
};

class AcceleratorConfig {
private:
    std::vector<MagneticElement> elements;
    double maxRadiusSq = 0;

public:
    void addElement(const MagneticElement& el) {
        elements.push_back(el);
        // Обновляем максимальный радиус для детектора вылета
        double r_sq = el.transform.tx * el.transform.tx + 
                      el.transform.ty * el.transform.ty + 
                      el.transform.tz * el.transform.tz;
        if (r_sq > maxRadiusSq) maxRadiusSq = r_sq;
    }

    double getMaxRadiusSq() const { return maxRadiusSq; }

    /**
     * @brief Расчет суммарных полей во всей установке
     */
    void getTotalFields(const FieldVector& pos, FieldVector& B, FieldVector& E) const {
        B = {0,0,0}; E = {0,0,0};
        for (const auto& el : elements) {
            el.addContribution(pos, B, E);
        }
    }
};

#endif // LATTICE_MANAGER_H
