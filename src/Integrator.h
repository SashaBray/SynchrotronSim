#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include "LatticeManager.h"
#include <cmath>

// Константы СИ
const double C_LIGHT = 299792458.0;
const double C_SQ = C_LIGHT * C_LIGHT;

// Расширим FieldVector базовой математикой для удобства RK4
inline FieldVector operator+(const FieldVector& a, const FieldVector& b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
inline FieldVector operator-(const FieldVector& a, const FieldVector& b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
inline FieldVector operator*(const FieldVector& v, double s) { return {v.x * s, v.y * s, v.z * s}; }
inline FieldVector operator*(double s, const FieldVector& v) { return v * s; }
inline FieldVector operator/(const FieldVector& v, double s) { return {v.x / s, v.y / s, v.z / s}; }

// Векторное произведение [A x B]
inline FieldVector cross(const FieldVector& a, const FieldVector& b) {
    return {
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    };
}

// Текущее состояние частицы: Координаты и Импульс
struct ParticleState {
    FieldVector pos; // (м)
    FieldVector mom; // Импульс p = m0 * gamma * v (кг*м/с)
};

// Производные для RK4
struct Derivatives {
    FieldVector d_pos; // Скорость v
    FieldVector d_mom; // Сила Лоренца F
};

class Integrator {
public:
    /**
     * @brief Вычисление производных в текущей точке
     */
    static Derivatives getDerivatives(const ParticleState& s, double m0, double q, const AcceleratorConfig& lattice) {
        // 1. Расчет Лоренц-фактора через импульс: gamma = sqrt(1 + p^2 / (m0^2 * c^2))
        double p_sq = s.mom.x * s.mom.x + s.mom.y * s.mom.y + s.mom.z * s.mom.z;
        double m0c_sq = (m0 * m0) * C_SQ;
        double gamma = std::sqrt(1.0 + p_sq / m0c_sq);

        // 2. Текущая скорость v = p / (m0 * gamma)
        FieldVector v = s.mom / (m0 * gamma);

        // 3. Получение суммарных полей
        FieldVector B, E;
        lattice.getTotalFields(s.pos, B, E);

        // 4. Уравнения движения: dx/dt = v, dp/dt = q(E + [v x B])
        Derivatives d;
        d.d_pos = v;
        d.d_mom = q * (E + cross(v, B));

        return d;
    }

    /**
     * @brief Один шаг методом Рунге-Кутты 4-го порядка
     */
    static void rk4Step(ParticleState& s, double m0, double q, const AcceleratorConfig& lattice, double dt) {
        Derivatives k1 = getDerivatives(s, m0, q, lattice);

        ParticleState s2 = { s.pos + k1.d_pos * (dt * 0.5), s.mom + k1.d_mom * (dt * 0.5) };
        Derivatives k2 = getDerivatives(s2, m0, q, lattice);

        ParticleState s3 = { s.pos + k2.d_pos * (dt * 0.5), s.mom + k2.d_mom * (dt * 0.5) };
        Derivatives k3 = getDerivatives(s3, m0, q, lattice);

        ParticleState s4 = { s.pos + k3.d_pos * dt, s.mom + k3.d_mom * dt };
        Derivatives k4 = getDerivatives(s4, m0, q, lattice);

        // Итоговое обновление состояния
        s.pos = s.pos + (dt / 6.0) * (k1.d_pos + 2.0 * k2.d_pos + 2.0 * k3.d_pos + k4.d_pos);
        s.mom = s.mom + (dt / 6.0) * (k1.d_mom + 2.0 * k2.d_mom + 2.0 * k3.d_mom + k4.d_mom);
    }

    /**
     * @brief Вспомогательная функция для пересчета начальной скорости в импульс
     */
    static FieldVector velocityToMomentum(FieldVector v, double m0) {
        double v_sq = v.x * v.x + v.y * v.y + v.z * v.z;
        double gamma = 1.0 / std::sqrt(1.0 - v_sq / C_SQ);
        return v * (m0 * gamma);
    }
};

#endif // INTEGRATOR_H
