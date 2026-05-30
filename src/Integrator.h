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
     * @brief Один шаг релятивистского Boris pusher (drift-kick-drift splitting).
     *
     * Преимущества по сравнению с RK4:
     *   1) ТОЧНО сохраняет |p|² для чистого магнитного поля (E=0).
     *      Это математическое тождество: rotation в шаге 5 ортогональна,
     *      она не может изменить норму вектора.  RK4 же дрейфует |p| как
     *      O(dt^5) на шаг, что за 10^6 шагов даёт систематический уход γ
     *      и радиуса Larmor -- т.е. орбита "уплывает".
     *   2) 1 вызов getTotalFields() вместо 4 в RK4 -> ~4× ускорение.
     *   3) Симплектический: сохраняет фазовый объём, орбиты замкнуты.
     *   4) Симметричный по времени (DKD), формальный порядок 2 по dt.
     *
     * Алгоритм (см. Birdsall & Langdon, "Plasma Physics via Computer
     *           Simulation", 1985, §4-3; Boris 1970):
     *
     *   1. drift на dt/2 ->  x_{n+1/2}
     *   2. поле в x_{n+1/2} -> B, E
     *   3. E half-kick:        p^- = p_n + qE·dt/2
     *   4. γ^- = sqrt(1 + |p^-|² / (m₀²c²))
     *   5. magnetic rotation:  t = qB/(γ^- m₀)·dt/2;   s = 2t/(1+|t|²)
     *                          p' = p^- + p^- × t
     *                          p^+ = p^- + p' × s     (|p^+| = |p^-|)
     *   6. E half-kick:        p_{n+1} = p^+ + qE·dt/2
     *   7. drift на dt/2 ->  x_{n+1}
     */
    static void borisStep(ParticleState& s, double m0, double q,
                          const AcceleratorConfig& lattice, double dt) {
        const double dt_half  = 0.5 * dt;
        const double inv_m0   = 1.0 / m0;
        const double inv_m0c2 = 1.0 / ((m0 * m0) * C_SQ);

        // 1. half-drift: x_n -> x_{n+1/2}
        double p2     = s.mom.x * s.mom.x + s.mom.y * s.mom.y + s.mom.z * s.mom.z;
        double gamma  = std::sqrt(1.0 + p2 * inv_m0c2);
        double inv_mg = inv_m0 / gamma;
        s.pos = s.pos + s.mom * (inv_mg * dt_half);

        // 2. fields at half-step position
        FieldVector B, E;
        lattice.getTotalFields(s.pos, B, E);

        // 3. half E-kick
        FieldVector p_minus = s.mom + E * (q * dt_half);

        // 4. γ at intermediate momentum
        p2 = p_minus.x * p_minus.x + p_minus.y * p_minus.y + p_minus.z * p_minus.z;
        double gamma_minus = std::sqrt(1.0 + p2 * inv_m0c2);

        // 5. Boris magnetic rotation  p^- -> p^+   (норма сохраняется точно)
        double t_coeff = q * dt_half / (m0 * gamma_minus);
        FieldVector t_vec = B * t_coeff;
        double t2 = t_vec.x * t_vec.x + t_vec.y * t_vec.y + t_vec.z * t_vec.z;
        FieldVector s_vec = t_vec * (2.0 / (1.0 + t2));
        FieldVector p_prime = p_minus + cross(p_minus, t_vec);
        FieldVector p_plus  = p_minus + cross(p_prime, s_vec);

        // 6. second half E-kick
        s.mom = p_plus + E * (q * dt_half);

        // 7. half-drift: x_{n+1/2} -> x_{n+1}
        p2     = s.mom.x * s.mom.x + s.mom.y * s.mom.y + s.mom.z * s.mom.z;
        gamma  = std::sqrt(1.0 + p2 * inv_m0c2);
        inv_mg = inv_m0 / gamma;
        s.pos = s.pos + s.mom * (inv_mg * dt_half);
    }

    /**
     * @brief Аналитический drift: x += v * dt, mom не меняется.
     *
     * Используется когда мы знаем, что в окрестности частицы НЕТ магнитов
     * (field = 0). Эквивалентно Boris-step с B=E=0, но в 3-5× дешевле
     * (1 sqrt вместо 3, в 5× меньше арифметики). Сохраняет |v| строго.
     */
    static void driftStep(ParticleState& s, double m0, double dt) {
        double p2     = s.mom.x * s.mom.x + s.mom.y * s.mom.y + s.mom.z * s.mom.z;
        double gamma  = std::sqrt(1.0 + p2 / ((m0 * m0) * C_SQ));
        double k      = dt / (m0 * gamma);
        s.pos.x += s.mom.x * k;
        s.pos.y += s.mom.y * k;
        s.pos.z += s.mom.z * k;
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
