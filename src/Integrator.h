#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include "LatticeManager.h"
#include <cmath>

// --- Физические константы СИ ---
const double C_LIGHT = 299792458.0;
const double C_SQ    = C_LIGHT * C_LIGHT;

// --- Векторная арифметика над FieldVector ---
inline FieldVector operator+(const FieldVector& a, const FieldVector& b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
inline FieldVector operator-(const FieldVector& a, const FieldVector& b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
inline FieldVector operator*(const FieldVector& v, double s)             { return {v.x * s, v.y * s, v.z * s}; }
inline FieldVector operator*(double s, const FieldVector& v)             { return v * s; }
inline FieldVector operator/(const FieldVector& v, double s)             { return {v.x / s, v.y / s, v.z / s}; }

// Векторное произведение [A x B]
inline FieldVector cross(const FieldVector& a, const FieldVector& b) {
    return {
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    };
}

// Текущее состояние частицы: координаты и импульс
struct ParticleState {
    FieldVector pos;  // (м)
    FieldVector mom;  // импульс p = m0 * gamma * v (кг*м/с)
};

class Integrator {
public:
    /**
     * @brief Один шаг релятивистского Boris pusher (drift-kick-drift).
     *
     * Свойства:
     *   - Симплектичен (сохраняет фазовый объём).
     *   - Точно сохраняет |p|² для чисто магнитного поля (E = 0):
     *     поворот в шаге 5 -- ортогональная операция.
     *   - Один вызов getTotalFields() на шаг (RK4 потребовал бы 4).
     *   - Симметричный во времени, формальный порядок O(dt²).
     *
     * Литература: Boris, Proc. 4th Conf. on Numerical Simulation of
     *             Plasmas, 1970; Birdsall & Langdon, "Plasma Physics via
     *             Computer Simulation", 1985, §4-3.
     *
     * Замечание: в текущей конфигурации все элементы магнитные, поэтому
     * E ≡ 0 и шаги 3, 6 фактически no-op. Код оставлен в общем виде на
     * случай, если в будущем появятся электрические элементы.
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

        // 5. Boris magnetic rotation p^- -> p^+ (норма сохраняется точно)
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
     * @brief Пересчёт начальной скорости в релятивистский импульс.
     */
    static FieldVector velocityToMomentum(FieldVector v, double m0) {
        double v_sq  = v.x * v.x + v.y * v.y + v.z * v.z;
        double gamma = 1.0 / std::sqrt(1.0 - v_sq / C_SQ);
        return v * (m0 * gamma);
    }
};

#endif // INTEGRATOR_H
