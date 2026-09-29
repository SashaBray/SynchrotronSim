// matrix_core.h — ФИЗИЧЕСКИЙ ДВИЖОК матричного ядра SynchrotronSim на C++.
//
// Дословный порт локальной физики Python-ядра matrix_tracker.py (симплектический секторный интегратор —
// порт elegant CSBEND integrate_csbend_ordn — + мультипольный кик KQUAD/KSEXT/KOCT + распределённое
// синхротронное излучение). Канонические координаты (x, px, y, py, t, delta): px,py нормированы на p0,
// delta=Δp/p0. Порядок операций сохранён как в Python, чтобы результат совпадал с эталоном
// (matrix_tracker.py) — см. регрессию test_matrix_core_cpp.py. ИСКЛЮЧЕНИЕ (оптимизация скорости):
// секторный дрейф использует ТОЖДЕСТВО сложения углов (asin/sin/cos → sqrt + прекэш sinA/cosA) — физически
// эквивалентно, но не побитово (расхождение ~1e-14/элемент, ~3e-12/оборот кольца; тесты проходят).
//
// Python-ядро остаётся эталоном для разработки/опытов; нововведения переносятся сюда и валидируются.
#pragma once
#include <cmath>
#include <string>
#include <vector>
#include <random>
#include <stdexcept>

namespace mc {

// ---- физические константы ПУЧКА (сурогат) — дословно как bip_matrix.py -------------------------
constexpr double C_LIGHT = 299792458.0;
constexpr double ME      = 9.1093837e-31;     // масса электрона, кг
constexpr double E_CH    = 1.602176634e-19;   // элементарный заряд, Кл
constexpr double Q_E     = -E_CH;             // заряд электрона
constexpr double V_BEAM  = 2.967e8;           // скорость сурогатного пучка, м/с (particle_groups.csv)

// q/p0 (опорная жёсткость). Считается ровно как bip_matrix.QP — beta,gamma,p_mom,QP той же формулой.
inline double qp() {
    double beta  = V_BEAM / C_LIGHT;
    double gamma = 1.0 / std::sqrt(1.0 - beta * beta);
    double p_mom = gamma * ME * V_BEAM;
    return Q_E / p_mom;
}

// ---- константы РАДИАЦИИ — как matrix_tracker.py ------------------------------------------------
constexpr double MC2_MEV       = 0.51099895069;              // энергия покоя электрона, МэВ
constexpr double R_E           = 2.8179403262e-15;           // классический радиус электрона, м
constexpr double RAD_K         = 2.0 / 3.0 * R_E;            // rad_coef = (2/3)·r_e·γ³/ρ²
// квантовое возбуждение: isr_coef = ISR_K·γ^2.5/|ρ|^1.5 — дословно matrix_tracker._ISR_K
const double ISR_K = R_E * std::sqrt(55.0 / (24.0 * std::sqrt(3.0)) * 137.0359895);
constexpr double YOSHIDA_BETA  = 1.25992104989487316477;     // 2^(1/3) — расщепление Йошиды 4-го порядка
constexpr int    SECTOR_NSLICES = 5;                        // срезов сектор/кик-тела (_nsl_sector, как elegant N_SLICES=5)

struct ParticleLost : std::runtime_error {
    explicit ParticleLost(const std::string& m) : std::runtime_error(m) {}
};

struct State { double x, px, y, py, t, delta; };

// Коэффициенты (дрейф, кик) симплектического расщепления Йошиды порядка 2/4/6 — дословно
// matrix_tracker._yoshida (а тот — из elegant integrate_csbend_ordn). Последний кик шага = 0.
struct Yoshida {
    std::vector<double> d, k;
    explicit Yoshida(int order = 4) {
        if (order == 2) { d = {0.5, 0.5}; k = {1.0, 0.0}; return; }
        if (order == 4) {
            double b = YOSHIDA_BETA;
            d = {0.5 / (2 - b), (1 - b) / (2 - b) / 2, (1 - b) / (2 - b) / 2, 0.5 / (2 - b)};
            k = {1.0 / (2 - b), -b / (2 - b), 1.0 / (2 - b), 0.0};
            return;
        }
        if (order == 6) {
            d = {0.39225680523878, 0.5100434119184585, -0.47105338540975655, 0.0687531682525181,
                 0.0687531682525181, -0.47105338540975655, 0.5100434119184585, 0.39225680523878};
            k = {0.784513610477560, 0.235573213359357, -1.17767998417887, 1.3151863206839063,
                 -1.17767998417887, 0.235573213359357, 0.784513610477560, 0.0};
            return;
        }
        throw std::runtime_error("порядок интегрирования — допустимы 2/4/6");
    }
};
using Yoshida4 = Yoshida;   // прежнее имя (порядок 4 по умолчанию)

// Генератор квантового возбуждения: у каждой частицы свой поток (детерминизм при OpenMP).
struct IsrRng {
    std::mt19937_64 gen;
    std::normal_distribution<double> norm{0.0, 1.0};
    explicit IsrRng(uint64_t seed) : gen(seed) {}
    double operator()() { return norm(gen); }
};

// Точный дрейф на длину ell (δ-корректный, канонический) — matrix_tracker._drift.
inline State drift(const State& st, double ell) {
    double arg = (1.0 + st.delta) * (1.0 + st.delta) - st.px * st.px - st.py * st.py;
    if (arg <= 0.0) throw ParticleLost("p_perp>p (drift)");
    double ps = std::sqrt(arg);
    return { st.x + ell * st.px / ps, st.px, st.y + ell * st.py / ps, st.py,
             st.t + ell * (1.0 + st.delta) / ps, st.delta };
}

// Нормированное поле (Fx,Fy) комбинированной функции до 4-го порядка — matrix_tracker._csbend_fields.
inline void csbend_fields(double x, double y, double b1, double b2, double b3, double b4, double h,
                          double& Fx, double& Fy) {
    double h2 = h * h, h3 = h2 * h;
    double x2 = x * x, y2 = y * y;
    Fy = (1.0
          + b1 * x + 0.5 * b2 * x2 + (b3 / 6.0) * x2 * x + (b4 / 24.0) * x2 * x2
          + (-(h * b1) / 2.0 - b2 / 2.0) * y2
          + ((h2 * b1) / 2.0 - (h * b2) / 2.0 - b3 / 2.0) * x * y2
          + (-(h3 * b1) / 2.0 + (h2 * b2) / 2.0 - (h * b3) / 4.0 - b4 / 4.0) * x2 * y2
          + ((h3 * b1) / 24.0 - (h2 * b2) / 24.0 + (h * b3) / 12.0 + b4 / 24.0) * y2 * y2);
    Fx = (b1 * y + b2 * x * y + (b3 / 2.0) * x2 * y + (b4 / 6.0) * x2 * x * y
          + ((h2 * b1) / 6.0 - (h * b2) / 6.0 - b3 / 6.0) * y2 * y
          + (-(h3 * b1) / 3.0 + (h2 * b2) / 3.0 - (h * b3) / 6.0 - b4 / 6.0) * x * y2 * y);
}

struct Element {
    std::string kind;
    double Leff = 0.0;
    int    nsl  = 5;
    double theta0 = 0, k1 = 0, k1v = 0, k2 = 0, k3 = 0;   // силы по опорной жёсткости
    double kEdge1 = 0, kEdge2 = 0;                        // краевые линзы E1/E2
    bool   sector = false;                                // изгиб (dipole/dq) с theta0≠0 → секторный интегратор
    bool   kick   = false;                                // прямой мультиполь (quad/sext/oct)
    double rad_coef = 0;                                  // распределённое СИ (классическое)
    double isr_coef = 0;                                  // квантовое возбуждение (ISR), нужен IsrRng
    double hkick = 0;                                     // однородный ГОР. кик корректора, рад (возмущение/коррекция)
    double vkick = 0;                                     // однородный ВЕРТ. кик корректора, рад
    int    nsl_box = 40;                                  // срезов «ящика» из БД (n_slices) — для кика корректора
    double rho0 = 0, h = 0, b1 = 0, b2 = 0, b3 = 0;       // ПРЕКЭШ сектора: ρ, 1/ρ, bₙ=kₙ·ρ (не пересчитывать в цикле)
    double s_slice = 0;                                   // длина среза Leff/nsl (сектор/кик)
    std::vector<double> sinA, cosA;                       // ПРЕКЭШ sin/cos угла стадии dsh/ρ (замена asin/sin/cos в дрейфе)
    Yoshida ys;

    // Сборка из физкоэф. mat_* (B,G,Gv,S,O — суммы обмоток) + Leff + проектной энергии — как Magnet.__init__.
    // order — порядок схемы Йошиды (2/4/6), isr — квантовое возбуждение (только вместе с sr).
    static Element make(const std::string& kind, double Leff, double B, double G, double Gv,
                        double S, double O, double edge1, double edge2, int nsl,
                        double design_mev, bool sr, int order = 4, bool isr = false, int nsl_box = 40) {
        Element e; e.kind = kind; e.Leff = Leff; e.nsl = nsl; e.ys = Yoshida(order); e.nsl_box = nsl_box;
        double b = B, g = G, gv = Gv, s = S, o = O;
        if (kind == "dipole") { g *= b; gv *= b; s *= b; o *= b; }   // градиенты диполя заданы относит. поля
        double QP = qp();
        e.theta0 = QP * b * Leff;
        e.k1 = QP * g; e.k1v = QP * gv; e.k2 = QP * s; e.k3 = QP * o;
        double rinv = (Leff > 0) ? e.theta0 / Leff : 0.0;
        e.kEdge1 = std::tan(edge1) * rinv;
        e.kEdge2 = std::tan(edge2) * rinv;
        bool bend = (kind == "dipole" || kind == "dq") && std::fabs(e.theta0) > 0.0;
        e.sector = bend;
        bool nonzero = e.theta0 || e.k1 || e.k1v || e.k2 || e.k3 || e.kEdge1 || e.kEdge2;
        e.kick = (!bend) && nonzero && (kind == "quad" || kind == "sext" || kind == "oct");
        if (e.sector || e.kick) e.s_slice = Leff / e.nsl;           // ПРЕКЭШ длины среза
        e.precache();
        if (e.sector && sr) {                                        // распределённое СИ (классическое)
            double gamma = design_mev / MC2_MEV;
            e.rad_coef = RAD_K * gamma * gamma * gamma / (e.rho0 * e.rho0);
            if (isr) e.isr_coef = ISR_K * std::pow(gamma, 2.5) / std::pow(std::fabs(e.rho0), 1.5);
        }
        return e;
    }

    // ПРЕКЭШ сектора (ρ, 1/ρ, bₙ, sin/cos стадий) — пересчитывается и у срезов изгиба (_bend_slices).
    void precache() {
        if (sector || kick) s_slice = Leff / nsl;
        if (sector) {
            rho0 = -Leff / theta0; h = 1.0 / rho0;                   // ρ>0 (конвенция elegant)
            b1 = k1 * rho0; b2 = k2 * rho0; b3 = k3 * rho0;
            sinA.assign(ys.d.size(), 0.0); cosA.assign(ys.d.size(), 0.0);
            for (size_t j = 0; j < ys.d.size(); ++j) {              // sin/cos угла стадии a=dsh/ρ=s·d[j]/ρ (постоянно)
                double a = s_slice * ys.d[j] / rho0;
                sinA[j] = std::sin(a); cosA[j] = std::cos(a);
            }
        }
    }

    // Диспетчер тела: изгиб → секторный интегратор; прямой мультиполь → кик; корректор с киком → «ящик»
    // (дрейф–кик–дрейф, n_slices из БД — как box-тело matrix_tracker для θ0 у неизгибающего); иначе дрейф.
    State body(const State& st, IsrRng* rng = nullptr) const {
        if (sector) return body_sector(st, rng);
        if ((hkick != 0.0 || vkick != 0.0) && !kick) return body_kicker(st);
        if (kick)   return body_kick(st);
        return drift(st, Leff);
    }

private:
    // Секторный интегратор изгиба (+ классическое СИ) — matrix_tracker._body_sector_symplectic.
    State body_sector(const State& st, IsrRng* rng) const {
        double x = st.x, px = st.px, y = st.y, py = st.py, t = st.t, d = st.delta;
        const double b4 = 0.0;                     // rho0/h/b1..b3/s_slice — ПРЕКЭШ (Element::make), не пересчёт
        int N = nsl; double s = s_slice;
        px += kEdge1 * x; py -= kEdge1 * y;        // вход: краевой фокус E1
        double X = x, QX = px, Y = y, QY = py, DP = d;
        double dist = 0.0;
        const int J = (int)ys.d.size();
        const double ds_isr = Leff / (N * std::max(1, J - 1));       // длина под-шага ISR (nSubsteps−1, как elegant)
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < J; ++j) {
                // -- точный секторный ДРЕЙФ (по кругу); asin/sin/cos → sqrt + сложение углов (прекэш sinA/cosA[j])
                double f = (1.0 + DP) * (1.0 + DP) - QY * QY;
                if (f <= 0.0) throw ParticleLost("sector: p_perp>p (drift)");
                f = std::sqrt(f);
                double sin_phi = QX / f;                              // φ=asin(QX/f)∈[−π/2,π/2]
                if (std::fabs(sin_phi) > 1.0) throw ParticleLost("sector: |QX/f|>1 (drift)");
                double cos_phi = std::sqrt(1.0 - sin_phi * sin_phi);  // =cos(asin(QX/f))
                double sine = sinA[j] * cos_phi + cosA[j] * sin_phi;  // sin(dsh/ρ+φ) — тождество сложения
                double cosi = cosA[j] * cos_phi - sinA[j] * sin_phi;  // cos(dsh/ρ+φ)
                if (cosi == 0.0) throw ParticleLost("sector: cos=0 (drift)");
                QX = f * sine;
                double factor = (rho0 + X) * cos_phi / f * (sine / cosi - sin_phi / cos_phi);
                Y += QY * factor;
                dist += factor * (1.0 + DP);
                double fr = cos_phi / cosi;
                X = rho0 * (fr - 1.0) + fr * X;
                if (ys.k[j] == 0.0) break;
                double ds = s * ys.k[j];           // -- КИК поля с фактором кривизны (1+X/ρ) --
                double Fx, Fy; csbend_fields(X, Y, b1, b2, b3, b4, h, Fx, Fy);
                QX += -ds * (1.0 + X / rho0) * Fy / rho0;
                QY +=  ds * (1.0 + X / rho0) * Fx / rho0;
                if (rad_coef != 0.0 || isr_coef != 0.0) {   // -- РАДИАЦИЯ: классическая + квантовая --
                    double denom = (1.0 + DP) * (1.0 + DP) - QX * QX - QY * QY;
                    if (denom <= 0.0) throw ParticleLost("sector: p_perp>p (radiation)");
                    double ff = (1.0 + X / rho0) / std::sqrt(denom);
                    double xp = QX * ff, yp = QY * ff;
                    double ds_factor = std::sqrt((1.0 + X / rho0) * (1.0 + X / rho0) + xp * xp + yp * yp);
                    double F2 = Fx * Fx + Fy * Fy;
                    double delta_factor = (1.0 + DP) * (1.0 + DP);
                    QX /= (1.0 + DP); QY /= (1.0 + DP);           // канонич. → наклон
                    if (rad_coef != 0.0) DP -= rad_coef * delta_factor * F2 * ds * ds_factor;
                    if (isr_coef != 0.0 && rng)                   // квантовый скачок энергии (гауссов шум)
                        DP -= isr_coef * delta_factor * std::pow(F2, 0.75)
                              * std::sqrt(ds_isr * ds_factor) * (*rng)();
                    QX *= (1.0 + DP); QY *= (1.0 + DP);           // наклон → канонич. при новой энергии
                }
            }
        }
        x = X; px = QX; y = Y; py = QY; d = DP;
        px += kEdge2 * x; py -= kEdge2 * y;        // выход: краевой фокус E2
        return { x, px, y, py, t + dist, d };
    }

    // Корректор (кик θ у неизгибающего элемента) — box-тело matrix_tracker при θ0≠0 и нулевых мультиполях:
    // n_slices раз дрейф(ds/2) → кик θ/N → дрейф(ds/2).
    State body_kicker(const State& st) const {
        State p = st;
        int N = nsl_box; double ds = Leff / N;
        for (int i = 0; i < N; ++i) {
            p = drift(p, ds / 2.0);
            p.px += hkick / N; p.py += vkick / N;
            p = drift(p, ds / 2.0);
        }
        return p;
    }

    // Прямой мультиполь (KQUAD/KSEXT/KOCT) — matrix_tracker._body_kick_symplectic.
    State body_kick(const State& st) const {
        State p = st;
        int N = nsl; double ds = s_slice;          // Leff/nsl — ПРЕКЭШ
        const int J = (int)ys.d.size();
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < J; ++j) {
                p = drift(p, ds * ys.d[j]);        // дрейф по слопу
                if (ys.k[j] == 0.0) break;
                double kf = ds * ys.k[j];
                double x = p.x, y = p.y;
                p.px += (-k1 * x - 0.5 * k2 * (x * x - y * y)
                         - (k3 / 6.0) * (x * x * x - 3.0 * x * y * y)) * kf;
                p.py += (k1v * y + k2 * x * y
                         + (k3 / 6.0) * (3.0 * x * x * y - y * y * y)) * kf;
            }
        }
        return p;
    }
};

} // namespace mc
