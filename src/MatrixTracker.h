#ifndef MATRIX_TRACKER_H
#define MATRIX_TRACKER_H
// ==========================================================================
// MatrixTracker.h — C++ порт matrix_tracker.py (альтернативный МАТРИЧНЫЙ движок
// SynchrotronSim). Каноническое фазовое пространство (x, pₓ, y, p_y, t, δ):
// pₓ,p_y — нормированные на p₀ поперечные импульсы, δ=Δp/p₀, t — путь (задел под
// ВЧ), s — независимая переменная. p_полный=p₀(1+δ) сохраняется в магнитах.
//
// Модель элемента: прямоугольный ящик, центр=(x,y,z) из lattice_configs, ось по
// yaw, торцы = центр ± L_phys/2·ось, плоскости входа/выхода ⊥ оси. Со-локованные
// обмотки → один магнит (сумма физкоэф.). Кики по опорной жёсткости q/p₀:
// θ₀=(q/p₀)·B·L_eff, k1=(q/p₀)G, k2=(q/p₀)S, k3=(q/p₀)O. Тело — N срезов
// drift-kick-drift: ДРЕЙФ по геометрии (L_phys), КИКИ по полю (L_eff). Модифиц.
// пустой промежуток — геом. пролёт от плоскости выхода пред. магнита до плоскости
// входа след. (произв. ориентации) + проверка направления импульса.
//
// Только стандартная библиотека C++17. Читает те же CSV, что и Python-версия.
// ==========================================================================
#include <vector>
#include <string>
#include <map>
#include <cmath>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <stdexcept>
#include <array>

namespace mtx {

// ---- Опорная жёсткость q/p₀ (константы bip_matrix; совпадает с bm.QP) --------
constexpr double C_LIGHT = 299792458.0;
constexpr double V_BEAM  = 2.967e8;            // м/с — наш пучок (particle_groups)
constexpr double ME      = 9.1093837e-31;      // масса электрона, кг
constexpr double Q_CHG   = -1.602176634e-19;   // заряд электрона −e (CODATA), как в bip_matrix
inline double computeQP() {
    double beta = V_BEAM / C_LIGHT;
    double gamma = 1.0 / std::sqrt(1.0 - beta * beta);
    return Q_CHG / (gamma * ME * V_BEAM);       // = -84.92538656332358 (== bm.QP)
}
static const double QP = computeQP();

// ---- малый 3-вектор (без внешних библиотек) ---------------------------------
struct Vec3 {
    double x = 0, y = 0, z = 0;
    Vec3 operator+(const Vec3& o) const { return {x + o.x, y + o.y, z + o.z}; }
    Vec3 operator-(const Vec3& o) const { return {x - o.x, y - o.y, z - o.z}; }
    Vec3 operator*(double s)      const { return {x * s, y * s, z * s}; }
    double dot(const Vec3& o)     const { return x * o.x + y * o.y + z * o.z; }
};

// ---- фазовое состояние -------------------------------------------------------
struct State { double x = 0, px = 0, y = 0, py = 0, t = 0, d = 0; };

// ---- магнит ------------------------------------------------------------------
struct Magnet {
    Vec3 center, axis, uhat, vhat, entry, exit;
    double yaw = 0, Leff = 0, L = 0, edge1 = 0, edge2 = 0;   // edge1/edge2 — краевые углы E1/E2 (.lte); физика края — Фаза 3
    double theta0 = 0, k1 = 0, k1v = 0, k2 = 0, k3 = 0, kEdge = 0;   // k1v — верт. градиент (∂By/∂v≠∂Bz/∂u); kEdge — краевой фриндж
    bool monitor = false;          // ДАТЧИК ПОЛОЖЕНИЯ ПУЧКА (BPM): пропускает частицы, регистрирует центроид
    std::string sensorId;          // = instance_id экземпляра датчика в решётке
    int nsl = 40;                  // число срезов тела (drift-kick-drift)
    bool nonzero = false;
    std::string dev, kind;         // kind — тип элемента (dipole/quad/dq/...) из lattice_configs
    std::vector<std::string> devs;   // типы устройств со-локованных обмоток
};

// ---- CSV-ридер (как в src/main.cpp): пропуск заголовка, ';', trim, EU-запятая → '.') ----
inline std::vector<std::vector<std::string>> readCSV(const std::string& path) {
    std::vector<std::vector<std::string>> data;
    std::ifstream f(path);
    if (!f.is_open()) return data;
    std::string line;
    std::getline(f, line);                       // заголовок
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        std::vector<std::string> row;
        std::stringstream ss(line);
        std::string val;
        while (std::getline(ss, val, ';')) {
            val.erase(0, val.find_first_not_of(" \t\r\n\""));
            size_t last = val.find_last_not_of(" \t\r\n\"");
            val.erase(last == std::string::npos ? 0 : last + 1);
            for (char& c : val) if (c == ',') c = '.';
            row.push_back(val);
        }
        if (!row.empty()) data.push_back(row);
    }
    return data;
}

inline double parseOr0(const std::vector<std::string>& r, size_t idx) {
    if (idx >= r.size() || r[idx].empty()) return 0.0;
    try { return std::stod(r[idx]); } catch (...) { return 0.0; }
}

// ---- локальные орты рамки по yaw --------------------------------------------
inline void frame(double yaw, Vec3& axis, Vec3& uhat, Vec3& vhat) {
    axis = { std::cos(yaw),  std::sin(yaw), 0.0};
    uhat = {-std::sin(yaw),  std::cos(yaw), 0.0};
    vhat = { 0.0, 0.0, 1.0};
}

// ---- база данных (те же таблицы) --------------------------------------------
struct DB {
    std::map<std::string, double> eff, edgeff;                  // Leff, форм-фактор края по типу (phys_length УДАЛЁН)
    std::map<std::string, std::vector<std::string>> trans;      // element_transforms строка
};

inline DB loadDB(const std::string& ws) {
    DB db;
    for (auto& r : readCSV(ws + "tables/hardware_library.csv"))
        if (r.size() > 10) {
            try { db.eff[r[0]] = std::stod(r[10]); } catch (...) {}
            if (r.size() > 11) { try { db.edgeff[r[0]] = std::stod(r[11]); } catch (...) {} }  // edge_ff ([11] после eff_length; default 0)
        }
    for (auto& r : readCSV(ws + "tables/element_transforms.csv"))
        if (r.size() > 4) db.trans[r[0]] = r;
    return db;
}

// ---- сборка магнита из группы со-локованных обмоток --------------------------
// group[i] = строка lattice_configs: [0]lattice_id [1]instance_id [2]device_type_id
//   [3]x [4]y [5]z [6]yaw [7]pitch [8]roll [9]arg_val [10]correct_coeff
//   [11]mat_B [12]mat_G [13]mat_S [14]mat_O [15]mat_Gv
//   [16]eff_length [17]kind [18]n_slices [19]edge1 [20]edge2   (структурные — ЕДИНЫЙ ИСТОЧНИК)
inline Magnet makeMagnet(const std::vector<std::vector<std::string>>& group,
                         const DB& db, const std::string& source) {
    Magnet m;
    const auto& r0 = group[0];
    m.center = { std::stod(r0[3]), std::stod(r0[4]), std::stod(r0[5]) };
    m.yaw    = std::stod(r0[6]);
    frame(m.yaw, m.axis, m.uhat, m.vhat);
    m.dev  = r0[2];
    m.sensorId = (r0.size() > 1) ? r0[1] : "";              // instance_id (для датчиков BPM)
    // ЕДИНЫЙ ИСТОЧНИК — структурные метаданные из lattice_configs (eff_length/kind/model/n_slices/edge1/edge2).
    // Fallback на hardware_library/element_transforms — только для неполной БД (те же значения).
    m.Leff = (r0.size() > 16 && !r0[16].empty()) ? parseOr0(r0, 16)
             : (db.eff.count(m.dev) ? db.eff.at(m.dev) : 0.0);   // эффективная (проектная .lte) длина
    m.L    = m.Leff;   // ОТКАЗ от металл-длины: кик работает ТОЛЬКО по L_eff
    m.entry = m.center - m.axis * (0.5 * m.L);
    m.exit  = m.center + m.axis * (0.5 * m.L);
    if (r0.size() > 17 && !r0[17].empty()) {                 // из lattice_configs
        m.kind = r0[17];
        if (r0.size() > 18 && !r0[18].empty()) { try { m.nsl = std::stoi(r0[18]); } catch (...) {} }
    } else if (db.trans.count(m.dev)) {                      // fallback: element_transforms [1]kind [3]n_slices
        const auto& tr = db.trans.at(m.dev);
        m.kind = (tr.size() > 1) ? tr[1] : "";
        if (tr.size() > 3) { try { m.nsl = std::stoi(tr[3]); } catch (...) {} }
    }
    m.monitor = (m.kind == "monitor");                       // датчик положения пучка (BPM)
    m.edge1 = (r0.size() > 19) ? parseOr0(r0, 19) : 0.0;     // краевые углы E1/E2 (.lte); физика края — Фаза 3
    m.edge2 = (r0.size() > 20) ? parseOr0(r0, 20) : 0.0;
    double B = 0, G = 0, Gv = 0, S = 0, O = 0;               // Gv — ВЕРТ. градиент (∂By/∂v), отдельно от G
    if (source != "physical")
        throw std::runtime_error("токовый режим матричного ядра архивирован — используйте source=\"physical\" "
                                 "(физкоэф. mat_* берутся из lattice_configs)");
    for (const auto& r : group) {
        m.devs.push_back(r[2]);                              // из mat_* столбцов lattice_configs
        double gcol = parseOr0(r, 12);                       // mat_G (гор. градиент)
        B += parseOr0(r, 11); G += gcol;
        double gvcol = gcol;                                 // mat_Gv (col 15); дефолт = mat_G ⇒ гармонич.
        if (r.size() > 15 && !r[15].empty()) { try { gvcol = std::stod(r[15]); } catch (...) {} }
        Gv += gvcol;
        S += parseOr0(r, 13); O += parseOr0(r, 14);
    }
    // ДИПОЛЬ: mat_G/Gv/S/O заданы ОТНОСИТЕЛЬНО поля (доли mat_B) → в абсолютные (1/f=θ₀·g, готово к излучению)
    bool isDip = db.trans.count(m.dev) && db.trans.at(m.dev).size() > 1 && db.trans.at(m.dev)[1] == "dipole";
    if (isDip) { G *= B; Gv *= B; S *= B; O *= B; }          // градиенты диполя заданы относительно поля
    // k1=k1u (гор. ∂Bz/∂u), k1v (верт. ∂By/∂v): неголоморфн. диполь Gu≠Gv ⇒ |f_x|≠|f_y|
    m.theta0 = QP * B * m.Leff; m.k1 = QP * G; m.k1v = QP * Gv; m.k2 = QP * S; m.k3 = QP * O;
    // краевой (фриндж) верт. фокус 2-го порядка: 1/f_y = FE·θ₀²/L_eff (∝поле², offset/угол-независим).
    // = −(q/p)²∫Bz²ds·v (интеграл углового члена Δv'=−(q/p)u'v∂ₓBz по частям). Диполи: θ₀≠0.
    double FE = db.edgeff.count(m.dev) ? db.edgeff.at(m.dev) : 0.0;
    m.kEdge = (m.Leff > 0.0) ? FE * m.theta0 * m.theta0 / m.Leff : 0.0;
    m.nonzero = (std::abs(m.theta0) > 0 || std::abs(m.k1) > 0 || std::abs(m.k1v) > 0 ||
                 std::abs(m.k2) > 0 || std::abs(m.k3) > 0 || std::abs(m.kEdge) > 0);
    return m;
}

// ---- сборка линейки магнитов решётки (группировка со-локованных обмоток) -----
inline std::vector<Magnet> buildLine(const std::string& ws, const DB& db, int lattice,
                                     int lo = -1, int hi = -1, const std::string& source = "physical") {
    std::map<int, std::vector<std::string>> conf;            // instance_id → строка (этой решётки)
    for (auto& r : readCSV(ws + "tables/lattice_configs.csv"))
        if (r.size() > 2) { try { if (std::stoi(r[0]) == lattice) conf[std::stoi(r[1])] = r; } catch (...) {} }
    std::vector<int> order;                                   // instance_id по порядку (connectivity)
    for (auto& r : readCSV(ws + "tables/element_connectivity.csv"))
        if (r.size() > 1) { try { if (std::stoi(r[0]) == lattice) order.push_back(std::stoi(r[1])); } catch (...) {} }
    std::sort(order.begin(), order.end());

    std::vector<Magnet> mags;
    size_t i = 0;
    while (i < order.size()) {
        const auto& r0 = conf[order[i]];
        std::vector<std::vector<std::string>> group{ r0 };
        while (i + 1 < order.size()) {
            const auto& rn = conf[order[i + 1]];
            if (rn[3] == r0[3] && rn[4] == r0[4] && rn[5] == r0[5]) { group.push_back(rn); ++i; }
            else break;
        }
        mags.push_back(makeMagnet(group, db, source));
        ++i;
    }
    if (lo >= 0) {                                            // срез [lo:hi]
        int n = (int)mags.size();
        int a = std::max(0, std::min(lo, n));
        int b = (hi < 0) ? n : std::max(a, std::min(hi, n));
        mags = std::vector<Magnet>(mags.begin() + a, mags.begin() + b);
    }
    return mags;
}

// ---- канонический дрейф (δ-корректный). false = частица потеряна (p⊥>p) ------
inline bool drift(State& s, double ell) {
    double arg = (1.0 + s.d) * (1.0 + s.d) - s.px * s.px - s.py * s.py;
    if (arg <= 0.0) return false;
    double ps = std::sqrt(arg);
    s.x += ell * s.px / ps;
    s.y += ell * s.py / ps;
    s.t += ell * (1.0 + s.d) / ps;
    return true;
}

// ---- тело магнита: N срезов drift-kick-drift. false = потеряна ---------------
inline bool body(const Magnet& m, State& s) {
    if (!m.nonzero) return drift(s, m.L);
    int N = m.nsl;
    double dsg = m.L / N, dsf = m.Leff / N;
    double th = m.theta0, k1 = m.k1, k1v = m.k1v, k2 = m.k2, k3 = m.k3, kE = m.kEdge;
    for (int i = 0; i < N; ++i) {
        if (!drift(s, dsg * 0.5)) return false;
        double x = s.x, y = s.y;
        // НОРМАЛЬНЫЕ мультиполи (Re Cₙ): изгиб θ0, квад −k1·x, секступоль, октуполь + краевая линза
        s.px += th / N - k1 * x * dsf - 0.5 * k2 * (x * x - y * y) * dsf
                - (k3 / 6.0) * (x * x * x - 3.0 * x * y * y) * dsf;
        s.py += k1v * y * dsf + k2 * x * y * dsf
                + (k3 / 6.0) * (3.0 * x * x * y - y * y * y) * dsf - (kE / N) * y;   // k1v — верт. градиент; kE — край
        if (!drift(s, dsg * 0.5)) return false;
    }
    return true;
}

// ---- модифицированный пустой промежуток: плоскость выхода A → входа B --------
// Переводит состояние в рамку B. false = не долетит / плоскость позади / p⊥>p.
inline bool modifiedDrift(State& s, const Magnet& A, const Magnet& B) {
    Vec3 P = A.exit + A.uhat * s.x + A.vhat * s.y;           // глоб. точка на выходной плоскости A
    double arg = (1.0 + s.d) * (1.0 + s.d) - s.px * s.px - s.py * s.py;
    if (arg <= 0.0) return false;
    double ps = std::sqrt(arg);
    Vec3 dir = (A.uhat * s.px + A.vhat * s.py + A.axis * ps) * (1.0 / (1.0 + s.d));
    double denom = B.axis.dot(dir);
    if (denom <= 1e-9) return false;                          // пучок не долетит до входа B
    double tt = B.axis.dot(B.entry - P) / denom;
    if (tt < 0.0) return false;                              // плоскость входа B позади частицы
    Vec3 Q = P + dir * tt;
    Vec3 rel = Q - B.entry;
    s.x  = rel.dot(B.uhat);       s.y  = rel.dot(B.vhat);
    s.px = (1.0 + s.d) * dir.dot(B.uhat);
    s.py = (1.0 + s.d) * dir.dot(B.vhat);
    s.t += tt;
    return true;
}

// ---- результат трека одной частицы ------------------------------------------
struct TrackResult {
    State fin;
    int lostAt = -1;            // индекс магнита потери, −1 если долетела
    std::string reason;
};

// Датчики: sensorId → [(x,px,y,py,t,δ) по частицам] (полный 6D-стейт: центроиды углов + σ)
using BPMAcc = std::map<std::string, std::vector<std::array<double, 6>>>;

inline TrackResult track(State s, const std::vector<Magnet>& mags, BPMAcc* monitors = nullptr) {
    TrackResult r;
    if (mags.empty()) { r.fin = s; return r; }
    if (monitors && mags[0].monitor) (*monitors)[mags[0].sensorId].push_back({s.x, s.px, s.y, s.py, s.t, s.d});
    if (!body(mags[0], s)) { r.fin = s; r.lostAt = 0; r.reason = "тело/p⊥>p"; return r; }
    for (size_t j = 1; j < mags.size(); ++j) {
        if (!modifiedDrift(s, mags[j - 1], mags[j])) { r.fin = s; r.lostAt = (int)j; r.reason = "модиф.дрейф"; return r; }
        if (monitors && mags[j].monitor) (*monitors)[mags[j].sensorId].push_back({s.x, s.px, s.y, s.py, s.t, s.d});
        if (!body(mags[j], s))                       { r.fin = s; r.lostAt = (int)j; r.reason = "тело/p⊥>p";   return r; }
    }
    r.fin = s;
    return r;
}

// ---- BPM: буфер в ОЗУ → сброс в БД ПОСЛЕ прогона; центроид+σ, turn-by-turn, модель несовершенства ----
struct BPMImperf { double dx = 0, dy = 0, tilt = 0, gx = 1, gy = 1, sx0 = 0, sy0 = 0; };  // смещ./наклон/усил./setpoint
using BPMConfig = std::map<std::string, BPMImperf>;
struct BPMRecord { std::string sid; int turn; double x, y, xp, yp, t, d, sigx, sigy; int n; };

// Несовершенства из bpm_config.csv: lattice_id;id_sensor;dx;dy;tilt;gain_x;gain_y;setpoint_x;setpoint_y
inline BPMConfig loadBPMConfig(const std::string& path, int lattice) {
    BPMConfig cfg; std::ifstream in(path); if (!in) return cfg;
    std::string line; bool first = true;
    while (std::getline(in, line)) {
        if (first) { first = false; continue; }
        if (line.empty()) continue;
        std::vector<std::string> f; std::stringstream ss(line); std::string tok;
        while (std::getline(ss, tok, ';')) f.push_back(tok);
        if (f.size() < 9) continue;
        if (!f[0].empty()) { try { if (std::stoi(f[0]) != lattice) continue; } catch (...) {} }   // фильтр по решётке
        auto num = [&](int i, double def) { try { return f[i].empty() ? def : std::stod(f[i]); } catch (...) { return def; } };
        cfg[f[1]] = { num(2,0), num(3,0), num(4,0), num(5,1), num(6,1), num(7,0), num(8,0) };
    }
    return cfg;
}

struct BPMRegistry {
    BPMConfig config;
    std::vector<BPMRecord> records;                        // буфер в ОЗУ (в БД только по flush)

    void addTurn(const BPMAcc& mon, int turn) {            // центроид + σ (RMS) + несовершенство → буфер
        for (const auto& kv : mon) {
            size_t n = kv.second.size(); if (!n) continue;
            double c[6] = {0, 0, 0, 0, 0, 0};
            for (const auto& a : kv.second) for (int k = 0; k < 6; ++k) c[k] += a[k];
            for (int k = 0; k < 6; ++k) c[k] /= (double)n;
            double vx = 0, vy = 0;
            for (const auto& a : kv.second) { vx += (a[0]-c[0])*(a[0]-c[0]); vy += (a[2]-c[2])*(a[2]-c[2]); }
            double sigx = std::sqrt(vx / n), sigy = std::sqrt(vy / n);
            BPMImperf im; auto it = config.find(kv.first); if (it != config.end()) im = it->second;
            double cx = im.gx * (c[0] - im.dx), cy = im.gy * (c[2] - im.dy);   // усиление·(центроид−смещение)
            double ct = std::cos(im.tilt), s3 = std::sin(im.tilt);            // поворот на tilt, −setpoint
            records.push_back({ kv.first, turn, cx*ct + cy*s3 - im.sx0, cy*ct - cx*s3 - im.sy0,
                                c[1], c[3], c[4], c[5], sigx, sigy, (int)n });
        }
    }

    void flush(const std::string& idExp, const std::string& path) {   // сброс буфера → bpm_readings.csv
        std::vector<std::string> keep;
        { std::ifstream in(path); std::string line; bool first = true;
          while (std::getline(in, line)) {
              if (first) { first = false; continue; }
              if (line.empty()) continue;
              size_t p = line.find(';'); std::string ex = (p == std::string::npos) ? line : line.substr(0, p);
              if (ex != idExp) keep.push_back(line);        // id_exp — 1-й столбец; чужие эксперименты сохранить
          } }
        std::ofstream out(path);
        out << "id_exp;id_sensor;turn;x;y;xp;yp;t;delta;sigma_x;sigma_y;n_particles\n";
        for (const auto& l : keep) out << l << "\n";
        for (const auto& r : records)
            out << idExp << ";" << r.sid << ";" << r.turn << ";" << r.x << ";" << r.y << ";" << r.xp << ";"
                << r.yp << ";" << r.t << ";" << r.d << ";" << r.sigx << ";" << r.sigy << ";" << r.n << "\n";
    }
};

} // namespace mtx
#endif
