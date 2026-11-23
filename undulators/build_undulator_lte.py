# -*- coding: utf-8 -*-
"""Генератор .lte-определения ОДНОГО ондулятора — линейка знакопеременных диполей.

ЗАЧЕМ. Ондулятор одинаков для всех мест установки, поэтому описывается ОДИН раз и вставляется
сборщиком (assemble_undulators.py) в нужные промежутки. Руками 2N полюсов набивать нельзя —
концевые полюса считаются по формуле на ЗАМЫКАНИЕ, поэтому устройство генерируется скриптом
(как и сам sila2026.lte — машинно-сгенерированный).

МОДЕЛЬ. N периодов → M = 2N полюсов со знакопеременным углом; каждый полупериод λ_u/2 — прямоугольный полюс-
ЭКВИВАЛЕНТ «ЖЁСТКИЙ КРАЙ» в его центре:  L' = (∫B)²/∫B² = (8/π²)·λ_u/2 (7,30 из 9 мм), поле (π/4)·B0, грани
E1 = E2 = ANGLE/2, промежутки λ_u/2 − L' (1,70 мм) между полюсами и половинки по торцам. Совпадают сразу:
  ∫B  полупериода синусоиды — угол θ_full = B0·λ_u/(π·Bρ) = 2K/γ (орбита-качка, K, замыкание);
  ∫B² — вертикальная фокусировка граней Σθ²/L' = ∫B²ds/(Bρ)² = ½(B0/Bρ)²·L и мощность СИ;
  горизонталь — дрейф: грани θ/2 гасят слабую фокусировку 1/ρ² секторного тела (прямоугольный магнит).
Отдельной линзы (бывший UFC) не нужно: фокусировка — свойство полюсов (почему не полюса вплотную с подогнанными
E: края лишь перекладывают фокусировку между плоскостями, сумма = ∫(B/Bρ)²ds, у ступеньки во весь полупериод —
81% синусоиды). Концевые полюса масштабируются на ЗАМЫКАНИЕ: I1 ∝ Σθ = 0 (выход параллельно оси) и I2 ∝ Σ i·θ = 0
(и на оси) — центры полюсов на (i+½)·λ_u/2, как раньше. Масштаб всех углов (и граней) общим множителем
k∈[0,1] («челюсти сведены→разведены») СОХРАНЯЕТ замыкание — сила ондулятора = ровно одна ручка; живой зазор в
трекинге — команда BSCALE (угол и грани полюсов × f, рамки заново), фокусировка граней ∝ f².
Концевые полюса вдвое слабее ⇒ их ∫B² меньше: суммарная вертикальная фокусировка на ~2,5% ниже ½(B0/Bρ)²·N·λ_u
(формулы «по полной длине»). См. [[undulator-parameters]], undulators/migrate_hardedge_poles.py.

ЗАПУСК
  py -3.14 undulators/build_undulator_lte.py                 # CPMU18 при зазоре 5 мм (K≈1.85)
  py -3.14 undulators/build_undulator_lte.py --gap 10        # другое поле через зазор
  py -3.14 undulators/build_undulator_lte.py --b0 1.0 --prefix U --name UND_CPMU18
"""
import os, sys, argparse, math
sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

# единый источник констант и физики поля — паспорт излучения
import undulator_sr as sr


def pole_angles(B0, p=sr.P):
    """Углы 2N полюсов [рад]: знакопеременные, концевые решены на I1=I2=0.

    Возвращает (angles, ell, theta_full, c0, c1, res_I1, res_I2), где ell — длина полюса (λ_u/2),
    c0/c1 — масштабы первого/последнего полюсов, res_* — остаточные интегралы после решения."""
    Brho = p.E / 0.299792458                       # Т·м (магнитная жёсткость при энергии пучка)
    ell = p.lam_u / 2.0
    theta_full = B0 * p.lam_u / (math.pi * Brho)    # = 2K/γ, полуинтеграл синусоиды → угол
    M = 2 * p.N

    # базовые знакопеременные углы (i=0 положительный)
    b = [((-1) ** i) * theta_full for i in range(M)]
    # концевые корректоры на полюсах 0 и M-1: решаем 2×2 на Σθ=0 и Σ i·θ=0
    #   усл. B (I2):  c1·b[M-1]·(M-1) + Σ_{1..M-2} i·b[i] = 0   (c0 не входит: i=0)
    #   усл. A (I1):  c0·b[0] + c1·b[M-1] + Σ_{1..M-2} b[i] = 0
    S1_mid = sum(b[i] for i in range(1, M - 1))
    Si_mid = sum(i * b[i] for i in range(1, M - 1))
    c1 = -Si_mid / (b[M - 1] * (M - 1))
    c0 = -(c1 * b[M - 1] + S1_mid) / b[0]
    ang = list(b)
    ang[0] = c0 * b[0]
    ang[M - 1] = c1 * b[M - 1]

    res_I1 = sum(ang)
    res_I2 = sum(i * ang[i] for i in range(M))
    return ang, ell, theta_full, c0, c1, res_I1, res_I2


HARD_EDGE = 8.0 / math.pi ** 2         # L'/(λ_u/2): эквивалент «жёсткий край» (совпадают ∫B и ∫B² полупериода)


def build_lte(B0, p=sr.P, prefix="U", name="UND_CPMU18", n_slices=1):
    """Текст .lte-определения ондулятора + отчёт-словарь."""
    ang, ell, theta_full, c0, c1, r1, r2 = pole_angles(B0, p)
    M = len(ang)
    K = sr.K_of_B0(B0, p)
    g = sr.gamma(p)
    Ldev = M * ell                                  # длина устройства = N·λ_u
    Lp = HARD_EDGE * ell                            # длина полюса-эквивалента
    dg = ell - Lp                                   # промежуток между полюсами

    def ename(i):
        if i == 0:
            return f"{prefix}E0"
        if i == M - 1:
            return f"{prefix}E1"
        return f"{prefix}PP" if ang[i] > 0 else f"{prefix}PM"

    defs, seen = [], {}
    for i in range(M):
        nm = ename(i)
        if nm not in seen:
            seen[nm] = ang[i]
            defs.append((nm, ang[i]))
    seq = [f"{prefix}DH"]                           # полупромежуток — торец устройства
    for i in range(M):
        seq += [ename(i)] + ([f"{prefix}DG"] if i < M - 1 else [f"{prefix}DH"])

    lines = [
        f"! сгенерировано build_undulator_lte.py: CPMU λ_u={p.lam_u*1e3:.0f}мм N={p.N} "
        f"L={Ldev:.4f}м B0={B0:.4f}Т K={K:.4f} γ={g:.0f}",
        f"! полюса-эквивалент «жёсткий край»: L'=(8/π²)·λ_u/2={Lp*1e3:.4f} мм, поле (π/4)·B0, грани E1=E2=ANGLE/2, "
        f"промежутки {dg*1e3:.4f} мм; концевые решены на замыкание I1=I2=0 (c0={c0:.4f}, c1={c1:.4f})",
        f"! сила ондулятора = общий множитель k∈[0,1] на ВСЕ углы и грани (замыкание сохраняется); "
        f"ручка на устройство — через assemble_undulators.py",
        f"! полный полюс θ_full={theta_full:.6e} рад = 2K/γ; фокусировка — гранями полюсов (без отдельной линзы)",
    ]
    for nm, a in defs:
        lines.append(f"{nm}: CSBEND, L={Lp:.12g}, ANGLE={a:.15e}, E1={0.5 * a:.15e}, E2={0.5 * a:.15e}, "
                     f"N_SLICES={n_slices}")
    lines.append(f"{prefix}DG: DRIF, L={dg:.12g}")
    lines.append(f"{prefix}DH: DRIF, L={0.5 * dg:.12g}")
    lines.append(f"{prefix}MK0: MARK")
    lines.append(f"{prefix}MK1: MARK")
    body = ",".join([f"{prefix}MK0"] + seq + [f"{prefix}MK1"])
    toks = body.split(",")
    rows = [",".join(toks[i:i + 12]) for i in range(0, len(toks), 12)]
    lines.append(f"{name}: LINE=(" + ",&\n".join(rows) + ")")

    # проверка эквивалента: вертикальная фокусировка граней и ∫B² против синусоиды
    Brho = p.E / 0.299792458
    ky = 0.5 * (B0 / Brho) ** 2                     # м^-2 (естественный k_y синусоиды)
    foc = sum(2.0 * math.tan(0.5 * t) * t / Lp for t in ang)          # Σ (грани) tan(E)/ρ по двум граням
    dnu_y = foc * p.beta_y / (4 * math.pi)
    rep = {
        "B0, Т": B0, "K": K, "θ_full, рад": theta_full, "полюсов (2N)": M,
        "длина устройства, м": Ldev, "полюс L', мм": Lp * 1e3, "промежуток, мм": dg * 1e3,
        "c0 (конц.)": c0, "c1 (конц.)": c1, "остаток I1=Σθ, рад": r1, "остаток I2=Σiθ": r2,
        "верт. фокус граней Σθ²/L', м^-1": foc, "½(B0/Bρ)²·L (синусоида), м^-1": ky * Ldev,
        "Δν_y/устр (β_y=%.1f доп.)" % p.beta_y: dnu_y,
    }
    return "\n".join(lines) + "\n", rep


def main():
    ap = argparse.ArgumentParser(description="сгенерировать .lte-определение одного ондулятора")
    ap.add_argument("--gap", type=float, default=None, help="зазор, мм (поле B0 через ролл-офф)")
    ap.add_argument("--b0", type=float, default=None, help="пиковое поле B0, Т (переопределяет --gap)")
    ap.add_argument("--name", default="UND_CPMU18", help="имя линии устройства")
    ap.add_argument("--prefix", default="U", help="префикс имён элементов полюсов")
    ap.add_argument("--n-slices", type=int, default=1, help="N_SLICES на полюс")
    ap.add_argument("--out", default=os.path.join(os.path.dirname(__file__), "undulator_cpmu18.lte"),
                    help="куда писать .lte")
    a = ap.parse_args()

    if a.b0 is not None:
        B0 = a.b0
    elif a.gap is not None:
        B0 = sr.B0_of_gap(a.gap, sr.P)
    else:
        B0 = sr.P.B0_ref                             # зазор по умолчанию (5 мм)

    txt, rep = build_lte(B0, sr.P, a.prefix, a.name, a.n_slices)
    for k, v in rep.items():
        print(f"  {k}: {v:.6g}" if isinstance(v, float) else f"  {k}: {v}")
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    open(a.out, "w", encoding="utf-8", newline="\n").write(txt)
    print(f"[OK] {a.out}")


if __name__ == "__main__":
    main()
