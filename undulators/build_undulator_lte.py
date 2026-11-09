# -*- coding: utf-8 -*-
"""Генератор .lte-определения ОДНОГО ондулятора — линейка знакопеременных диполей.

ЗАЧЕМ. Ондулятор одинаков для всех мест установки, поэтому описывается ОДИН раз и вставляется
сборщиком (assemble_undulators.py) в нужные промежутки. Руками 2N полюсов набивать нельзя —
концевые полюса считаются по формуле на ЗАМЫКАНИЕ, поэтому устройство генерируется скриптом
(как и сам sila2026.lte — машинно-сгенерированный).

МОДЕЛЬ. N периодов → M = 2N полюсов, каждый — короткий сектор CSBEND длиной λ_u/2 со знакопеременным
углом. Полный (внутренний) полюс: θ_full = B0·λ_u/(π·Bρ) = 2·K/γ (полуинтеграл синусоиды, переведённый
в угол отклонения). Концевые полюса масштабируются так, чтобы занулить ПЕРВЫЙ и ВТОРОЙ интегралы поля
(I1 ∝ Σθ = 0 — пучок выходит параллельно оси; I2 ∝ Σ i·θ = 0 — и на оси), т.е. устройство не даёт
ни кика, ни смещения замкнутой орбиты. Масштаб всех углов общим множителем k∈[0,1] («челюсти сведены→
разведены») СОХРАНЯЕТ замыкание (линейно по θ) — поэтому сила ондулятора = ровно одна ручка.

ЧЕГО НЕ ДЕЛАЕТ (осознанно). Плоская линейка секторных диполей воспроизводит ОРБИТУ (горизонтальную
качку), ДЛИНУ ПУТИ и орбитное ВОЗМУЩЕНИЕ при сканировании зазора, но НЕ естественную вертикальную
фокусировку ондулятора k_y=½(B0/Bρ)² (эффект 2-го порядка от y-зависимости поля, его в теле CSBEND
нет). Скрипт ПЕЧАТАЕТ целевые k_y·L и Δν_y, чтобы было видно, что именно не покрыто телом; способ
представления (отдельная вигглерная линза / распределённые тонкие линзы / краевые углы) калибруется
отдельно под Δν_y. См. [[undulator-parameters]].

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


def build_lte(B0, p=sr.P, prefix="U", name="UND_CPMU18", n_slices=1):
    """Текст .lte-определения ондулятора + отчёт-словарь."""
    ang, ell, theta_full, c0, c1, r1, r2 = pole_angles(B0, p)
    M = len(ang)
    K = sr.K_of_B0(B0, p)
    g = sr.gamma(p)
    Ldev = M * ell                                  # магнитная длина = N·λ_u

    # уникальные типы полюсов: внутренние ±, и отдельные концевые (их углы отличаются)
    # большинство полюсов — ровно ±theta_full; концевые 0 и M-1 — свои. Соберём компактно.
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
    seq = [ename(i) for i in range(M)]

    lines = [
        f"! сгенерировано build_undulator_lte.py: CPMU λ_u={p.lam_u*1e3:.0f}мм N={p.N} "
        f"L={Ldev:.4f}м B0={B0:.4f}Т K={K:.4f} γ={g:.0f}",
        f"! линейка знакопеременных секторных диполей; концевые решены на замыкание "
        f"I1=I2=0 (c0={c0:.4f}, c1={c1:.4f})",
        f"! сила ондулятора = общий множитель k∈[0,1] на ВСЕ углы (замыкание сохраняется); "
        f"ручка на устройство — через assemble_undulators.py",
        f"! полный полюс θ_full={theta_full:.6e} рад = 2K/γ; длина полюса λ_u/2={ell:.4f} м",
    ]
    for nm, a in defs:
        lines.append(f"{nm}: CSBEND, L={ell:.12g}, ANGLE={a:.15e}, N_SLICES={n_slices}")
    # линия устройства: вход-маркер + полюса + выход-маркер
    lines.append(f"{prefix}MK0: MARK")
    lines.append(f"{prefix}MK1: MARK")
    body = ",".join([f"{prefix}MK0"] + seq + [f"{prefix}MK1"])
    # перенос по 12 токенов (стиль проекта)
    toks = body.split(",")
    rows = [",".join(toks[i:i + 12]) for i in range(0, len(toks), 12)]
    lines.append(f"{name}: LINE=(" + ",&\n".join(rows) + ")")

    # целевая естественная вертикальная фокусировка (НЕ в теле CSBEND — для информации)
    Brho = p.E / 0.299792458
    ky = 0.5 * (B0 / Brho) ** 2                      # м^-2 (= Ky фокусировки, верт.)
    kyL = ky * Ldev
    dnu_y = kyL * p.beta_y / (4 * math.pi)           # Δν_y ≈ β_y·∫k ds /4π (β_y=допущение; реальная из оптики)
    # НОСИТЕЛЬ фокусировки в киковом ядре: расщеплённый верт. градиент k1v (НЕ edge_ff — та линза
    # не применяется в секторном теле). Ky=−k1v ⇒ k1v=−k_y (одинаков по полюсам, ∝B0² ⇒ ∝k² при зазоре).
    # В lattice_configs: kind=dq (абсолютный градиент), mat_G=0, mat_Gv=k1v/QP. Ставится apply_focus.
    k1v = -ky
    rep = {
        "B0, Т": B0, "K": K, "θ_full, рад": theta_full, "полюсов (2N)": M,
        "длина устройства, м": Ldev, "c0 (конц.)": c0, "c1 (конц.)": c1,
        "остаток I1=Σθ, рад": r1, "остаток I2=Σiθ": r2,
        "естеств. k_y, м^-2": ky, "∫k_y ds, м^-1": kyL, "k1v на полюс, м^-2": k1v,
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
