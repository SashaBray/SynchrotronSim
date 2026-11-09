# -*- coding: utf-8 -*-
"""ЭТАП B — спектр и ИНТЕРФЕРЕНЦИОННАЯ картина на ФОТОННОМ ЭКРАНЕ (когерентный интеграл поля).

Для экрана (замороженного) и его КОГЕРЕНТНОГО осветителя (привязанный ондулятор) считает:
  1) СПЕКТР в центре экрана — гребёнка гармоник (= интерференционная решётка): одноэлектронный (острый)
     и РАЗМЫТЫЙ пучком (угловая расходимость σx'/σy' + энергоразброс σ_E);
  2) 2D-КАРТУ потока на экране при заданной энергии (центральный конус на гармонике / кольцо при отстройке).
Геометрия — та же, что на этапе A: экран неподвижен (глобал), источник берётся по ФАКТической геометрии
(--ws возмущённой машины) ⇒ перекос устройства двигает/крутит картину по экрану. Несколько устройств на
один экран складываются ПО ИНТЕНСИВНОСТИ (взаимно некогерентны); полюса ОДНОГО — когерентно (внутри ядра).

Пучковые параметры (σx', σy', σ_E) — из undulator_sr.P (ПРЕДПОЛОЖЕНИЯ); β можно задать реальные из реш.27.

ЗАПУСК
  py -3.14 undulators/screen_spectrum.py --lattice 27 --screen SCR_U01 --gap 5
  py -3.14 undulators/screen_spectrum.py --lattice 27 --screen SCR_U19 --ws <папка> --harmonic 1 --detune 1.004
"""
import os, sys, csv, math, argparse
sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
os.environ.setdefault("SS_SECTOR", "1")

import numpy as np
import undulator_sr as sr
from undulator_field import trajectory, amplitude, calibrate, _dirs, apply_phase_errors
from screen_analysis import load_screens, load_lattice, load_groups, emitter, _unit, RESDIR


def beam_params(p, beta_x, beta_y, kappa, sigE, etaxp=0.0):
    """Расходимости пучка σx',σy' (с дисперсионной добавкой η'x·σ_E в горизонтали)."""
    ey = kappa * p.emit_x
    sxp = math.sqrt(p.emit_x / beta_x + (etaxp * sigE) ** 2)
    return dict(sxp=sxp, syp=math.sqrt(ey / beta_y), sigE=sigE, bx=beta_x, by=beta_y, ey=ey)


def energy_smear(E, S, rel_sigma):
    """Свёртка спектра по ОТНОСИТЕЛЬНОЙ энергии (энергоразброс: δE/E=2σ_E постоянен) через ln-сетку."""
    from scipy.ndimage import gaussian_filter1d
    if rel_sigma <= 0:
        return S
    u = np.log(E); uu = np.linspace(u[0], u[-1], u.size)
    Su = np.interp(uu, u, S)
    sig_pix = rel_sigma / (uu[1] - uu[0])
    Sf = gaussian_filter1d(Su, sig_pix, mode="nearest")
    return np.interp(u, uu, Sf)


def center_spectrum(tr, scale, thx0, thy0, Es, beam, nang=5, nsig=2.5):
    """Спектр в точке (thx0,thy0): одноэлектронный и размытый пучком."""
    I_se = amplitude(tr, Es, _dirs(thx0, thy0)[None, :])[0] * scale
    if beam is None:
        return I_se, None
    offs = np.linspace(-nsig, nsig, nang)
    acc = np.zeros_like(Es); W = 0.0
    for ax in offs:
        for ay in offs:
            w = math.exp(-(ax ** 2 + ay ** 2) / 2)
            acc += w * amplitude(tr, Es, _dirs(thx0 - ax * beam["sxp"], thy0 - ay * beam["syp"])[None, :])[0]
            W += w
    I_beam = energy_smear(Es, acc / W * scale, 2 * beam["sigE"])
    return I_se, I_beam


def screen_map(tr, scale, P, triad, dc, D, E0, mesh, window, beam):
    """2D карта потока [ф/с/мрад²/0.1%BW] в УГЛОВОМ ОКНЕ ±window (рад) вокруг центра экрана.

    Интерференционная структура (конус/кольца) живёт в центральном конусе ~σ_r' ≪ размера экрана, поэтому
    карта зумирует туда (± несколько σ_r'); оси — смещение по экрану в мм (= угол·D). Окно центрировано на
    НАПРАВЛЕНИИ на центр экрана dc (при перекосе смещается)."""
    axis, uhat, vhat = triad
    thx0 = float(np.dot(dc, uhat)); thy0 = float(np.dot(dc, vhat))
    tx = np.linspace(thx0 - window, thx0 + window, mesh)
    ty = np.linspace(thy0 - window, thy0 + window, mesh)
    TX, TY = np.meshgrid(tx, ty, indexing="ij")
    dirs = _dirs(TX, TY)                                            # (m,m,3)
    I = amplitude(tr, E0, dirs)[..., 0] * scale
    if beam is not None:
        from scipy.ndimage import gaussian_filter
        dth = tx[1] - tx[0]
        I = gaussian_filter(I, sigma=(beam["sxp"] / dth, beam["syp"] / dth), mode="nearest")
    ext = window * D * 1e3                                          # полуразмер окна на экране, мм
    return dict(I=I, extent=(-ext, ext, -ext, ext), peak=float(I.max()), D=D,
                thx0=thx0, thy0=thy0, window=window)


def find_harmonics(tr, scale, K, nmax=7):
    rows = []
    for n in range(1, nmax + 1, 2):
        En = sr.harmonic_energy_keV(n, K, 0.0)
        Es = np.linspace(0.97 * En, 1.03 * En, 200)
        s = center_spectrum(tr, scale, 0.0, 0.0, Es, None)[0]
        j = int(np.argmax(s)); rows.append((n, float(Es[j]), float(s[j])))
    return rows


def plot_spectrum(Es, Ise, Ibeam, scr, harm, path, gap):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6.2, 3.8))
    ax.semilogy(Es, np.maximum(Ise, 1), color="#9fb4c8", lw=0.9, label="один электрон (острый)")
    if Ibeam is not None:
        ax.semilogy(Es, np.maximum(Ibeam, 1), color="#185fa5", lw=1.6, label="пучок (σx',σy',σ_E)")
    for n, En, Fn in harm:
        ax.annotate(f"n={n}", (En, Fn), fontsize=8, ha="center", va="bottom", color="#333")
    ax.set_xlabel("энергия фотона, кэВ"); ax.set_ylabel("ф/с/мрад²/0.1%BW")
    ax.set_title(f"{scr['sid']} ← {scr['src']}: спектр на оси (зазор {gap:g} мм)", fontsize=9)
    ax.legend(fontsize=8); ax.grid(alpha=0.25, which="both")
    fig.tight_layout(); fig.savefig(path, dpi=120); plt.close(fig)


def plot_map(mp, scr, E0, path, gap):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(4.4, 3.8))
    im = ax.imshow(mp["I"].T, origin="lower", extent=mp["extent"], aspect="equal", cmap="inferno")
    ax.set_xlabel("u − центр конуса, мм"); ax.set_ylabel("v − центр конуса, мм")
    ax.set_title(f"{scr['sid']} ← {scr['src']}: поток при {E0:.2f} кэВ (зазор {gap:g} мм, D={mp['D']:.0f} м)\n"
                 f"пик {mp['peak']:.2e} ф/с/мрад²/0.1%BW", fontsize=8.5)
    cb = fig.colorbar(im, ax=ax, shrink=0.85); cb.set_label("ф/с/мрад²/0.1%", fontsize=8)
    fig.tight_layout(); fig.savefig(path, dpi=120); plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description="спектр и интерференционная картина на экране (этап B)")
    ap.add_argument("--lattice", type=int, default=27)
    ap.add_argument("--screen", required=True, help="screen_id из screens.csv")
    ap.add_argument("--ws", default=None, help="папка возмущённой машины (иначе номинал)")
    ap.add_argument("--gap", type=float, default=5.0)
    ap.add_argument("--harmonic", type=int, default=1, help="гармоника для 2D-карты")
    ap.add_argument("--detune", type=float, default=1.0, help="множитель энергии карты (>1 — кольцо)")
    ap.add_argument("--beta-x", type=float, default=None, help="β_x, м (по умолч. РЕАЛЬНАЯ из twiss реш.)")
    ap.add_argument("--beta-y", type=float, default=None, help="β_y, м (по умолч. РЕАЛЬНАЯ из twiss реш.)")
    ap.add_argument("--kappa", type=float, default=sr.P.kappa)
    ap.add_argument("--sigE", type=float, default=sr.P.sig_E)
    ap.add_argument("--ppp", type=int, default=64)
    ap.add_argument("--mesh", type=int, default=121)
    ap.add_argument("--window-urad", type=float, default=None, help="полуокно карты, мкрад (по умолч. ~8·σ_r')")
    ap.add_argument("--phase-error-deg", type=float, default=0.0, help="RMS фазовая ошибка ИД, град (качество)")
    ap.add_argument("--phase-seed", type=int, default=0)
    ap.add_argument("--no-beam", action="store_true", help="без размытия пучком (чистый один электрон)")
    a = ap.parse_args()

    scr = next((s for s in load_screens(a.lattice) if s["sid"] == a.screen), None)
    if scr is None:
        raise SystemExit(f"экран {a.screen} не найден в screens.csv решётки {a.lattice}")
    groups = load_groups(a.lattice); ci, byid = load_lattice(a.lattice, a.ws)
    insts = groups.get(scr["src"])
    P, axis = emitter(insts, ci, byid)
    uhat = _unit(np.cross([0, 0, 1.0], axis)); vhat = _unit(np.cross(axis, uhat))

    B0 = sr.B0_of_gap(a.gap); K = sr.K_of_B0(B0); g = sr.gamma()
    tr = trajectory(K, sr.P.lam_u, sr.P.N, g, a.ppp)
    scale, E1 = calibrate(tr)                                   # калибровка по ИДЕАЛЬНОЙ траектории
    if a.phase_error_deg:                                       # качество ИД: фазовые ошибки по периодам
        tr = apply_phase_errors(tr, a.phase_error_deg, seed=a.phase_seed)

    # РЕАЛЬНЫЕ β,η из twiss реш. в месте ондулятора (если не заданы явно)
    etaxp = 0.0; bsrc = "задано"
    if (a.beta_x is None or a.beta_y is None) and not a.no_beam:
        from lattice_optics import undulator_optics
        o = undulator_optics(a.lattice, scr["src"])
        if a.beta_x is None: a.beta_x = o["betax"]
        if a.beta_y is None: a.beta_y = o["betay"]
        etaxp = o["etaxp"]; bsrc = f"twiss реш.{a.lattice} (s={o['s']:.1f} м)"
    beam = None if a.no_beam else beam_params(sr.P, a.beta_x, a.beta_y, a.kappa, a.sigE, etaxp)

    # направление на центр экрана в рамке устройства (перекос ⇒ ненулевые углы)
    dc = _unit(scr["C"] - P); thx0 = float(np.dot(dc, uhat)); thy0 = float(np.dot(dc, vhat))
    tag = "номинал" if not a.ws else os.path.basename(a.ws.rstrip("/\\"))
    print(f"решётка {a.lattice} [{tag}] экран {a.screen} ← {scr['src']}: зазор {a.gap:g} мм, "
          f"B0={B0:.3f} Т, K={K:.3f}; центр экрана в рамке устройства θ=({thx0*1e6:+.1f},{thy0*1e6:+.1f}) мкрад")
    if beam:
        print(f"пучок: σx'={beam['sxp']*1e6:.2f} мкрад, σy'={beam['syp']*1e6:.2f} мкрад "
              f"(β=({a.beta_x:.3f},{a.beta_y:.3f}) м [{bsrc}], εx={sr.P.emit_x*1e12:.0f} пм, "
              f"εy={beam['ey']*1e12:.2f} пм, связь {a.kappa*100:.0f}% [предпол.]), σ_E={a.sigE*100:.2f}% [предпол.]")

    # спектр на оси (в точке центра экрана)
    Emax = sr.harmonic_energy_keV(9, K, 0.0) * 1.05
    Es = np.linspace(0.5 * E1, Emax, 3000)
    Ise, Ibeam = center_spectrum(tr, scale, thx0, thy0, Es, beam)
    harm = find_harmonics(tr, scale, K)
    os.makedirs(RESDIR, exist_ok=True)
    sfx = "" if not a.ws else "_pert"
    sp_png = os.path.join(RESDIR, f"L{a.lattice}_{a.screen}_spectrum{sfx}.png")
    plot_spectrum(Es, Ise, Ibeam, scr, harm, sp_png, a.gap)

    print("\n  гармоники на ОСИ УСТРОЙСТВА (из интеграла, независимо):")
    for n, En, Fn in harm:
        print(f"    n={n}: E={En:7.3f} кэВ (sr {sr.harmonic_energy_keV(n,K,0.0):7.3f})  поток {Fn:.3e} ф/с/мрад²/0.1%")
    # 1-я гармоника В ТОЧКЕ ЦЕНТРА ЭКРАНА (перекос уводит конус с точки → красное смещение и падение потока)
    band = Es < 1.12 * E1
    j = int(np.argmax(Ise[band])); Epk = float(Es[band][j]); Fpk = float(Ise[band][j])
    jb = int(np.argmax(Ibeam[band])) if Ibeam is not None else j
    Fb = float(Ibeam[band][jb]) if Ibeam is not None else float("nan")
    print(f"  1-я гарм. В ЦЕНТРЕ ЭКРАНА (θ={thx0*1e6:+.0f},{thy0*1e6:+.0f} мкрад): "
          f"E={Epk:.3f} кэВ, поток 1эл={Fpk:.3e}" + (f", пучок={Fb:.3e}" if Ibeam is not None else "")
          + "  ← падает/смещается при перекосе")

    # 2D-карта: угловое окно вокруг центрального конуса
    E0 = sr.harmonic_energy_keV(a.harmonic, K, 0.0) * a.detune
    srp = sr.sigma_photon(E0)[1]                                    # σ_r' фотонного конуса
    window = a.window_urad * 1e-6 if a.window_urad else 8.0 * srp
    D = float(np.linalg.norm(scr["C"] - P))
    mp = screen_map(tr, scale, P, (axis, uhat, vhat), dc, D, np.array([E0]), a.mesh, window, beam)
    mp_png = os.path.join(RESDIR, f"L{a.lattice}_{a.screen}_map_{E0:.2f}keV{sfx}.png")
    plot_map(mp, scr, E0, mp_png, a.gap)
    print(f"\n  2D-карта при {E0:.2f} кэВ (гарм.{a.harmonic}×{a.detune:g}): пик {mp['peak']:.3e}; "
          f"окно ±{window*1e6:.0f} мкрад (±{window*D*1e3:.2f} мм на экране, D={D:.1f} м); "
          f"центр конуса θ=({mp['thx0']*1e6:+.1f},{mp['thy0']*1e6:+.1f}) мкрад")
    print(f"\n[OK] {sp_png}\n[OK] {mp_png}")


if __name__ == "__main__":
    main()
