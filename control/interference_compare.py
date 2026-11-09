# -*- coding: utf-8 -*-
"""Стадия B: интерференционная картина излучения ондулятора при СУ ВКЛ vs ВЫКЛ.

Физика: в дальней зоне поперечный УГОЛ пучка x′ жёстко сдвигает всю интерференционную картину по экрану
на x′·D. Пучок, блуждающий от возмущений (СУ выкл), за время накопления детектора РАЗМАЗЫВАЕТ картину;
орбитный FOFB (СУ вкл) гасит это блуждание → картина резче. Поэтому:
  • чистую когерентную картину одного ондулятора считаем ОДИН раз (screen_spectrum.screen_map, этап B);
  • угловой джиттер пучка в источнике снимаем из feedback-сима (control.tbt_plant) ВЫКЛ и ВКЛ: photon.x/y
    = центроид пятна на экране (x+D·x′) по тактам; σ_угл = СКО(пятна)/D;
  • картину сворачиваем с гауссовым ядром этого джиттера (через штатное размытие screen_map: складываем
    джиттер в квадратуре к собственной расходимости пучка σx′/σy′). Узкое ядро (ВКЛ) → резко; широкое (ВЫКЛ).

Картина берётся на оси (идеальный источник), параметры ондулятора — undulator_sr.P (ПРЕДПОЛОЖЕНИЯ, CPMU-18);
джиттер — с реальной машины реш.26 (контрольная платформа). Рендер — 3 карты (идеал / FOFB-вкл / FOFB-выкл)
+ метрика резкости в один HTML (screen_report).

ЗАПУСК: py -3.14 control/interference_compare.py [smoke] [--no-open]
"""
import argparse
import math
import os
import sys

os.environ.setdefault("SS_SECTOR", "1")
sys.stdout.reconfigure(encoding="utf-8")
_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, _ROOT)
sys.path.insert(0, os.path.join(_ROOT, "undulators"))

import numpy as np

import undulator_sr as sr
from undulator_field import trajectory, calibrate
from screen_spectrum import beam_params, screen_map
import screen_report as rep

from control.ring_model import RingModel

GAP_MM = 5.0
D_SCREEN = 15.0                                               # плечо источник→экран, м
BETA_X, BETA_Y = 4.0, 2.0                                     # β прямолинейного промежутка, м (ПРЕДПОЛОЖ.)
DIST_AMP = 5e-6                                               # СКО сдвига квадруполей в ансамбле, м (масштаб блюра)
N_ENS = 200                                                   # реализаций возмущения в ансамбле
RCOND = 1e-2


def source_jitter(n_ens=N_ENS, amp=DIST_AMP, seed=0):
    """Угловой разброс пучка σx′,σy′ в источнике ондулятора ПРИ СУ ВЫКЛ и ВКЛ — по линейной модели орбиты
    реш.26 (детерминированно, без шума трекинга): ансамбль сдвигов квадруполей (feed-down), орбитный FOFB
    зануляет орбиту в мониторах (SVD), источник едет следом. Возврат (σx′,σy′) выкл, (σx′,σy′) вкл."""
    rm = RingModel(lattice=26, nsl=8, n_photon=1)
    src = int(rm.photon_rows[0])
    ax, ay = f"plane.{src}.px", f"plane.{src}.py"             # угол пучка в источнике (x′,y′)
    Hc_m = rm.response(rm.hk_ids, ["ebpm.x"])["ebpm.x"]       # корректор hkick → орбита мониторов (для FOFB)
    Hcy_m = rm.response(rm.vk_ids, ["ebpm.y"])["ebpm.y"]
    Hc_s = rm.response(rm.hk_ids, [ax])[ax][0]               # hkick → угол источника x′ (строка)
    Hcy_s = rm.response(rm.vk_ids, [ay])[ay][0]
    Gx_m = rm.response(rm.qx_ids, ["ebpm.x"])["ebpm.x"]       # сдвиг квадруполя → орбита мониторов
    Gy_m = rm.response(rm.qy_ids, ["ebpm.y"])["ebpm.y"]
    Gx_s = rm.response(rm.qx_ids, [ax])[ax][0]               # сдвиг квадруполя → угол источника x′
    Gy_s = rm.response(rm.qy_ids, [ay])[ay][0]
    rm.close()
    Px = np.linalg.pinv(Hc_m, rcond=RCOND); Py = np.linalg.pinv(Hcy_m, rcond=RCOND)
    rng = np.random.default_rng(seed); nq = Gx_s.size
    off = np.zeros((n_ens, 2)); on = np.zeros((n_ens, 2))
    for k in range(n_ens):
        dx = rng.standard_normal(nq) * amp; dy = rng.standard_normal(nq) * amp
        off[k] = (Gx_s @ dx, Gy_s @ dy)                      # угол источника без коррекции
        ux = -Px @ (Gx_m @ dx); uy = -Py @ (Gy_m @ dy)      # FOFB: зануляем орбиту мониторов
        on[k] = (Gx_s @ dx + Hc_s @ ux, Gy_s @ dy + Hcy_s @ uy)   # остаточный угол источника
    return (float(np.std(off[:, 0])), float(np.std(off[:, 1])),
            float(np.std(on[:, 0])), float(np.std(on[:, 1])))


def main():
    ap = argparse.ArgumentParser(description="интерференция ондулятора: СУ вкл vs выкл")
    ap.add_argument("mode", nargs="?", default="")
    ap.add_argument("--out", default=os.path.join("workspace", "results", "interference_fb.html"))
    ap.add_argument("--no-open", action="store_true")
    a = ap.parse_args()
    smoke = a.mode == "smoke"
    n_ens = 80 if smoke else 200
    mesh = 101 if smoke else 161

    print("=== ИНТЕРФЕРЕНЦИЯ ОНДУЛЯТОРА: СУ вкл vs выкл ===")
    # --- чистая когерентная картина одного ондулятора (1-я гармоника, на оси) ---
    B0 = sr.B0_of_gap(GAP_MM); K = sr.K_of_B0(B0); g = sr.gamma()
    tr = trajectory(K, sr.P.lam_u, sr.P.N, g, 64)
    scale, E1 = calibrate(tr)
    E0 = sr.harmonic_energy_keV(1, K, 0.0)
    srp = float(sr.sigma_photon(E0)[1])                       # σ_r′ фотонного конуса
    beam0 = beam_params(sr.P, BETA_X, BETA_Y, sr.P.kappa, sr.P.sig_E)
    print(f"ондулятор CPMU-18: зазор {GAP_MM} мм, K={K:.3f}, 1-я гарм {E0:.3f} кэВ, σ_r′={srp*1e6:.2f} мкрад")
    print(f"собственная расходимость пучка σx′={beam0['sxp']*1e6:.2f}, σy′={beam0['syp']*1e6:.2f} мкрад "
          f"(β=({BETA_X},{BETA_Y}) м предпол.)\n", flush=True)

    # --- угловой разброс орбиты в источнике: ВЫКЛ и ВКЛ (линейная модель, ансамбль) ---
    D = D_SCREEN
    print(f"ансамбль {n_ens} реализаций возмущения (сдвиг квадруполей СКО {DIST_AMP*1e6:.0f} мкм); "
          f"орбитный FOFB по SVD …", flush=True)
    jx_off, jy_off, jx_on, jy_on = source_jitter(n_ens=n_ens, amp=DIST_AMP)
    print(f"  угол источника СУ ВЫКЛ: σx′={jx_off*1e6:.3f}, σy′={jy_off*1e6:.3f} мкрад")
    print(f"  угол источника СУ ВКЛ:  σx′={jx_on*1e6:.3f}, σy′={jy_on*1e6:.3f} мкрад")
    print(f"  подавление орбиты в источнике FOFB: X ×{jx_off/max(jx_on,1e-30):.1f}, Y ×{jy_off/max(jy_on,1e-30):.1f}\n", flush=True)

    # --- три картины: идеал (без джиттера) / FOFB-вкл / FOFB-выкл ---
    window = max(8.0 * srp, 4.0 * math.hypot(jx_off, jy_off))  # окно вмещает блюр СУ-выкл
    axis = np.array([0.0, 0.0, 1.0]); uhat = np.array([1.0, 0.0, 0.0]); vhat = np.array([0.0, 1.0, 0.0])
    P = np.zeros(3); dc = axis

    def pattern(jx, jy):
        b = dict(beam0)
        b["sxp"] = math.hypot(beam0["sxp"], jx); b["syp"] = math.hypot(beam0["syp"], jy)
        return screen_map(tr, scale, P, (axis, uhat, vhat), dc, D, np.array([E0]), mesh, window, b)

    mp_ideal = pattern(0.0, 0.0)
    mp_on = pattern(jx_on, jy_on)
    mp_off = pattern(jx_off, jy_off)
    pk_i, pk_on, pk_off = mp_ideal["peak"], mp_on["peak"], mp_off["peak"]
    print(f"пик потока (резкость): идеал {pk_i:.3e} | FOFB-вкл {pk_on:.3e} ({pk_on/pk_i*100:.0f}% идеала) | "
          f"FOFB-выкл {pk_off:.3e} ({pk_off/pk_i*100:.0f}% идеала)")
    print(f"восстановление резкости при включении СУ: ×{pk_on/max(pk_off,1e-30):.2f} по пику")

    # --- HTML: три карты + KPI ---
    def kpi(label, val, unit):
        return f'<div class="kpi"><div class="l">{label}</div><div class="v">{val}</div><div class="u">{unit}</div></div>'
    zl = "ф/с/мрад²/0.1%BW"
    panels = []
    for title, mp in (("Идеал (стабильный пучок)", mp_ideal), ("СУ ВКЛ (FOFB)", mp_on), ("СУ ВЫКЛ", mp_off)):
        panels.append(f'<h2>{title}</h2><div class="pic">{rep.fig_div(rep.heatmap(mp["I"][..., 0] if mp["I"].ndim==3 else mp["I"], mp["extent"], zl))}'
                      f'<div class="cap">пик {mp["peak"]:.3e} {zl}; окно ±{window*D*1e3:.2f} мм на экране (D={D:.0f} м)</div></div>')
    kpis = "".join([
        kpi("джиттер σx′ ВЫКЛ", f"{jx_off*1e6:.2f}", "мкрад"),
        kpi("джиттер σx′ ВКЛ", f"{jx_on*1e6:.2f}", "мкрад"),
        kpi("подавление X", f"×{jx_off/max(jx_on,1e-30):.1f}", "FOFB"),
        kpi("пик ВКЛ/идеал", f"{pk_on/pk_i*100:.0f}%", "резкость"),
        kpi("пик ВЫКЛ/идеал", f"{pk_off/pk_i*100:.0f}%", "резкость"),
        kpi("резкость ВКЛ/ВЫКЛ", f"×{pk_on/max(pk_off,1e-30):.2f}", "по пику"),
    ])
    body = (f'<div class="kpis">{kpis}</div>'
            f'<div class="unitbox">Картина — 1-я гармоника CPMU-18 на оси; размытие = собственная расходимость '
            f'пучка ⊕ разброс угла пучка в источнике. Разброс — из линейной модели орбиты реш.26: ансамбль '
            f'сдвигов квадруполей, орбитный FOFB (SVD) зануляет орбиту в мониторах, источник едет следом. '
            f'СУ держит орбиту → интерференционный конус резче. Параметры ондулятора — предположительные; '
            f'амплитуда возмущения подобрана под видимый блюр (важен коэффициент подавления, не абсолют).</div>'
            + "".join(panels))
    html = (rep._TPL.replace("/*TITLE*/", "Интерференция ондулятора: СУ вкл vs выкл")
            .replace("/*BODY*/", body)
            .replace("/*SUB*/", f"CPMU-18, зазор {GAP_MM} мм, {E0:.2f} кэВ · орбита из линейной модели реш.26 · "
                                f"графики интерактивны")
            .replace("/*PLOTLYJS*/", rep.get_plotlyjs()))
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    open(a.out, "w", encoding="utf-8").write(html)
    print(f"\n[OK] {a.out}")
    if not a.no_open:
        try:
            import webbrowser
            webbrowser.open("file://" + os.path.abspath(a.out))
        except Exception:
            pass


if __name__ == "__main__":
    main()
