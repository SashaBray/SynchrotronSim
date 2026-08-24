# -*- coding: utf-8 -*-
"""field_quality.py — оценка качества полевой карты магнита (для проверки экспортов ANSYS).

Принимает карту (.bin проекта или .fld ANSYS) и считает физические индикаторы, по которым
видно, насколько хорошо сошлось FEM-решение и где «пол ошибки»:

  1. ∇·B / |∇B|       — дивергенция ОБЯЗАНА быть 0 у любого поля. Чистый индикатор ошибки.
                         Меряется на нескольких радиусах: падает с r => шум; держится => невязка.
  2. ∇×B / |∇B|       — ротор обязан быть 0 в области без токов (тело магнита). Меряется в теле.
  3. поле на оси       — у мультиполя (N≥2) на оси поле = 0. Остаток = прямая мера ошибки.
  4. симметрия Bx(0)   — если карта симметрична (x≥0+зеркало), на плоскости x=mnx должно быть Bx=0.
                         Скачок = 2·Bx(0). Ненулевой => асимметрия / несошедшееся решение.
  5. чистота спектра   — доля главной мультипольной гармоники vs «неидеальные» гармоники.

Использование:
  from field_quality import assess_field
  m = assess_field("workspace/fields_bin/SP_I1_1A.bin")     # печатает отчёт, возвращает dict
  python field_quality.py path/to/field.bin [order]          # из командной строки
"""
import struct, os, re, math, sys
import numpy as np

HDR = 96


# ----------------------------- загрузка карты -----------------------------
def load(path):
    """Чтение .bin (проект) или .fld (ANSYS, текст). Возврат (header, B[nx,ny,nz,3])."""
    if path.lower().endswith(".bin"):
        with open(path, "rb") as f:
            b = f.read(HDR)
            mnx, mny, mnz, mxx, mxy, mxz, sx, sy, sz = struct.unpack_from("<9d", b, 0)
            nx, ny, nz = struct.unpack_from("<3i", b, 72)
            B = np.frombuffer(f.read(nx * ny * nz * 24), dtype="<f8").reshape(nx, ny, nz, 3)
    else:  # .fld
        with open(path, "r", encoding="utf-8", errors="ignore") as f:   # явная кодировка: на Linux дефолт иной
            hdr = f.readline()
        nums = [float(x) for x in re.findall(r"-?\d+\.?\d*", hdr)]
        # заголовок .fld в ММ, а данные (столбцы X,Y,Z) в МЕТРАХ — приводим заголовок к метрам
        mnx, mny, mnz, mxx, mxy, mxz, sx, sy, sz = [v / 1000.0 for v in nums[:9]]
        nx = round((mxx - mnx) / sx) + 1; ny = round((mxy - mny) / sy) + 1; nz = round((mxz - mnz) / sz) + 1
        B = np.loadtxt(path, skiprows=2, usecols=(3, 4, 5)).reshape(nx, ny, nz, 3)
    return dict(mnx=mnx, mny=mny, mnz=mnz, mxx=mxx, sx=sx, sy=sy, sz=sz, nx=nx, ny=ny, nz=nz), B


def _body_slice(h, B):
    """Срез тела магнита — максимум ВНЕОСЕВого поперечного поля."""
    j0, k0 = h["ny"] // 2, h["nz"] // 2
    R = max(1, round(3e-3 / h["sy"]))
    reg = B[:, j0 - R:j0 + R + 1, k0 - R:k0 + R + 1, 1:3]
    mag = np.sqrt(np.mean(reg ** 2, axis=(1, 2, 3)))
    return int(np.argmax(mag)) if mag.max() > 0 else h["nx"] // 2


def _circle(h, B, ix, r, nth=64):
    """By, Bz на окружности радиуса r (м) в срезе ix — билинейно. Возврат (By[θ], Bz[θ])."""
    j0, k0 = h["ny"] // 2, h["nz"] // 2
    by, bz = [], []
    for t in np.linspace(0, 2 * math.pi, nth, endpoint=False):
        fy = r * math.cos(t) / h["sy"] + j0; fz = r * math.sin(t) / h["sz"] + k0
        j, k = int(fy), int(fz); ty, tz = fy - j, fz - k
        vals = []
        for c in (1, 2):
            g = B[ix, :, :, c]
            vals.append(g[j, k] * (1 - ty) * (1 - tz) + g[j + 1, k] * ty * (1 - tz) +
                        g[j, k + 1] * (1 - ty) * tz + g[j + 1, k + 1] * ty * tz)
        by.append(vals[0]); bz.append(vals[1])
    return np.array(by), np.array(bz)


def detect_order(h, B, ix):
    """Порядок мультиполя по РАДИАЛЬНОМУ масштабу |B⊥| ∝ r^(N−1): диполь N=1, квадруполь 2, …"""
    def bperp(r):
        by, bz = _circle(h, B, ix, r)
        return math.sqrt(np.mean(by ** 2 + bz ** 2))
    b1, b2 = bperp(2e-3), bperp(8e-3)
    if b1 < 1e-12:
        return 2
    ratio = b2 / b1                          # = (8/2)^(N−1) = 4^(N−1)
    N = 1 + round(math.log(ratio) / math.log(4.0)) if ratio > 1e-9 else 1
    return int(max(1, min(4, N)))


# ----------------------------- метрики -----------------------------
def _maxwell_at(h, B, ix, R):
    """∇·B/|∇B| и |∇×B|/|∇B| в блоке (несколько срезов по x, ±R узлов по y,z) — np.gradient."""
    j0, k0 = h["ny"] // 2, h["nz"] // 2
    xa, xb = max(1, ix - 2), min(h["nx"] - 1, ix + 3)
    sub = B[xa:xb, j0 - R:j0 + R + 1, k0 - R:k0 + R + 1, :]
    if min(sub.shape[:3]) < 3:
        return float("nan"), float("nan")
    gx = [np.gradient(sub[..., c], h["sx"], axis=0) for c in range(3)]
    gy = [np.gradient(sub[..., c], h["sy"], axis=1) for c in range(3)]
    gz = [np.gradient(sub[..., c], h["sz"], axis=2) for c in range(3)]
    s = (slice(1, -1),) * 3                                       # внутренние точки
    div = (gx[0] + gy[1] + gz[2])[s]
    cx = (gy[2] - gz[1])[s]; cy = (gz[0] - gx[2])[s]; cz = (gx[1] - gy[0])[s]
    gnorm = np.sqrt(sum(g[s] ** 2 for g in (gx[0], gx[1], gx[2], gy[0], gy[1], gy[2], gz[0], gz[1], gz[2])))
    sc = np.sqrt(np.mean(gnorm ** 2)) + 1e-30
    return (np.sqrt(np.mean(div ** 2)) / sc, np.sqrt(np.mean(cx ** 2 + cy ** 2 + cz ** 2)) / sc)


def assess_field(path, order=None, beam_r=2e-3, verbose=True):
    """Полная оценка качества карты. Возвращает dict метрик; печатает отчёт при verbose."""
    h, B = load(path)
    ix = _body_slice(h, B)
    j0, k0 = h["ny"] // 2, h["nz"] // 2
    N = order if order else detect_order(h, B, ix)
    name = os.path.basename(path)

    # масштаб поля (внеосевое в теле)
    Rb = max(1, round(beam_r / h["sy"]))
    scale = np.sqrt(np.mean(B[ix, j0 - Rb:j0 + Rb + 1, k0 - Rb:k0 + Rb + 1, 1:3] ** 2)) + 1e-30

    # 1-2. Максвелл по радиусам
    radii = [2, 5, 10, 20]
    mw = {r: _maxwell_at(h, B, ix, max(1, round(r * 1e-3 / h["sy"]))) for r in radii}
    div2, curl2 = mw[2]

    # 3. поле на оси (для мультиполей N≥2 должно быть 0)
    onaxis = float(np.linalg.norm(B[ix, j0, k0, 1:3]))
    edge = float(np.linalg.norm(B[ix, j0 + Rb, k0, 1:3])) + 1e-30
    onaxis_rel = onaxis / edge if N >= 2 else float("nan")

    # 4. симметрия: Bx на плоскости x=mnx
    near = (np.add.outer((np.arange(h["ny"]) - j0) ** 2 * h["sy"] ** 2,
                         (np.arange(h["nz"]) - k0) ** 2 * h["sz"] ** 2) <= beam_r ** 2)
    bx0 = float(np.sqrt(np.mean(B[0, :, :, 0][near] ** 2)))
    jump_rel = 2 * bx0 / scale

    # 5. чистота спектра в зоне пучка
    by, _ = _circle(h, B, ix, max(beam_r, 4e-3))
    c = np.abs(np.fft.rfft(by)); c = c / (np.sum(c[1:7]) + 1e-30)
    kdes = max(0, N - 1)
    purity = float(c[kdes]) if kdes < len(c) else 0.0
    nonideal = float(np.sqrt(np.sum([c[k] ** 2 for k in range(1, 7) if k != kdes])))

    m = dict(name=name, order=N, field_scale_uT=scale * 1e6, grid=(h["nx"], h["ny"], h["nz"]),
             div_beam=div2, curl_beam=curl2, div_by_radius={r: mw[r][0] for r in radii},
             onaxis_uT=onaxis * 1e6, onaxis_rel=onaxis_rel,
             symmetry_jump_uT=2 * bx0 * 1e6, symmetry_jump_rel=jump_rel,
             purity=purity, nonideal=nonideal,
             noise_floor_uT=(onaxis * 1e6 if N >= 2 else float("nan")))

    if verbose:
        _report(m)
    return m


def _flag(val, good, warn):
    if val != val:  # nan
        return "  —  "
    return "  ✓  " if val < good else (" ⚠  " if val < warn else " ✗✗ ")


def _report(m):
    print(f"\n=== Качество поля: {m['name']}  (мультиполь N={m['order']}, "
          f"сетка {m['grid'][0]}×{m['grid'][1]}×{m['grid'][2]}, |B⊥|тело={m['field_scale_uT']:.2f} мкТл) ===")
    print(f"  1. ∇·B/|∇B| у пучка (r=2мм):    {m['div_beam']*100:6.1f}% {_flag(m['div_beam'],0.02,0.10)}"
          f"  (по радиусам 2/5/10/20мм: " + "/".join(f"{m['div_by_radius'][r]*100:.0f}%" for r in (2,5,10,20)) + ")")
    print(f"     ↳ падает с r => шум слабого поля; держится => невязка решения")
    print(f"  2. ∇×B/|∇B| у пучка (r=2мм):    {m['curl_beam']*100:6.1f}% {_flag(m['curl_beam'],0.05,0.20)}")
    if m['order'] >= 2:
        print(f"  3. поле на оси (должно=0):      {m['onaxis_uT']:6.3f} мкТл = {m['onaxis_rel']*100:5.1f}% "
              f"от поля на пучке {_flag(m['onaxis_rel'],0.05,0.20)}")
    else:
        print(f"  3. поле на оси:                 {m['onaxis_uT']:.3f} мкТл (диполь — поле на оси штатно ≠0)")
    print(f"  4. скачок симметрии 2·Bx(0):    {m['symmetry_jump_uT']:6.3f} мкТл = {m['symmetry_jump_rel']*100:5.2f}% "
          f"от поля {_flag(m['symmetry_jump_rel'],0.01,0.05)}")
    print(f"  5. чистота гл. гармоники:       {m['purity']*100:6.1f}%   неидеальность: {m['nonideal']*100:.1f}%")
    if m['order'] == 1:
        verdict = ("ДИПОЛЬ: ∇·B/∇×B/спектр ненадёжны (поле ~однородно, |∇B|≈0). Ключевое — симметрия: "
                   + ("в норме." if m['symmetry_jump_rel'] < 0.01 else "ПРОВЕРИТЬ."))
    elif m['onaxis_rel'] > 0.2:
        verdict = (f"пол ошибки ~{m['noise_floor_uT']:.2f} мкТл, СЛАБОЕ поле близко к полу — "
                   "считать при бо́льшем токе / мельче сетка / проверить сходимость солвера.")
    elif m['onaxis_rel'] < 0.05 and m['div_beam'] < 0.05:
        verdict = f"пол ошибки ~{m['noise_floor_uT']:.2f} мкТл; поле ≫ пола — карта пригодна."
    else:
        verdict = f"пол ошибки ~{m['noise_floor_uT']:.2f} мкТл; промежуточный случай — см. ∇·B по радиусам."
    print(f"  → {verdict}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        # демо: все карты проекта
        import glob
        for p in sorted(glob.glob("workspace/fields_bin/*.bin")):
            assess_field(p)
    else:
        order = int(sys.argv[2]) if len(sys.argv) > 2 else None
        assess_field(sys.argv[1], order=order)
