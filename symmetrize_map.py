# -*- coding: utf-8 -*-
"""Симметризация карты поля секступоля по чётностям (денойз с сохранением симметрии).
Нормальный секступоль (ориентирован по сетке: Bz max на осях, By(y,0)≡0) имеет строгие чётности:
  Bz — ЧЁТНО по y и по z            (Bz=K(y²−z²))
  By — НЕЧЁТНО по y и по z          (By=2Kyz)
  Bx — ЧЁТНО по y, НЕЧЁТНО по z     (из ∇·B=0, ∇×B=0)
Усредняем 4 отражённых точки (±y,±z) с этими знаками ⇒ проекция на секступольно-симметричное
поле: убирает косые/чётнопорядковые (октуполь) шумовые компоненты, но оставляет остаточный шум
разрешённой чётности (поле «шумное, но симметричное»). Точные повороты 120° (полная D₃) в узлы
прямоугольной сетки не попадают — используем грид-точные отражения (4 партнёра, 2 на осях).
РАБОТАЕТ НА КОПИИ массива; оригинальные файлы не трогаются. Пишет .bin-копию."""
import os, struct, sys
import numpy as np
from scipy.ndimage import gaussian_filter
from field_quality import load
from _common import ansys_dir

HDR = 96


def smooth(B, sigma=1.0):
    """Сглаживание по соседним узлам (свёртка Гаусса, σ в узлах) с СОХРАНЕНИЕМ полной энергии.
    Линейная свёртка коммутирует с производными ⇒ ∇·B=0, ∇×B=0 сохраняются; симметричное ядро
    сохраняет симметрию поля. Глобальный ресейл возвращает суммарную энергию к исходной."""
    E0 = float(np.sum(B ** 2))
    S = np.empty_like(B)
    for c in range(3):
        S[..., c] = gaussian_filter(B[..., c], sigma=sigma, mode="nearest")
    E1 = float(np.sum(S ** 2))
    a = np.sqrt(E0 / E1) if E1 > 0 else 1.0
    S *= a
    print(f"  сглаживание σ={sigma} узлов: энергия сохранена (ресейл ×{a:.5f}); "
          f"убрано ВЧ-шума {(1-E1/E0)*100:.2f}% (до ресейла)")
    return S


# чётности (sy, sz) компонент (Bx, By, Bz) по ОТРАЖЕНИЯМ (грид-точны): +1 чёт, −1 нечёт, None — не проецируем.
# Поле мультиполя n: W=Bz+iBy=Cₙ(y+iz)^{n−1} (аналитично).
#   Квадруполь n=2: Bz=Gy, By=Gz.          Секступоль n=3: Bz=K(y²−z²), By=2Kyz.
#   Октуполь n=4:   Bz=y³−3yz², By=3y²z−z³.
# ВАЖНО: КВАДРУПОЛЬ и ОКТУПОЛЬ имеют ОДИНАКОВЫЕ чётности отражений (обе чёт.порядка) ⇒ отражения их НЕ
# разделяют. Разделяет грид-точный ПОВОРОТ 90° (на квадратной сетке): под ним квадруполь→−B (n=2),
# октуполь→+B (n=4). Секступоль (3-складчатая, поворот 120° не грид-точен) — только отражения.
# Диполь-квадруполь (комбинир.-функц.): чётность по y СМЕШАНА ⇒ по y НЕ проецируем, только средняя плоскость z.
PARITIES = {
    "quad": ((-1, -1), (+1, -1), (-1, +1)),         # как у окт по отражениям; разделяет поворот 90°
    "sext": ((+1, -1), (-1, -1), (+1, +1)),
    "oct":  ((-1, -1), (+1, -1), (-1, +1)),
    "dq":   ((None, -1), (None, -1), (None, +1)),   # только средняя плоскость (z)
}
ROT90 = {"quad": -1.0, "oct": +1.0}                 # знак поля под поворотом 90° (n чёт.); sext/dq — без поворота


def rot90(B):
    """Поворот поля на 90° CCW вокруг продольной оси X (квадратная сетка Y-Z): перестановка узлов +
    поворот вектора (By,Bz)→(−Bz,By), Bx неизм. R90(B)[p]=Rᵥ·B[R⁻¹p]."""
    R = np.rot90(B, k=1, axes=(1, 2))               # переставить узлы плоскости Y-Z (Y=ось1, Z=ось2)
    out = np.empty_like(R)
    out[..., 0] = R[..., 0]                          # Bx вдоль оси поворота
    out[..., 1] = -R[..., 2]                         # By' = −Bz
    out[..., 2] = R[..., 1]                          # Bz' =  By
    return out


def symmetrize(h, B, kind="sext"):
    """B[nx,ny,nz,3] → проекция на симметрию мультиполя (kind: quad|sext|oct|dq). Отражения по чётностям
    + грид-точный поворот 90° для quad/oct (разделяет n=2 и n=4). Печатает долю сохранённой энергии."""
    ymin, ymax = h["mny"], h["mny"] + (h["ny"] - 1) * h["sy"]
    zmin, zmax = h["mnz"], h["mnz"] + (h["nz"] - 1) * h["sz"]
    assert abs(zmin + zmax) < 1e-9, "сетка не симметрична по z"
    par = PARITIES[kind]
    fy = np.flip(B, axis=1); fz = np.flip(B, axis=2); fyz = np.flip(fy, axis=2)
    S = np.empty_like(B)
    for c in range(3):
        sy, sz = par[c]
        if sy is None:                                     # только z (комбинир.-функц.)
            S[..., c] = 0.5 * (B[..., c] + sz * fz[..., c])
        else:
            assert abs(ymin + ymax) < 1e-9, "сетка не симметрична по y"
            S[..., c] = 0.25 * (B[..., c] + sy * fy[..., c] + sz * fz[..., c] + sy * sz * fyz[..., c])
    if kind in ROT90:                                      # разделить n=2/n=4 поворотом 90° (√сетка)
        assert h["ny"] == h["nz"] and abs(h["sy"] - h["sz"]) < 1e-12, "поворот 90° требует квадратной сетки"
        S = 0.5 * (S + ROT90[kind] * rot90(S))
    for i, nm in enumerate(("Bx", "By", "Bz")):
        e_tot = float(np.sum(B[..., i] ** 2)) + 1e-30; e_kept = float(np.sum(S[..., i] ** 2))
        e_rem = float(np.sum((B[..., i] - S[..., i]) ** 2))
        print(f"  {nm}: сохранено {e_kept/e_tot*100:5.1f}% энергии, убрано {e_rem/e_tot*100:5.1f}%")
    return S


def write_bin(h, B, path, nominal_arg=1.0):
    """nominal_arg — ток снятия карты (offset 88): C++ currentScale=arg_val/nominalArg; при 0 → NaN.
    Для 100А-карт nominal_arg=100 (как у корректоров); PM-диполи — 1.0 (масштаб через correct_coeff)."""
    mxy = h["mny"] + (h["ny"] - 1) * h["sy"]; mxz = h["mnz"] + (h["nz"] - 1) * h["sz"]
    with open(path, "wb") as f:
        f.write(struct.pack("<9d", h["mnx"], h["mny"], h["mnz"], h["mxx"], mxy, mxz,
                            h["sx"], h["sy"], h["sz"]))              # 72 байта
        f.write(struct.pack("<3i", h["nx"], h["ny"], h["nz"]))      # → 84
        f.write(b"\x00" * (88 - 84))                                # паддинг до смещения 88
        f.write(struct.pack("<d", float(nominal_arg)))              # nominalArg на 88 (→ 96)
        f.write(np.ascontiguousarray(B, dtype="<f8").tobytes())


def main():
    src = sys.argv[1] if len(sys.argv) > 1 else str(ansys_dir() / "KD" / "SP" / "SP_I2_100A_v4.fld")
    kind = sys.argv[2] if len(sys.argv) > 2 else "sext"
    sigma = float(sys.argv[3]) if len(sys.argv) > 3 else 1.0
    print(f"вход:  {os.path.basename(src)}  (мультиполь: {kind})")
    h, B = load(src)
    print(f"сетка: {h['nx']}×{h['ny']}×{h['nz']}, шаг {h['sy']*1e3:.1f}мм")
    S = symmetrize(h, B.copy(), kind)                      # КОПИЯ, оригинал не трогаем
    write_bin(h, S, os.path.splitext(src)[0] + "_sym.bin")
    print(f"  → {os.path.splitext(os.path.basename(src))[0]}_sym.bin")
    SS = symmetrize(h, smooth(S, sigma), kind)             # сгладить + пере-симметризовать
    write_bin(h, SS, os.path.splitext(src)[0] + "_symsm.bin")
    print(f"  → {os.path.splitext(os.path.basename(src))[0]}_symsm.bin (симметризов. + сглажено)")


if __name__ == "__main__":
    main()
