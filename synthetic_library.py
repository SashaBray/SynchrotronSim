# -*- coding: utf-8 -*-
"""Менеджер библиотеки элементов: авто-добавление СИНТЕТИЧЕСКИХ элементов под эффективные длины из LTE.

Технологи задают в LTE эффективные длины магнитов, которых нет в библиотеке (она — из старой КД).
Этот модуль при разборе LTE находит для каждого элемента библиотечный эталон нужного типа; если его
L_eff не совпадает с требуемой (в пределах допуска), создаёт НОВЫЙ элемент, масштабируя эталон:
  L_eff = L из LTE;  L_phys = phys_эталона + (L_eff − eff_эталона)   — «сплайс тела» (фриндж сохранён).
Многообмоточные (DQ, SP) масштабируются на общую Δ (обмотки сохраняют относительные смещения длин).
Новый элемент помечается synthetic=True, scaled_from=<эталон>; kind и матрица связи копируются с эталона.
Опционально (gen_assets) — синтетическая карта .fld (гибридный ресайз тела) и STL (продольн. масштаб).

Идемпотентно: повторный разбор LTE переиспользует ранее созданные синтетические элементы из БД."""
import os
import csv

TB = os.path.join("workspace", "tables")


def _read(name):
    rows = list(csv.reader(open(os.path.join(TB, name), encoding="utf-8"), delimiter=";"))
    return rows[0], rows[1:]


def _write(name, header, rows):
    with open(os.path.join(TB, name), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter=";"); w.writerow(header); w.writerows(rows)


class LibraryManager:
    """Загружает hardware_library / element_transforms / coupling_matrix, добавляет синтетические
    элементы под требуемые L_eff, сохраняет обратно (только если что-то добавлено)."""

    def __init__(self, tol_m=1e-3, gen_assets=False, verbose=True):
        self.tol = tol_m; self.gen_assets = gen_assets; self.verbose = verbose
        self.hw_head, self.hw_rows = _read("hardware_library.csv")
        self.tr_head, self.tr_rows = _read("element_transforms.csv")
        self.cp_head, self.cp_rows = _read("coupling_matrix.csv")
        self.i = {c: k for k, c in enumerate(self.hw_head)}
        for req in ("synthetic", "scaled_from"):
            if req not in self.i:
                raise RuntimeError(f"hardware_library.csv нет столбца {req} — добавьте synthetic/scaled_from")
        self.hw = {r[0]: r for r in self.hw_rows}
        self.tr = {r[0]: r for r in self.tr_rows}
        self.cp = {r[0]: r for r in self.cp_rows}
        self.created = []; self.dirty = False

    def eff(self, did):
        return float(self.hw[did][self.i["eff_length"]])

    def phys(self, did):
        return float(self.hw[did][self.i["phys_length"]])

    def resolve(self, ref_ids, target_eff):
        """ref_ids — обмотки-эталоны (из map_to_hw); target_eff — L из LTE (м).
        Возврат: список фактических device_type_id (эталон, если совпал по L_eff; иначе синтетика)."""
        if not ref_ids:
            return []
        primary = ref_ids[0]
        if primary not in self.hw:
            return list(ref_ids)                          # эталон не в библиотеке — не трогаем
        if abs(self.eff(primary) - target_eff) <= self.tol:
            return list(ref_ids)                          # существующий подходит
        delta = target_eff - self.eff(primary)            # общий сдвиг длины (сплайс тела)
        suffix = f"L{round(target_eff * 1000)}"           # общий для всех обмоток магнита
        out = []
        for ref in ref_ids:
            nid = f"{ref}_{suffix}"
            if nid not in self.hw:
                self._create(ref, nid, self.eff(ref) + delta, target_eff)
            out.append(nid)
        return out

    def _create(self, ref, nid, new_eff, target_eff):
        i = self.i; row = list(self.hw[ref]); row[0] = nid
        row[i["eff_length"]] = f"{new_eff:.4f}"
        row[i["phys_length"]] = f"{self.phys(ref) + (new_eff - self.eff(ref)):.4f}"
        row[i["map_bin_path"]] = "-"; row[i["map_raw_path"]] = "-"   # без карты: матр. модель работает по L_eff
        row[i["synthetic"]] = "True"; row[i["scaled_from"]] = ref
        if self.gen_assets:                              # флаг --gen-assets: синтетические карта .bin и STL
            import synthetic_assets as sa
            binp = sa.gen_field_map(ref, nid, new_eff, self.hw[ref], i)
            if binp:
                row[i["map_bin_path"]] = binp
            stlp = sa.gen_stl(ref, nid, new_eff, self.hw[ref], i)
            if stlp:
                row[i["stl_path"]] = stlp                # синтетическая 3D-модель (иначе — STL эталона)
        self.hw[nid] = row; self.hw_rows.append(row)
        if ref in self.tr:                                # element_transforms: копия kind/model
            tr = list(self.tr[ref]); tr[0] = nid
            if len(tr) > 5:
                tr[5] = f"синтетика L_eff={new_eff*1000:.0f}мм (масштаб от {ref})"
            self.tr[nid] = tr; self.tr_rows.append(tr)
        if ref in self.cp:                                # coupling_matrix: per-амперные коэф. длино-независимы
            cp = list(self.cp[ref]); cp[0] = nid; self.cp[nid] = cp; self.cp_rows.append(cp)
        self.created.append((nid, ref, new_eff)); self.dirty = True
        if self.verbose:
            asset = " +карта+STL" if self.gen_assets else ""
            print(f"    + синтетика {nid}: L_eff={new_eff*1000:.1f}мм (от {ref}){asset}")

    def save(self):
        if not self.dirty:
            return 0
        _write("hardware_library.csv", self.hw_head, self.hw_rows)
        _write("element_transforms.csv", self.tr_head, self.tr_rows)
        _write("coupling_matrix.csv", self.cp_head, self.cp_rows)
        return len(self.created)
