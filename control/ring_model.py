# -*- coding: utf-8 -*-
"""Единый источник ring-производных матриц отклика — одно место, знающее, как зондировать кольцо.

Оборачивает ring_cpp.CppRing (через facility_report) + карты индексов (корректоры/мониторы/квадруполи)
и отдаёт applier/reader для response.build_response. Через него проходят ВСЕ матрицы отклика, так что
россыпь построителей (orbit_correct/regulator/xbpm_svd/shift_basis/loco) сводится к одному механизму.

Актуаторы: corr.NNN.hkick / corr.NNN.vkick (кики корректоров), quad.NNN.dx / quad.NNN.dy (сдвиг
квадруполя через feed-down: кик = K1L·сдвиг). Наблюдаемые: ebpm.x / ebpm.y (орбита в мониторах),
plane.<row>.<x|px|y|py> (любая плоскость — основа XBPM-селектора). Добавить тип = одна ветка.
"""
import os

os.environ.setdefault("SS_SECTOR", "1")

import numpy as np

import errors as er
import facility_report as fr
import orbit_correct as oc
import orbit_infer as oi

from . import response as resp


class RingModel:
    def __init__(self, lattice=26, nsl=8, n_photon=0, photon_d=15.0):
        self.lattice = lattice
        self.nsl = nsl
        inf = dict(fr.info_of(lattice)); inf["_id"] = lattice
        self.ring = fr._ring(inf, nsl=nsl)
        self.corr = oc.corrector_indices(self.ring)
        self.mon = oi.monitor_rows(self.ring)
        self.quads = oi.quad_indices(self.ring, kinds=("quad",))
        self.K1L = np.array([float(self.ring.mags[q].k1) * float(self.ring.mags[q].Leff) for q in self.quads])
        self.hk_ids = [f"corr.{i:03d}.hkick" for i in range(len(self.corr))]
        self.vk_ids = [f"corr.{i:03d}.vkick" for i in range(len(self.corr))]
        self.qx_ids = [f"quad.{i:03d}.dx" for i in range(len(self.quads))]
        self.qy_ids = [f"quad.{i:03d}.dy" for i in range(len(self.quads))]
        self._act = {}
        for i, c in enumerate(self.corr):
            self._act[self.hk_ids[i]] = ("hkick", c, 1.0)
            self._act[self.vk_ids[i]] = ("vkick", c, 1.0)
        for i, q in enumerate(self.quads):
            self._act[self.qx_ids[i]] = ("hkick", q, float(self.K1L[i]))   # сдвиг dx → гор. feed-down кик K1L·dx
            self._act[self.qy_ids[i]] = ("vkick", q, -float(self.K1L[i]))  # сдвиг dy → верт. feed-down (знак как в regulator)
        # источники для XBPM: центроид фотона = x + d·x' в центре ондулятора (плечо d до экрана)
        self.photon_rows, self.photon_d = self._source_rows(n_photon, photon_d) if n_photon else (np.array([], int), photon_d)

    def applier(self, ring, aid, eps):
        kind, idx, scale = self._act[aid]
        if kind == "hkick":
            ring.apply_perturb(kicks=[(idx, scale * eps)])
        else:
            ring.apply_perturb(vkicks=[(idx, scale * eps)])

    def _source_rows(self, n_photon, d):
        """Строки track_elements в центрах первых n_photon ондуляторов (ближайший mag.center к pivot группы)."""
        cen = np.array([m.center for m in self.ring.mags])
        groups = sorted(er.load_groups(self.lattice), key=lambda x: x["group_id"])[:n_photon]
        rows = [int(np.argmin(np.linalg.norm(cen - np.asarray(g["pivot"], float), axis=1))) + 1 for g in groups]
        return np.array(rows, int), float(d)

    def reader(self, oid, orbit):
        if oid == "ebpm.x":
            return orbit[self.mon, 0]
        if oid == "ebpm.y":
            return orbit[self.mon, 2]
        if oid == "photon.x":                                # центроид фотона: x + d·x' в источниках
            return orbit[self.photon_rows, 0] + self.photon_d * orbit[self.photon_rows, 1]
        if oid == "photon.y":
            return orbit[self.photon_rows, 2] + self.photon_d * orbit[self.photon_rows, 3]
        if oid.startswith("plane."):                         # plane.<row>.<x|px|y|py> — произвольная плоскость
            _, row, comp = oid.split(".")
            col = {"x": 0, "px": 1, "y": 2, "py": 3}[comp]
            return np.array([orbit[int(row), col]])
        raise ValueError(f"RingModel: наблюдаемая {oid} не поддержана")

    def response(self, actuators, observables, rebuild=False):
        key = (self.lattice, self.nsl, tuple(int(r) for r in self.photon_rows), self.photon_d)
        return resp.cached_build(key, self.ring, list(actuators), list(observables),
                                 self.applier, self.reader, rebuild=rebuild)

    def close(self):
        self.ring.apply_perturb()
        self.ring.close()
