"""Core 2D GeoEnthalpy-Delta model 
Source: https://github.com/GeoJorge/GeoEnthalpy-Delta/blob/main/python/geoenthalpy2d.py

This module is shared by the 2D figure scripts for the manuscript. All
variables are nondimensional and beta-normalized as described in Sect. 4:

    x,y -> x/l, y/l
    H, eta, E, Z -> H/(beta*l), eta/(beta*l), ...
    t -> nu_top^x t/l^2
    q -> q/(beta*nu_top^x)

With this scaling, the single-slope non-erodible basement is E(x,y) = -x.
Thresholds C are physical slopes divided by beta, and the total 2D input
Q_in is scaled by beta*nu_top^x*l.

The diagnostic routines use a true four-neighbour connected-component filter.
The retained active deposit is the component associated with the feeder: the
component containing the most landward active cell within the feeder's lateral
footprint. If no component intersects that footprint, the component nearest to
the upstream feeder segment is selected, with sediment volume used as the tie
breaker. Detached components remain in H and in the sediment budget, but are
excluded from ABT, shoreline, toe, and planform-geometry diagnostics.
"""
from __future__ import annotations

__version__ = "1.0.0"

from collections import deque
from dataclasses import dataclass
from typing import List, Optional

import numpy as np

try:
    from numba import njit
    NUMBA_AVAILABLE = True
except Exception:  # pragma: no cover
    NUMBA_AVAILABLE = False

    def njit(*args, **kwargs):
        def wrapper(func):
            return func
        return wrapper


@dataclass(frozen=True)
class SeaLevel2D:
    kind: str = "fixed"     # "fixed" or "linear"
    value0: float = 0.0
    rate: float = 0.0       # Z(t)=value0+rate*t

    def __call__(self, t: float | np.ndarray) -> float | np.ndarray:
        tt = np.asarray(t, dtype=float)
        if self.kind == "fixed":
            out = self.value0 + np.zeros_like(tt)
        elif self.kind == "linear":
            out = self.value0 + self.rate * tt
        else:
            raise ValueError(f"Unknown sea-level kind: {self.kind}")
        return float(out) if np.isscalar(t) else out


@dataclass(frozen=True)
class Params2D:
    name: str = "case"
    Qin: float = 0.5
    feeder_width: float = 1.0
    Ctop_x: float = 0.0
    Ctop_y: float = 0.0
    Cfore_x: float = 2.0
    Cfore_y: float = 2.0
    Dtop_x: float = 1.0
    Dtop_y: float = 1.0
    Dfore_x: float = 1.0
    Dfore_y: float = 1.0
    sea_level: SeaLevel2D = SeaLevel2D()


@dataclass(frozen=True)
class Grid2D:
    x_min: float = -5.0
    x_max: float = 5.0
    y_min: float = -3.0
    y_max: float = 3.0
    dx: float = 0.08
    dy: float = 0.08
    t_end: Optional[float] = None
    target_volume: Optional[float] = None
    cfl: float = 0.40
    save_dt: float = 0.10
    h_tol: float = 1.0e-7


@dataclass
class Boundaries2D:
    abt: np.ndarray
    shoreline: np.ndarray
    toe: np.ndarray
    active_mask: np.ndarray
    component_count: int
    selected_component_cells: int
    selected_volume: float
    detached_volume: float
    detached_fraction: float
    row_gap_count: int
    selection_mode: str


@dataclass
class Result2D:
    p: Params2D
    g: Grid2D
    x: np.ndarray
    y: np.ndarray
    basement: np.ndarray
    times: np.ndarray
    H: np.ndarray                 # final field
    H_history: Optional[np.ndarray]
    boundaries: List[Boundaries2D]
    center_abt: np.ndarray
    center_shoreline: np.ndarray
    center_toe: np.ndarray
    volume: np.ndarray
    mass_error: np.ndarray
    metrics: dict
    diagnostics: dict
    dt: float


def make_xy(g: Grid2D) -> tuple[np.ndarray, np.ndarray]:
    x = np.arange(g.x_min, g.x_max + 0.5 * g.dx, g.dx)
    y = np.arange(g.y_min, g.y_max + 0.5 * g.dy, g.dy)
    return x, y


def planar_basement_2d(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Beta-normalized, fixed, non-erodible planar basement E=-x."""
    return -x[:, None] * np.ones((1, y.size))


def _interp_level(v0: float, v1: float, x0: float, x1: float, target: float) -> float:
    if abs(v1 - v0) < 1.0e-30:
        return float(0.5 * (x0 + x1))
    return float(x0 + (target - v0) * (x1 - x0) / (v1 - v0))


@njit(cache=False)
def _simulate_2d_core(basement, nt, dt, save_steps, feeder_mask, Qin,
                      sea_kind_code, sea_value0, sea_rate,
                      Ctop_x, Ctop_y, Cfore_x, Cfore_y,
                      Dtop_x, Dtop_y, Dfore_x, Dfore_y, dx, dy):
    """Time loop with positivity-limited lateral fluxes and sequential x export."""
    nx, ny = basement.shape
    H = np.zeros((nx, ny))
    Hnew = np.zeros((nx, ny))
    Hsave = np.zeros((save_steps.size, nx, ny))
    qx = np.zeros((nx + 1, ny))
    qy = np.zeros((nx, ny + 1))
    eta = np.zeros((nx, ny))
    fac_y = np.ones((nx, ny))

    nfeed = 0
    for j in range(ny):
        if feeder_mask[j]:
            nfeed += 1
    qdens = Qin / (nfeed * dy)
    ptr = 0

    for step in range(nt + 1):
        tnow = step * dt
        if sea_kind_code == 0:
            Z = sea_value0
        else:
            Z = sea_value0 + sea_rate * tnow

        for i in range(nx):
            for j in range(ny):
                eta[i, j] = basement[i, j] + H[i, j]

        while ptr < save_steps.size and step == save_steps[ptr]:
            for i in range(nx):
                for j in range(ny):
                    Hsave[ptr, i, j] = H[i, j]
            ptr += 1

        if step == nt:
            break

        qx[:, :] = 0.0
        qy[:, :] = 0.0
        fac_y[:, :] = 1.0

        for j in range(ny):
            if feeder_mask[j]:
                qx[0, j] = qdens

        # Candidate signed lateral fluxes, positive in +y.
        for i in range(nx):
            for j in range(ny - 1):
                slope = (eta[i, j] - eta[i, j + 1]) / dy
                donor_top = eta[i, j] >= Z if slope >= 0.0 else eta[i, j + 1] >= Z
                D = Dtop_y if donor_top else Dfore_y
                C = Ctop_y if donor_top else Cfore_y
                excess = abs(slope) - C
                if excess > 0.0:
                    qy[i, j + 1] = D * excess if slope >= 0.0 else -D * excess

        # Lateral donor limiter.
        tiny = 1.0e-30
        for i in range(nx):
            for j in range(ny):
                out_rate = 0.0
                if qy[i, j] < 0.0:
                    out_rate += -qy[i, j] / dy
                if qy[i, j + 1] > 0.0:
                    out_rate += qy[i, j + 1] / dy
                if out_rate > 0.0:
                    f = (H[i, j] / dt) / (out_rate + tiny)
                    if f < 1.0:
                        fac_y[i, j] = max(f, 0.0)
        for i in range(nx):
            for j in range(ny - 1):
                face = qy[i, j + 1]
                qy[i, j + 1] = fac_y[i, j] * face if face >= 0.0 else fac_y[i, j + 1] * face

        # Sequential downstream limiter, including incoming and lateral flux.
        for i in range(nx - 1):
            for j in range(ny):
                slope = (eta[i, j] - eta[i + 1, j]) / dx
                donor_top = eta[i, j] >= Z
                D = Dtop_x if donor_top else Dfore_x
                C = Ctop_x if donor_top else Cfore_x
                candidate = D * (slope - C) if slope > C else 0.0
                lateral_net = (qy[i, j] - qy[i, j + 1]) / dy
                available = qx[i, j] + dx * (H[i, j] / dt + lateral_net)
                qx[i + 1, j] = min(candidate, max(available, 0.0))

        # Conservative update. Basement is fixed and H is clipped only at zero.
        for i in range(nx):
            for j in range(ny):
                rate = ((qx[i, j] - qx[i + 1, j]) / dx
                        + (qy[i, j] - qy[i, j + 1]) / dy)
                Hnew[i, j] = max(H[i, j] + dt * rate, 0.0)
        H, Hnew = Hnew, H

    return Hsave


def _sea_kind_code(kind: str) -> int:
    if kind == "fixed":
        return 0
    if kind == "linear":
        return 1
    raise ValueError(f"Unknown sea-level kind: {kind}")


def _label_components_4(mask: np.ndarray, H: np.ndarray, feeder_mask: np.ndarray) -> tuple[np.ndarray, list[dict]]:
    """Label four-neighbour components and return deterministic component stats."""
    nx, ny = mask.shape
    labels = np.zeros((nx, ny), dtype=np.int32)
    feeder_js = np.flatnonzero(feeder_mask)
    stats: list[dict] = []
    label = 0

    for i0 in range(nx):
        for j0 in range(ny):
            if not mask[i0, j0] or labels[i0, j0] != 0:
                continue
            label += 1
            q: deque[tuple[int, int]] = deque([(i0, j0)])
            labels[i0, j0] = label
            cells: list[tuple[int, int]] = []
            volume_weight = 0.0
            min_i = nx
            min_feeder_i = nx + 1
            source_distance = np.inf

            while q:
                i, j = q.popleft()
                cells.append((i, j))
                volume_weight += float(H[i, j])
                min_i = min(min_i, i)
                if feeder_mask[j]:
                    min_feeder_i = min(min_feeder_i, i)
                if feeder_js.size:
                    lateral_distance = int(np.min(np.abs(feeder_js - j)))
                    source_distance = min(source_distance, float(i + lateral_distance))

                if i > 0 and mask[i - 1, j] and labels[i - 1, j] == 0:
                    labels[i - 1, j] = label; q.append((i - 1, j))
                if i + 1 < nx and mask[i + 1, j] and labels[i + 1, j] == 0:
                    labels[i + 1, j] = label; q.append((i + 1, j))
                if j > 0 and mask[i, j - 1] and labels[i, j - 1] == 0:
                    labels[i, j - 1] = label; q.append((i, j - 1))
                if j + 1 < ny and mask[i, j + 1] and labels[i, j + 1] == 0:
                    labels[i, j + 1] = label; q.append((i, j + 1))

            stats.append({
                "label": label,
                "cells": len(cells),
                "volume_weight": volume_weight,
                "min_i": min_i,
                "min_feeder_i": min_feeder_i,
                "intersects_feeder_footprint": min_feeder_i <= nx,
                "source_distance": source_distance,
            })
    return labels, stats


def feeder_associated_component(
    H: np.ndarray,
    h_tol: float,
    feeder_mask: np.ndarray,
    dx: float,
    dy: float,
) -> tuple[np.ndarray, dict]:
    """Return the feeder-associated four-neighbour active component.

    The finite-thickness field itself is not altered. The returned mask is used
    only for active-boundary and planform-geometry diagnostics.
    """
    active = H > h_tol
    total_active_volume = float(np.sum(H[active]) * dx * dy)
    if not np.any(active):
        return np.zeros_like(active, dtype=bool), {
            "component_count": 0,
            "selected_component_cells": 0,
            "selected_volume": 0.0,
            "detached_volume": 0.0,
            "detached_fraction": 0.0,
            "selection_mode": "empty",
        }

    labels, stats = _label_components_4(active, H, feeder_mask)
    feeder_candidates = [s for s in stats if s["intersects_feeder_footprint"]]
    if feeder_candidates:
        best = min(feeder_candidates, key=lambda s: (s["min_feeder_i"], -s["volume_weight"], s["label"]))
        mode = "feeder-footprint"
    else:
        best = min(stats, key=lambda s: (s["source_distance"], -s["volume_weight"], s["label"]))
        mode = "nearest-feeder"

    main_mask = labels == int(best["label"])
    selected_volume = float(np.sum(H[main_mask]) * dx * dy)
    detached_volume = max(total_active_volume - selected_volume, 0.0)
    detached_fraction = detached_volume / total_active_volume if total_active_volume > 0.0 else 0.0
    return main_mask, {
        "component_count": len(stats),
        "selected_component_cells": int(np.count_nonzero(main_mask)),
        "selected_volume": selected_volume,
        "detached_volume": detached_volume,
        "detached_fraction": detached_fraction,
        "selection_mode": mode,
    }


def diagnose_boundaries_2d(
    H: np.ndarray,
    eta: np.ndarray,
    Z: float,
    x: np.ndarray,
    y: np.ndarray,
    h_tol: float,
    feeder_mask: np.ndarray,
    dx: float,
    dy: float,
) -> Boundaries2D:
    """Diagnose ABT, SH, and TOE from the feeder-associated active deposit."""
    main_mask, cstats = feeder_associated_component(H, h_tol, feeder_mask, dx, dy)
    ny = H.shape[1]
    abt = np.full(ny, np.nan)
    shoreline = np.full(ny, np.nan)
    toe = np.full(ny, np.nan)
    row_gap_count = 0

    for j in range(ny):
        active = np.flatnonzero(main_mask[:, j])
        if active.size == 0:
            continue
        if active.size > 1 and np.any(np.diff(active) > 1):
            row_gap_count += 1
        i0 = int(active[0])
        i1 = int(active[-1])
        abt[j] = x[0] if i0 == 0 else _interp_level(H[i0 - 1, j], H[i0, j], x[i0 - 1], x[i0], h_tol)
        toe[j] = x[-1] if i1 == x.size - 1 else _interp_level(H[i1, j], H[i1 + 1, j], x[i1], x[i1 + 1], h_tol)

        phi = eta[:, j] - Z
        # Prefer a genuine adjacent crossing inside the selected component.
        crossings = [i for i in range(i0, i1)
                     if main_mask[i, j] and main_mask[i + 1, j]
                     and phi[i] >= 0.0 and phi[i + 1] < 0.0]
        if crossings:
            ish = int(crossings[-1])
            shoreline[j] = _interp_level(phi[ish], phi[ish + 1], x[ish], x[ish + 1], 0.0)
        else:
            wet = active[phi[active] >= 0.0]
            if wet.size:
                ish = int(wet[-1])
                shoreline[j] = x[-1] if ish == x.size - 1 else _interp_level(
                    phi[ish], phi[ish + 1], x[ish], x[ish + 1], 0.0
                )

    return Boundaries2D(
        abt=abt,
        shoreline=shoreline,
        toe=toe,
        active_mask=main_mask,
        component_count=int(cstats["component_count"]),
        selected_component_cells=int(cstats["selected_component_cells"]),
        selected_volume=float(cstats["selected_volume"]),
        detached_volume=float(cstats["detached_volume"]),
        detached_fraction=float(cstats["detached_fraction"]),
        row_gap_count=int(row_gap_count),
        selection_mode=str(cstats["selection_mode"]),
    )


def _range_or_nan(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    return float(np.ptp(values)) if values.size >= 2 else np.nan


def geometry_metrics(
    p: Params2D,
    g: Grid2D,
    x: np.ndarray,
    y: np.ndarray,
    H: np.ndarray,
    eta: np.ndarray,
    b: Boundaries2D,
    time: float,
) -> dict:
    """Compute all final-state diagnostics used by Figs. 6, 8, and 9."""
    j0 = int(np.argmin(np.abs(y)))
    volume = float(np.sum(H) * g.dx * g.dy)
    expected = p.Qin * time
    mass_error = abs(volume - expected) / max(expected, np.finfo(float).eps) if time > 0 else 0.0
    active_y = np.any(b.active_mask, axis=0)
    fan_width = _range_or_nan(y[active_y])
    sh_width = _range_or_nan(y[np.isfinite(b.shoreline)])
    toe_width = _range_or_nan(y[np.isfinite(b.toe)])
    Z = float(p.sea_level(time))

    center_abt = float(b.abt[j0]) if np.isfinite(b.abt[j0]) else np.nan
    center_sh = float(b.shoreline[j0]) if np.isfinite(b.shoreline[j0]) else np.nan
    center_toe = float(b.toe[j0]) if np.isfinite(b.toe[j0]) else np.nan
    center_sep = center_toe - center_sh if np.isfinite(center_toe) and np.isfinite(center_sh) else np.nan

    valid_abt = np.isfinite(b.abt)
    landward_extent = float(np.nanmin(b.abt[valid_abt])) if np.any(valid_abt) else np.nan
    fan_length = center_toe - landward_extent if np.isfinite(center_toe) and np.isfinite(landward_extent) else np.nan
    source_to_toe = center_toe - x[0] if np.isfinite(center_toe) else np.nan
    width_length_ratio = fan_width / fan_length if np.isfinite(fan_width) and fan_length > 0 else np.nan
    opening_angle_deg = float(np.degrees(2.0 * np.arctan2(0.5 * fan_width, source_to_toe))) \
        if np.isfinite(fan_width) and source_to_toe > 0 else np.nan
    deposit_opening_angle_deg = float(np.degrees(2.0 * np.arctan2(0.5 * fan_width, fan_length))) \
        if np.isfinite(fan_width) and fan_length > 0 else np.nan

    # Sediment budget includes the complete H field, including detached remnants.
    topset_mask = (H > g.h_tol) & (eta >= Z)
    topset_volume = float(np.sum(H[topset_mask]) * g.dx * g.dy)
    foreset_volume = max(volume - topset_volume, 0.0)

    return {
        "time": float(time),
        "volume": volume,
        "expected_volume": expected,
        "mass_error": mass_error,
        "final_sea_level": Z,
        "center_abt": center_abt,
        "center_shoreline": center_sh,
        "center_toe": center_toe,
        "center_shore_to_toe": center_sep,
        "fan_width": fan_width,
        "shoreline_width": sh_width,
        "toe_width": toe_width,
        "landward_extent": landward_extent,
        "fan_length": fan_length,
        "width_length_ratio": width_length_ratio,
        "opening_angle_deg": opening_angle_deg,
        "deposit_opening_angle_deg": deposit_opening_angle_deg,
        "max_thickness": float(np.max(H)),
        "topset_volume": topset_volume,
        "foreset_volume": foreset_volume,
        "foreset_fraction": float(foreset_volume / volume) if volume > 0 else np.nan,
        "component_count": b.component_count,
        "selected_component_cells": b.selected_component_cells,
        "selected_volume": b.selected_volume,
        "detached_volume": b.detached_volume,
        "detached_fraction": b.detached_fraction,
        "row_gap_count": b.row_gap_count,
        "selection_mode": b.selection_mode,
    }


def smooth_gradient(values: np.ndarray, times: np.ndarray, half_window: int = 15) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    good = np.isfinite(values)
    if good.sum() < 3:
        return np.full_like(values, np.nan)
    filled = np.interp(times, times[good], values[good])
    out = np.full_like(values, np.nan)
    for i in range(values.size):
        lo = max(0, i - half_window)
        hi = min(values.size, i + half_window + 1)
        if hi - lo >= 3:
            out[i] = np.polyfit(times[lo:hi], filled[lo:hi], 1)[0]
    return out


def curve_rms(curve: np.ndarray, ref: np.ndarray) -> float:
    valid = np.isfinite(curve) & np.isfinite(ref)
    if valid.sum() < 3:
        return np.nan
    return float(np.sqrt(np.mean((curve[valid] - ref[valid]) ** 2)))


def _first_persistent_plateau(values: np.ndarray, start: int, tolerance: float, tail_reference: float) -> int:
    """First index after start from which values remain within tolerance."""
    n = values.size
    for i in range(max(start, 0), n):
        tail = values[i:]
        finite = np.isfinite(tail)
        if np.any(finite) and np.nanmax(np.abs(tail[finite] - tail_reference)) <= tolerance:
            return i
    return n - 1


def constant_slr_diagnostics(result: Result2D) -> dict:
    """Centralized phase, translation, and collapse diagnostics for Figure 7."""
    t = result.times
    r = result.center_abt
    s = result.center_shoreline
    u = result.center_toe
    a = result.p.sea_level.rate
    if t.size < 3 or not np.any(np.isfinite(s)) or not np.any(np.isfinite(u)):
        return {}

    i_sr = int(np.nanargmax(s))
    # The practical toe-abandonment time is the first time at which the
    # centerline TOE reaches its final maximum position. This reproduces the
    # phase definition used in the original Figure 7 workflow while keeping
    # the diagnostic centralized in the core.
    i_toe = int(np.nanargmax(u))
    tail_n = max(5, t.size // 10)
    toe_ref = float(np.nanmedian(u[-tail_n:]))
    tol = max(2.0 * result.g.dx, 0.03)

    rxi = r + a * t
    sxi = s + a * t
    refs = (
        float(np.nanmedian(rxi[-tail_n:])),
        float(np.nanmedian(sxi[-tail_n:])),
        toe_ref,
    )
    i_post = t.size - 1
    for i in range(max(i_toe, i_sr + 1), t.size - 2):
        good = True
        for k in range(i, min(i + 3, t.size)):
            vals = (rxi[k], sxi[k], u[k])
            if not all(np.isfinite(vals)) or max(abs(vals[0]-refs[0]), abs(vals[1]-refs[1]), abs(vals[2]-refs[2])) > tol:
                good = False
                break
        if good:
            i_post = i
            break

    i_auto = int(round(0.5 * (i_sr + i_toe))) if i_toe > i_sr else min(i_sr + 1, t.size - 1)
    i_auto = max(i_sr, min(i_auto, i_post))
    stage_indices = [i_sr, i_auto, i_post]

    late_start = max(i_post, int(0.75 * (t.size - 1)))
    late = np.arange(late_start, t.size)
    def slope(v: np.ndarray) -> float:
        good = np.isfinite(v[late])
        return float(np.polyfit(t[late][good], v[late][good], 1)[0]) if good.sum() >= 3 else np.nan

    speed_abt = slope(r)
    speed_sh = slope(s)
    speed_toe = slope(u)

    late_indices = np.unique(np.linspace(i_post, t.size - 1, 4).round().astype(int))
    b_ref = result.boundaries[-1]
    collapse_values: list[float] = []
    for idx in late_indices:
        b = result.boundaries[idx]
        collapse_values.extend([
            curve_rms(b.abt + a*t[idx], b_ref.abt + a*t[-1]),
            curve_rms(b.shoreline + a*t[idx], b_ref.shoreline + a*t[-1]),
        ])
    collapse_values = [v for v in collapse_values if np.isfinite(v)]
    collapse_rms = float(np.sqrt(np.mean(np.square(collapse_values)))) if collapse_values else np.nan
    normalized_collapse = collapse_rms / result.metrics["fan_width"] \
        if np.isfinite(collapse_rms) and result.metrics["fan_width"] > 0 else np.nan

    return {
        "shoreline_reversal_index": i_sr,
        "shoreline_reversal_time": float(t[i_sr]),
        "toe_abandonment_index": i_toe,
        "toe_abandonment_time": float(t[i_toe]),
        "post_autobreak_index": i_post,
        "post_autobreak_time": float(t[i_post]),
        "stage_indices": stage_indices,
        "stage_times": [float(t[i]) for i in stage_indices],
        "stage_names": ["progradation", "autoretreat", "post-autobreak"],
        "expected_translation_speed": float(-a),
        "late_abt_speed": speed_abt,
        "late_shoreline_speed": speed_sh,
        "late_toe_speed": speed_toe,
        "late_abt_speed_error": float(abs(speed_abt + a)) if np.isfinite(speed_abt) else np.nan,
        "late_shoreline_speed_error": float(abs(speed_sh + a)) if np.isfinite(speed_sh) else np.nan,
        "late_toe_speed_error": float(abs(speed_toe)) if np.isfinite(speed_toe) else np.nan,
        "late_planform_collapse_rms": collapse_rms,
        "late_planform_collapse_normalized": normalized_collapse,
        "late_indices": late_indices.tolist(),
        "late_times": [float(t[i]) for i in late_indices],
    }


def choose_constant_slr_stages(result: Result2D) -> list[int]:
    """Return the centralized representative stage indices used by Figure 7."""
    diag = result.diagnostics.get("constant_slr", {})
    if diag:
        return list(diag["stage_indices"])
    return [0, max(0, len(result.times)//2), len(result.times)-1]


def metrics_record(result: Result2D) -> dict:
    """Flatten parameters and centralized metrics for reproducible CSV output."""
    p = result.p
    row = {
        "case": p.name,
        "Qin": p.Qin,
        "feeder_width": p.feeder_width,
        "Ctop_x": p.Ctop_x,
        "Ctop_y": p.Ctop_y,
        "Cfore_x": p.Cfore_x,
        "Cfore_y": p.Cfore_y,
        "Dtop_x": p.Dtop_x,
        "Dtop_y": p.Dtop_y,
        "Dfore_x": p.Dfore_x,
        "Dfore_y": p.Dfore_y,
        "sea_level_kind": p.sea_level.kind,
        "sea_level_rate": p.sea_level.rate,
    }
    row.update(result.metrics)
    return row


def simulate_2d(p: Params2D, g: Grid2D, save_history: bool = True) -> Result2D:
    x, y = make_xy(g)
    basement = planar_basement_2d(x, y)
    t_end = g.t_end if g.t_end is not None else (g.target_volume / p.Qin)
    if t_end is None:
        raise ValueError("Either g.t_end or g.target_volume must be specified.")
    Dmax = max(p.Dtop_x, p.Dtop_y, p.Dfore_x, p.Dfore_y, 1.0e-12)
    dt_stable = g.cfl / (2.0 * Dmax * (1.0 / g.dx**2 + 1.0 / g.dy**2))
    nt = int(np.ceil(t_end / dt_stable))
    dt = t_end / nt
    if save_history:
        req = np.arange(0.0, t_end + 0.5 * g.save_dt, g.save_dt)
        if req.size == 0 or req[-1] < t_end - 0.5*dt:
            req = np.r_[req, t_end]
        save_steps = np.unique(np.clip(np.rint(req / dt).astype(int), 0, nt))
    else:
        save_steps = np.array([nt], dtype=int)
    times = save_steps * dt
    feeder = np.abs(y) <= 0.5 * p.feeder_width + 1.0e-12
    if not np.any(feeder):
        raise ValueError(f"{p.name}: feeder_width is smaller than grid spacing.")

    Hsave = _simulate_2d_core(
        basement, nt, dt, save_steps.astype(np.int64), feeder, p.Qin,
        _sea_kind_code(p.sea_level.kind), p.sea_level.value0, p.sea_level.rate,
        p.Ctop_x, p.Ctop_y, p.Cfore_x, p.Cfore_y,
        p.Dtop_x, p.Dtop_y, p.Dfore_x, p.Dfore_y, g.dx, g.dy,
    )

    boundaries: List[Boundaries2D] = []
    center_abt: list[float] = []
    center_sh: list[float] = []
    center_toe: list[float] = []
    volumes: list[float] = []
    mass_errors: list[float] = []
    j0 = int(np.argmin(np.abs(y)))
    for Hk, tt in zip(Hsave, times):
        eta = basement + Hk
        b = diagnose_boundaries_2d(
            Hk, eta, float(p.sea_level(tt)), x, y, g.h_tol, feeder, g.dx, g.dy
        )
        boundaries.append(b)
        center_abt.append(b.abt[j0])
        center_sh.append(b.shoreline[j0])
        center_toe.append(b.toe[j0])
        vol = float(np.sum(Hk) * g.dx * g.dy)
        volumes.append(vol)
        mass_errors.append(abs(vol - p.Qin*tt) / max(p.Qin*tt, np.finfo(float).eps) if tt > 0 else 0.0)

    H_final = Hsave[-1]
    eta_final = basement + H_final
    result = Result2D(
        p=p, g=g, x=x, y=y, basement=basement, times=times,
        H=H_final, H_history=Hsave if save_history else None,
        boundaries=boundaries,
        center_abt=np.asarray(center_abt),
        center_shoreline=np.asarray(center_sh),
        center_toe=np.asarray(center_toe),
        volume=np.asarray(volumes),
        mass_error=np.asarray(mass_errors),
        metrics={}, diagnostics={}, dt=float(dt),
    )
    result.metrics = geometry_metrics(p, g, x, y, H_final, eta_final, boundaries[-1], float(times[-1]))
    result.diagnostics = {
        "connectivity": {
            "max_component_count": int(max(b.component_count for b in boundaries)),
            "max_detached_fraction": float(max(b.detached_fraction for b in boundaries)),
            "max_row_gap_count": int(max(b.row_gap_count for b in boundaries)),
            "final_selection_mode": boundaries[-1].selection_mode,
        }
    }
    if save_history and p.sea_level.kind == "linear" and p.sea_level.rate > 0.0:
        result.diagnostics["constant_slr"] = constant_slr_diagnostics(result)
    return result
