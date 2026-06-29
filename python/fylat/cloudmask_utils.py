"""Cloud mask utility functions — exact Fortran algorithm replicas.

Ported from Fortran source files:
  conf_test.f, conf_test_2val.f, tview.f, get_sg_thresholds.f90
"""

import os
from typing import Dict, List, Optional, Tuple

import numpy as np
import yaml


# ===========================================================================
# conf_test — exact Fortran replica (conf_test.f)
# ===========================================================================

def _conf_test_scalar(val, locut, hicut, power, midpt):
    """Scalar conf_test matching Fortran logic exactly.

    coeff = 2**(power-1)
    S-curve split at midpt (beta):
      below beta:  s1 = (val-alpha)/(2*(beta-alpha)), c = coeff * s1^power
      above beta:  s1 = (val-gamma)/(2*(beta-gamma)), c = 1 - coeff * s1^power
    If flipped (locut > hicut): swap alpha/gamma and invert.
    """
    coeff = 2.0 ** (power - 1.0)

    if hicut > locut:
        alpha, gamma, flipped = locut, hicut, False
    else:
        alpha, gamma, hicut, flipped = hicut, locut, True

    beta = midpt

    if not flipped:
        if val > gamma:
            return 1.0
        elif val < alpha:
            return 0.0
    else:
        if val > gamma:
            return 0.0
        elif val < alpha:
            return 1.0

    # val within [alpha, gamma]
    if val <= beta:
        rng = 2.0 * (beta - alpha)
        if abs(rng) < 1e-9:
            return 0.5
        s1 = (val - alpha) / rng
        c = coeff * s1 ** power if not flipped else 1.0 - coeff * s1 ** power
    else:
        rng = 2.0 * (beta - gamma)
        if abs(rng) < 1e-9:
            return 0.5
        s1 = (val - gamma) / rng
        c = 1.0 - coeff * s1 ** power if not flipped else coeff * s1 ** power

    return max(0.0, min(1.0, c))


def conf_test(val, locut, hicut, power, midpt, nmval):
    """Vectorized conf_test — exact Fortran semantics."""
    val = np.asarray(val, dtype=np.float64)
    locut = np.asarray(locut, dtype=np.float64)
    hicut = np.asarray(hicut, dtype=np.float64)
    midpt = np.asarray(midpt, dtype=np.float64)

    coeff = 2.0 ** (power - 1.0)

    # Determine flipped state
    if np.isscalar(locut) and np.isscalar(hicut):
        flipped = locut > hicut
        alpha = min(locut, hicut) if not isinstance(locut, np.ndarray) else np.minimum(locut, hicut)
        gamma = max(locut, hicut) if not isinstance(locut, np.ndarray) else np.maximum(locut, hicut)
    else:
        alpha = np.minimum(locut, hicut)
        gamma = np.maximum(locut, hicut)
        flipped = locut > hicut

    beta = midpt

    # Initialize result
    c = np.full_like(val, 0.5, dtype=np.float64)

    # Beyond range — assign directly to avoid np.where order dependencies
    out_clear = ~flipped & (val > gamma)
    out_cloudy = ~flipped & (val < alpha)
    flip_cloudy = flipped & (val > gamma)
    flip_clear = flipped & (val < alpha)
    c[out_clear] = 1.0
    c[out_cloudy] = 0.0
    c[flip_cloudy] = 0.0
    c[flip_clear] = 1.0

    # Within range: compute S-curve (only on unclassified pixels)
    classified = out_clear | out_cloudy | flip_cloudy | flip_clear
    within = ~classified & (val >= alpha) & (val <= gamma)

    # Below midpt
    below = within & (val <= beta)
    rng_lo = 2.0 * (beta - alpha)
    s1_lo = np.where(np.abs(rng_lo) > 1e-9,
                     (val - alpha) / np.maximum(np.abs(rng_lo), 1e-9), 0.0)
    c[below] = np.where(~flipped[below] if isinstance(flipped, np.ndarray) else ~flipped,
                        coeff * s1_lo[below] ** power,
                        1.0 - coeff * s1_lo[below] ** power)

    # Above midpt
    above = within & (val > beta)
    rng_hi = 2.0 * (beta - gamma)
    s1_hi = np.where(np.abs(rng_hi) > 1e-9,
                     (val - gamma) / np.maximum(np.abs(rng_hi), 1e-9), 0.0)
    c[above] = np.where(~flipped[above] if isinstance(flipped, np.ndarray) else ~flipped,
                        1.0 - coeff * s1_hi[above] ** power,
                        coeff * s1_hi[above] ** power)

    return np.clip(c, 0.0, 1.0).astype(np.float32)


# ===========================================================================
# conf_test_2val — exact Fortran replica (conf_test_2val.f)
# ===========================================================================

def conf_test_2val(val, locut, hicut, power, midpt, nmval):
    """Two-sided confidence test for band ratio (exact Fortran logic)."""
    val = np.asarray(val, dtype=np.float64)
    alpha1, alpha2 = np.asarray(locut[0], dtype=np.float64), np.asarray(locut[1], dtype=np.float64)
    gamma1, gamma2 = np.asarray(hicut[0], dtype=np.float64), np.asarray(hicut[1], dtype=np.float64)
    beta1, beta2 = np.asarray(midpt[0], dtype=np.float64), np.asarray(midpt[1], dtype=np.float64)
    coeff = 2.0 ** (power - 1.0)

    c = np.full_like(val, 0.5, dtype=np.float64)

    if (alpha1 - gamma1) > 0.0:
        # Inner region fails test (cloudy) — Fortran if/else priority
        inner_fail = (val > alpha1) & (val < alpha2)
        c[inner_fail] = 0.0
        outer_pass = ~inner_fail & ((val < gamma1) | (val > gamma2))
        c[outer_pass] = 1.0

        # Remaining: in transition zones
        remaining = ~inner_fail & ~outer_pass

        # Lower set of limits: [gamma1, alpha1]
        lower = remaining & (val <= alpha1) & (val >= gamma1)
        lo_below = lower & (val >= beta1)
        lo_above = lower & (val < beta1)
        rng_lo = 2.0 * (beta1 - alpha1)
        s1_lo = np.where(np.abs(rng_lo) > 1e-9, (val - alpha1) / np.abs(rng_lo), 0.0)
        c[lo_below] = (coeff * s1_lo[lo_below] ** power)
        rng_lo2 = 2.0 * (beta1 - gamma1)
        s1_lo2 = np.where(np.abs(rng_lo2) > 1e-9, np.abs(val - gamma1) / np.abs(rng_lo2), 0.0)
        c[lo_above] = (1.0 - coeff * s1_lo2[lo_above] ** power)

        # Upper set of limits: [alpha2, gamma2]
        upper = remaining & (val >= alpha2) & (val <= gamma2)
        up_below = upper & (val <= beta2)
        up_above = upper & (val > beta2)
        rng_up = 2.0 * (beta2 - alpha2)
        s1_up = np.where(np.abs(rng_up) > 1e-9, (val - alpha2) / np.abs(rng_up), 0.0)
        c[up_below] = (coeff * s1_up[up_below] ** power)
        rng_up2 = 2.0 * (beta2 - gamma2)
        s1_up2 = np.where(np.abs(rng_up2) > 1e-9, (val - gamma2) / np.abs(rng_up2), 0.0)
        c[up_above] = (1.0 - coeff * s1_up2[up_above] ** power)
    else:
        # Inner region passes test (clear) — Fortran if/else priority
        inner_pass = (val > gamma1) & (val < gamma2)
        c[inner_pass] = 1.0
        outer_fail = ~inner_pass & ((val < alpha1) | (val > alpha2))
        c[outer_fail] = 0.0

        # Remaining: in transition zones
        remaining = ~inner_pass & ~outer_fail

        lower = remaining & (val <= gamma1) & (val >= alpha1)
        lo_below = lower & (val <= beta1)
        lo_above = lower & (val > beta1)
        rng_lo = 2.0 * (beta1 - alpha1)
        s1_lo = np.where(np.abs(rng_lo) > 1e-9, (val - alpha1) / np.abs(rng_lo), 0.0)
        c[lo_below] = (coeff * s1_lo[lo_below] ** power)
        rng_lo2 = np.abs(2.0 * (beta1 - gamma1))
        s1_lo2 = np.where(rng_lo2 > 1e-9, np.abs(val - gamma1) / rng_lo2, 0.0)
        c[lo_above] = (1.0 - coeff * s1_lo2[lo_above] ** power)

        upper = remaining & (val >= gamma2) & (val <= alpha2)
        up_below = upper & (val >= beta2)
        up_above = upper & (val < beta2)
        rng_up = 2.0 * (beta2 - alpha2)
        s1_up = np.where(np.abs(rng_up) > 1e-9, (val - alpha2) / np.abs(rng_up), 0.0)
        c[up_below] = (coeff * s1_up[up_below] ** power)
        rng_up2 = 2.0 * (beta2 - gamma2)
        s1_up2 = np.where(np.abs(rng_up2) > 1e-9, (val - gamma2) / np.abs(rng_up2), 0.0)
        c[up_above] = (1.0 - coeff * s1_up2[up_above] ** power)

    return np.clip(c, 0.0, 1.0).astype(np.float32)


# ===========================================================================
# tview — APOLLO 11-12um BTD LUT (tview.f)
# ===========================================================================

# LUT data from tview.f
_TVIEW_UTAB = np.array([2.00, 1.75, 1.5, 1.25, 1.00], dtype=np.float64)
_TVIEW_TTAB = np.array([190., 200., 210., 220., 230., 240., 250.,
                         260., 270., 280., 290., 300., 310.], dtype=np.float64)
_TVIEW_TAB = np.array([
    [0.559, 0.542, 0.520, 0.491, 0.450],
    [0.424, 0.416, 0.405, 0.391, 0.370],
    [0.286, 0.294, 0.305, 0.319, 0.340],
    [0.137, 0.162, 0.194, 0.238, 0.300],
    [0.123, 0.156, 0.199, 0.257, 0.340],
    [0.198, 0.240, 0.294, 0.367, 0.470],
    [0.333, 0.366, 0.409, 0.467, 0.550],
    [0.696, 0.704, 0.715, 0.729, 0.750],
    [1.217, 1.184, 1.141, 1.083, 1.000],
    [3.184, 2.926, 2.591, 2.140, 1.500],
    [5.178, 4.854, 4.433, 3.866, 3.060],
    [8.269, 7.885, 7.389, 6.720, 5.770],
    [12.452, 11.985, 11.381, 10.567, 9.410],
], dtype=np.float64).T  # (5, 13) — (u_idx, t_idx)


def tview(key: int, xmu: float, bt11: float) -> float:
    """APOLLO 11-12um BTD thin cirrus threshold LUT (exact Fortran replica).

    Args:
        key: 1=linear interpolation, 2=quadratic.
        xmu: secant of view zenith angle.
        bt11: 11um brightness temperature (K).

    Returns:
        11-12um BTD threshold for thin cirrus detection.
    """
    u = min(max(xmu, _TVIEW_UTAB[4]), _TVIEW_UTAB[0])
    t = min(max(bt11, _TVIEW_TTAB[0]), _TVIEW_TTAB[12])

    # Find u-index
    for i in range(1, 5):
        if u >= _TVIEW_UTAB[i]:
            if key == 1:
                i0, i1 = i - 1, i
            else:
                if i == 4:
                    i0, i1, i2 = i - 2, i - 1, i
                else:
                    i0, i1, i2 = i - 1, i, i + 1
            break

    # Find t-index
    for j in range(1, 13):
        if t <= _TVIEW_TTAB[j]:
            if key == 1:
                j0, j1 = j - 1, j
            else:
                if j == 12:
                    j0, j1, j2 = j - 2, j - 1, j
                else:
                    j0, j1, j2 = j - 1, j, j + 1
            break

    if key == 1:
        u0, u1 = _TVIEW_UTAB[i0], _TVIEW_UTAB[i1]
        t0, t1 = _TVIEW_TTAB[j0], _TVIEW_TTAB[j1]
        lu0 = (u - u1) / (u0 - u1)
        lu1 = (u - u0) / (u1 - u0)
        lt0 = (t - t1) / (t0 - t1)
        lt1 = (t - t0) / (t1 - t0)
        p0 = _TVIEW_TAB[i0, j0] * lu0 + _TVIEW_TAB[i1, j0] * lu1
        p1 = _TVIEW_TAB[i0, j1] * lu0 + _TVIEW_TAB[i1, j1] * lu1
        return p0 * lt0 + p1 * lt1
    else:
        u0, u1, u2 = _TVIEW_UTAB[i0], _TVIEW_UTAB[i1], _TVIEW_UTAB[i2]
        t0, t1, t2 = _TVIEW_TTAB[j0], _TVIEW_TTAB[j1], _TVIEW_TTAB[j2]
        lu0 = (u - u1) * (u - u2) / (u0 - u1) / (u0 - u2)
        lu1 = (u - u0) * (u - u2) / (u1 - u0) / (u1 - u2)
        lu2 = (u - u0) * (u - u1) / (u2 - u0) / (u2 - u1)
        lt0 = (t - t1) * (t - t2) / (t0 - t1) / (t0 - t2)
        lt1 = (t - t0) * (t - t2) / (t1 - t0) / (t1 - t2)
        lt2 = (t - t0) * (t - t1) / (t2 - t0) / (t2 - t1)
        p0 = (_TVIEW_TAB[i0, j0] * lu0 + _TVIEW_TAB[i1, j0] * lu1 +
              _TVIEW_TAB[i2, j0] * lu2)
        p1 = (_TVIEW_TAB[i0, j1] * lu0 + _TVIEW_TAB[i1, j1] * lu1 +
              _TVIEW_TAB[i2, j1] * lu2)
        p2 = (_TVIEW_TAB[i0, j2] * lu0 + _TVIEW_TAB[i1, j2] * lu1 +
              _TVIEW_TAB[i2, j2] * lu2)
        return p0 * lt0 + p1 * lt1 + p2 * lt2


def tview_array(key: int, xmu: np.ndarray, bt11: np.ndarray) -> np.ndarray:
    """Vectorized tview for arrays."""
    result = np.zeros_like(bt11, dtype=np.float64)
    for idx in np.ndindex(bt11.shape):
        result[idx] = tview(key, float(xmu[idx]), float(bt11[idx]))
    return result


# ===========================================================================
# Sun glint threshold interpolation (get_sg_thresholds.f90)
# ===========================================================================

def get_sg_thresholds(refang: float, snglnt_thr: Dict[str, list]) -> Tuple[float, float, float, float]:
    """Sun glint threshold interpolation (exact Fortran replica).

    Args:
        refang: Reflectance angle.
        snglnt_thr: Sun glint thresholds dict from YAML.

    Returns:
        (locut, hicut, midpt, power) tuple.
    """
    bounds = snglnt_thr["snglnt_bounds"]
    snglnt0 = snglnt_thr["snglnt0"]
    snglnt10 = snglnt_thr["snglnt10"]
    snglnt20 = snglnt_thr["snglnt20"]

    if refang <= bounds[1]:
        return snglnt0[0], snglnt0[2], snglnt0[1], snglnt0[3]

    if bounds[1] < refang <= bounds[2]:
        lo_ang, hi_ang = bounds[1], bounds[2]
        lo_val, hi_val = snglnt10[0], snglnt10[1]
        power = snglnt10[3]
        conf_range = snglnt10[2]
    else:  # bounds[2] < refang <= bounds[3]
        lo_ang, hi_ang = bounds[2], bounds[3]
        lo_val, hi_val = snglnt20[0], snglnt20[1]
        power = snglnt20[3]
        conf_range = snglnt20[2]

    a = (refang - lo_ang) / (hi_ang - lo_ang)
    midpt = lo_val + a * (hi_val - lo_val)
    hicut = midpt - conf_range
    locut = midpt + conf_range

    return locut, hicut, midpt, power


# ===========================================================================
# trispc — 8-11um clear-sky BTD regression
# ===========================================================================

def trispc(btd_11_12: np.ndarray) -> np.ndarray:
    """8-11 um clear-sky BTD regression from 11-12 um BTD."""
    x = np.asarray(btd_11_12, dtype=np.float64)
    return (2.7681 - 3.729 * x + 1.054 * x**2 - 0.102 * x**3).astype(np.float32)


# ===========================================================================
# Threshold loading
# ===========================================================================

def load_thresholds_yaml(yaml_path: Optional[str] = None) -> Dict:
    """Load cloud mask thresholds from YAML."""
    if yaml_path is None:
        project_root = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
        yaml_path = os.path.join(project_root, "coeff", "thresholds_mersi_ii.yaml")
    with open(yaml_path) as f:
        data = yaml.safe_load(f)
    return data["scenes"]


def get_thresholds(scene: str, yaml_path: Optional[str] = None) -> Dict[str, List[float]]:
    """Get thresholds for a specific scene type, including shared."""
    scenes = load_thresholds_yaml(yaml_path)
    result = {}
    if scene in scenes:
        result.update(scenes[scene].get("thresholds", {}))
    if "shared" in scenes:
        result.update(scenes["shared"].get("thresholds", {}))
    return result
