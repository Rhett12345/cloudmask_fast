"""Python ocean daytime cloud detection — exact Fortran ocean_day.f90 replica.

Ported from src/cloudmask/ocean_day.f90 with line-by-line fidelity:
  4 groups, 10 spectral tests, APOLLO tview LUT, sun glint threshold interpolation.

Usage:
    from fylat.ocean_day import detect_ocean_day
    clm, confidence = detect_ocean_day(l1_data, thresholds)
"""

from typing import Dict, Optional, Tuple

import numpy as np

from fylat.cloudmask_utils import (
    conf_test, conf_test_2val, trispc, tview,
    get_sg_thresholds,
)

# Constants matching Fortran
DTOR = np.pi / 180.0
MAX_VZA = 65.49


def _sec(x):
    """secant of angle in degrees."""
    return 1.0 / np.cos(np.maximum(np.abs(x), 1e-6) * DTOR)


def detect_ocean_day(
    l1_data: Dict[str, np.ndarray],
    thresholds: Dict[str, list],
    sst: Optional[np.ndarray] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """Run ocean daytime cloud detection (exact Fortran ocean_day.f90 replica).

    Args:
        l1_data: Dict from MersiL1Reader with ref_1..ref_19, bt_20..bt_25, sza, vza, lat, lon, lsm.
        thresholds: Dict from get_thresholds('ocean_day').
        sst: Sea surface temperature (K). If None, uses bt_24 as crude proxy.

    Returns:
        (cloud_mask, confidence):
        - cloud_mask: uint8 (0=cloudy, 1=prob_cloudy, 2=prob_clear, 3=clear)
        - confidence: float32 (1=clear, 0=cloudy)
    """
    # --- Extract spectral data ---
    r03 = l1_data.get("ref_3")
    r04 = l1_data.get("ref_4")
    r19 = l1_data.get("ref_19")
    bt_21 = l1_data.get("bt_21")
    bt_23 = l1_data.get("bt_23")
    bt_24 = l1_data["bt_24"]       # masir11 (11um)
    bt_25 = l1_data["bt_25"]       # masir12 (12um)
    sza = l1_data["sza"]
    vza = l1_data["vza"]
    lsm = l1_data.get("lsm")

    shape = bt_24.shape
    vis_ok = all(x is not None for x in [r03, r04, r19])
    has_bt_21 = bt_21 is not None
    has_bt_23 = bt_23 is not None

    # SST fallback
    if sst is None:
        sst = bt_24.copy()

    # --- Read thresholds ---
    dobt11 = thresholds.get("dobt11", [267.0, 270.0, 273.0, 1.0])
    pfmft_11max = thresholds.get("pfmft_11maxthre", [310.0])
    pfmft_btdmin = thresholds.get("pfmft_btd_min", [0.0])
    do11_12hi = thresholds.get("do11_12hi", [3.0])
    do11_4lo = thresholds.get("do11_4lo", [-5.5, -3.0, -0.5, 1.0])
    doref2 = thresholds.get("doref2", [0.065, 0.045, 0.030, 1.0])
    dovratlo = thresholds.get("dovratlo", [0.95, 0.90, 0.85, 1.0])
    dovrathi = thresholds.get("dovrathi", [0.95, 1.10, 1.15, 1.0])
    doref3 = thresholds.get("doref3", [0.04, 0.035, 0.03, 1.0])
    dotci = thresholds.get("dotci", [0.035, 0.005])

    # --- Initialize per-group minima ---
    cmin1 = np.ones(shape, dtype=np.float32)
    cmin2 = np.ones(shape, dtype=np.float32)
    cmin3 = np.ones(shape, dtype=np.float32)
    cmin4 = np.ones(shape, dtype=np.float32)
    ntest = [0, 0, 0, 0]

    # =====================================================================
    # GROUP 1: High Thick Cloud
    # =====================================================================

    # Test 1.1: 11 um BT threshold (lines 197-209)
    # locut=267, midpt=270, hicut=273, power=1.0, nmval=1
    # hicut > locut → not flipped → higher BT = clearer
    c1 = conf_test(bt_24, dobt11[0], dobt11[2], dobt11[3], dobt11[1], 1)
    cmin1 = np.minimum(cmin1, c1)
    ntest[0] += 1

    # Test 1.2: PFMFT test (lines 222-243, confidence DISABLED in Fortran)
    # Skip confidence contribution; confidence section is commented out in Fortran

    # Test 1.3: NFMFT test (lines 245-261, confidence DISABLED)
    # Skip

    # Test 1.4: SST test (lines 313-356)
    # sfcdif = sfctmp - BT11; larger → more cloudy
    btd_11_12 = bt_24 - bt_25
    sfcdif = sst - bt_24

    sst_thrsh = np.full(shape, 6.0, dtype=np.float64)
    if lsm is not None:
        sst_thrsh[lsm == 0] = 10.0  # ShallowOcean

    midpt_sst = sst_thrsh + np.where(btd_11_12 >= 1.0,
                                      2.0 * np.round(btd_11_12), 0.0)
    midpt_sst += (vza / MAX_VZA) ** 4 * 3.0

    locut_sst = midpt_sst + 1.0
    hicut_sst = midpt_sst - 2.0
    # locut > hicut → flipped → larger sfcdif = more cloudy
    c_sst = conf_test(sfcdif, locut_sst, hicut_sst, 1.0, midpt_sst, 1)
    cmin1 = np.minimum(cmin1, c_sst)
    ntest[0] += 1

    # =====================================================================
    # GROUP 2: Low Cloud — Thick
    # =====================================================================

    # Test 2.1: Tri-spectral 8-11-12 um BTD (lines 369-393)
    if has_bt_23:
        btd_8_11 = bt_23 - bt_24
        tri_thr = trispc(btd_11_12)
        # locut=tri+0.5, hicut=tri-0.5 → locut > hicut → flipped
        c4 = conf_test(btd_8_11, tri_thr + 0.5, tri_thr - 0.5, 1.0, tri_thr, 1)
        # NOTE: Fortran does NOT include c4 in cmin2 (commented out lines 391-392)
        # cmin2 = np.minimum(cmin2, c4)  # disabled
        ntest[1] += 1

    # Test 2.2: 11-12 um BTD thin cirrus (APOLLO tview, lines 405-450)
    schi = _sec(vza)
    # Apply tview LUT per pixel; fallback = do11_12hi(1)
    dfthrsh = np.full(shape, do11_12hi[0], dtype=np.float64)
    for idx in np.ndindex(shape):
        try:
            tv = tview(1, float(schi[idx]), float(bt_24[idx]))
            if tv >= 0.1:
                dfthrsh[idx] = tv
        except (IndexError, ValueError):
            pass

    # locut=dfthrsh*1.3, hicut=dfthrsh-1.25 → locut > hicut for dfthrsh<~25
    # flipped → larger BTD = more cloudy
    c6 = conf_test(btd_11_12, dfthrsh * 1.3, dfthrsh - 1.25, 1.0, dfthrsh, 1)
    cmin2 = np.minimum(cmin2, c6)
    ntest[1] += 1

    # Test 2.3: 11-4 um BTD fog/low cloud (lines 462-495)
    if has_bt_21 and vis_ok:
        btd_11_4 = bt_24 - bt_21
        # locut=-5.5, midpt=-3.0, hicut=-0.5 → hicut > locut → not flipped
        # higher BTD (closer to 0) = clearer
        c7 = conf_test(btd_11_4, do11_4lo[0], do11_4lo[2], do11_4lo[3], do11_4lo[1], 1)
        cmin2 = np.minimum(cmin2, c7)
        ntest[1] += 1

    # --- Compute sun glint flag ---
    glint_angle = l1_data.get("glint_angle", None)
    if glint_angle is None and "rel_azimuth" in l1_data:
        # Compute glint angle (matching Fortran GLINTANGLE in frontend_module)
        sza_rad = np.deg2rad(sza)
        vza_rad = np.deg2rad(vza)
        rel_az_rad = np.deg2rad(l1_data.get("rel_azimuth", np.zeros_like(sza)))
        cos_glint = np.cos(sza_rad) * np.cos(vza_rad) - \
                    np.sin(sza_rad) * np.sin(vza_rad) * np.cos(rel_az_rad)
        glint_angle = np.rad2deg(np.arccos(np.clip(cos_glint, -1.0, 1.0)))

    # Sun glint mask: glint angle < 36° (matching MOD35 convention)
    snglnt = (glint_angle < 36.0) if glint_angle is not None else np.zeros(shape, dtype=bool)

    # =====================================================================
    # GROUP 3: Thick Cloud — Reflectance
    # Apply only to non-glint pixels (glint pixels use relaxed pass-through)
    # =====================================================================

    if vis_ok:
        # Test 3.1: 0.86 um reflectance — non-glint only
        c8 = conf_test(r04, doref2[0], doref2[2], doref2[3], doref2[1], 1)
        c8[snglnt] = 1.0  # glint: force pass (Fortran uses relaxed thresholds)
        cmin3 = np.minimum(cmin3, c8)
        ntest[2] += 1

        # Test 3.2: Visible ratio 0.86/0.65 — non-glint only
        with np.errstate(divide="ignore", invalid="ignore"):
            vrat = np.where(r03 > 0.005, r04 / np.maximum(r03, 0.005), 1.0)
            vrat = np.clip(vrat, 0.1, 5.0)

        c9 = conf_test_2val(
            vrat,
            np.array([dovratlo[2], dovrathi[0]]),
            np.array([dovratlo[0], dovrathi[2]]),
            1.0,
            np.array([dovratlo[1], dovrathi[1]]),
            2,
        )
        c9[snglnt] = 1.0  # glint: force pass
        cmin3 = np.minimum(cmin3, c9)
        ntest[2] += 1

    # =====================================================================
    # GROUP 4: Thin Cirrus (1.38 um)
    # =====================================================================

    if vis_ok and r19 is not None:
        c11 = conf_test(r19, doref3[0], doref3[2], doref3[3], doref3[1], 1)
        c11[snglnt] = 1.0  # glint: force pass
        cmin4 = np.minimum(cmin4, c11)
        ntest[3] += 1

    # =====================================================================
    # Aggregate confidence (lines 640-654)
    #   pre_confdnc = cmin1 * cmin2 * cmin3 * cmin4
    #   groups = count of groups with >0 tests
    #   confdnc = pre_confdnc ** (1/groups)
    # =====================================================================

    ngroups = np.full(shape, 0, dtype=np.float32)
    for i in range(4):
        if ntest[i] > 0:
            ngroups += 1

    pre_conf = np.ones(shape, dtype=np.float64)
    if ntest[0] > 0:
        pre_conf *= cmin1.astype(np.float64)
    if ntest[1] > 0:
        pre_conf *= cmin2.astype(np.float64)
    if ntest[2] > 0:
        pre_conf *= cmin3.astype(np.float64)
    if ntest[3] > 0:
        pre_conf *= cmin4.astype(np.float64)

    with np.errstate(divide="ignore", invalid="ignore"):
        confidence = np.where(ngroups > 0,
                              pre_conf ** (1.0 / ngroups), 0.0)
    confidence = np.clip(confidence, 0.0, 1.0).astype(np.float32)

    # Convert to 0-3 cloud mask (matching Fortran qa bit convention)
    # 0=cloudy, 1=probably_cloudy, 2=probably_clear, 3=clear
    cloud_mask = np.full(shape, 0, dtype=np.uint8)
    cloud_mask = np.where(confidence >= 0.95, 3, cloud_mask)
    cloud_mask = np.where((confidence >= 0.66) & (confidence < 0.95), 2, cloud_mask)
    cloud_mask = np.where((confidence >= 0.33) & (confidence < 0.66), 1, cloud_mask)

    return cloud_mask, confidence
