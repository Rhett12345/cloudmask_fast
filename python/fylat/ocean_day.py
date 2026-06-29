"""Python implementation of ocean daytime cloud detection (Fortran ocean_day.f90).

Uses YAML thresholds and replicates the 4-group decision tree structure.
This is the Phase 2 prototype for validating the YAML threshold approach.

Usage:
    from fylat.ocean_day import detect_ocean_day
    clm, confidence = detect_ocean_day(l1_data, thresholds)
"""

from typing import Dict, Optional, Tuple

import numpy as np

from fylat.cloudmask_utils import conf_test, conf_test_2val, trispc


def detect_ocean_day(
    l1_data: Dict[str, np.ndarray],
    thresholds: Dict[str, list],
    sst: Optional[np.ndarray] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """Run ocean daytime cloud detection.

    Args:
        l1_data: Dict from MersiL1Reader.read_all() with keys:
            ref_1..ref_19, bt_20..bt_25, sza, vza, lat, lon, lsm.
        thresholds: Dict from cloudmask_utils.get_thresholds('ocean_day').
        sst: Optional sea surface temperature field (K). If None, uses
             bt_24 as proxy.

    Returns:
        (cloud_mask, confidence) tuple:
        - cloud_mask: uint8 array (0=cloud, 1=probably_cloudy,
                       2=probably_clear, 3=clear)
        - confidence: float32 array (0=cloudy, 1=clear)
    """
    # --- Extract spectral data ---
    r03 = l1_data.get("ref_3", None)     # 0.65 um
    r04 = l1_data.get("ref_4", None)     # 0.86 um
    r19 = l1_data.get("ref_19", None)    # 1.38 um
    bt_21 = l1_data.get("bt_21", None)   # 4.05 um
    bt_23 = l1_data.get("bt_23", None)   # 8.5 um
    bt_24 = l1_data["bt_24"]             # 10.8 um (masir11)
    bt_25 = l1_data["bt_25"]             # 12.0 um (masir12)
    sza = l1_data["sza"]
    vza = l1_data["vza"]

    shape = bt_24.shape
    vis_ok = (r03 is not None) and (r04 is not None) and (r19 is not None)

    # SST fallback
    if sst is None:
        sst = bt_24.copy()  # crude fallback

    # --- Initialize confidence ---
    # Per-group minimum confidences
    cmin1 = np.ones(shape, dtype=np.float32)  # Group 1
    cmin2 = np.ones(shape, dtype=np.float32)  # Group 2
    cmin3 = np.ones(shape, dtype=np.float32)  # Group 3
    cmin4 = np.ones(shape, dtype=np.float32)  # Group 4

    ntest1 = ntest2 = ntest3 = ntest4 = 0

    # =====================================================================
    # GROUP 1: High Thick Cloud
    # =====================================================================

    # Test 1.1: 11 um BT threshold
    dobt11 = thresholds.get("dobt11", [267.0, 270.0, 273.0, 1.0])
    c1 = conf_test(bt_24, dobt11[0], dobt11[2], dobt11[3], dobt11[1], 1)
    cmin1 = np.minimum(cmin1, c1)
    ntest1 += 1

    # Test 1.2: SST test (sfc_temp - BT11)
    # Dynamic thresholds simplified to scalar for Python prototype
    sst_thrsh_land = 10.0
    sst_thrsh_ocean = 6.0

    # Detect shallow water pixels using LSM if available
    lsm = l1_data.get("lsm", None)
    sst_thrsh_arr = np.full(shape, sst_thrsh_ocean, dtype=np.float64)
    if lsm is not None:
        sst_thrsh_arr[lsm == 0] = sst_thrsh_land  # ShallowOcean

    btd_11_12 = bt_24 - bt_25
    sfcdif_arr = sst - bt_24

    # Per-pixel SST test with dynamic threshold
    sst_midpt = sst_thrsh_arr + 2.0 * np.where(btd_11_12 >= 1.0, np.round(btd_11_12), 0)
    sst_midpt += (vza / 65.49) ** 4 * 3.0

    # Simplified confidence: sfcdif < midpt = clear, >= midpt+3 = cloudy
    sst_diff = sfcdif_arr - sst_midpt
    # Scale to [0,1]: negative diff (clear) -> high confidence, positive diff > 3 -> low confidence
    c_sst = np.clip(1.0 - sst_diff / 3.0, 0.0, 1.0).astype(np.float32)
    cmin1 = np.minimum(cmin1, c_sst)
    ntest1 += 1

    # =====================================================================
    # GROUP 2: Low Cloud — Thick
    # =====================================================================

    # Test 2.1: Tri-spectral 8-11-12 um BTD
    if bt_23 is not None:
        btd_8_11 = bt_23 - bt_24
        tri_thr = trispc(btd_11_12)
        c4 = conf_test(btd_8_11, tri_thr + 0.5, tri_thr - 0.5, 1.0, tri_thr, 1)
        cmin2 = np.minimum(cmin2, c4)
        ntest2 += 1

    # Test 2.2: 11-12 um BTD thin cirrus
    # Fallback threshold (without APOLLO tview LUT)
    do11_12hi = thresholds.get("do11_12hi", [3.0])
    dfthrsh_val = float(do11_12hi[0])
    c6 = conf_test(btd_11_12, dfthrsh_val * 1.3, dfthrsh_val - 1.25, 1.0, dfthrsh_val, 1)
    cmin2 = np.minimum(cmin2, c6)
    ntest2 += 1

    # Test 2.3: 11-4 um BTD fog/low cloud (visible only, no glint)
    if bt_21 is not None and vis_ok:
        do11_4lo = thresholds.get("do11_4lo", [-5.5, -3.0, -0.5, 1.0])
        btd_11_4 = bt_24 - bt_21
        # Skip sun glint regions (simplified: skip high SZA where glint likely)
        c7 = conf_test(btd_11_4, do11_4lo[0], do11_4lo[2], do11_4lo[3], do11_4lo[1], 1)
        cmin2 = np.minimum(cmin2, c7)
        ntest2 += 1

    # =====================================================================
    # GROUP 3: Reflectance Tests (visible only)
    # =====================================================================

    if vis_ok:
        # Test 3.1: 0.86 um reflectance
        doref2 = thresholds.get("doref2", [0.065, 0.045, 0.030, 1.0])
        c8 = conf_test(r04, doref2[0], doref2[2], doref2[3], doref2[1], -1)
        cmin3 = np.minimum(cmin3, c8)
        ntest3 += 1

        # Test 3.2: Visible ratio 0.86/0.65
        with np.errstate(divide="ignore", invalid="ignore"):
            vrat = np.where(r03 > 0.01, r04 / r03, 1.0)
            vrat = np.clip(vrat, 0.1, 5.0)

        dovrathi = thresholds.get("dovrathi", [0.95, 1.10, 1.15, 1.0])
        dovratlo = thresholds.get("dovratlo", [0.95, 0.90, 0.85, 1.0])
        c9 = conf_test_2val(
            vrat,
            np.array([dovratlo[0], dovrathi[0]]),
            np.array([dovratlo[2], dovrathi[2]]),
            1.0,
            np.array([dovratlo[1], dovrathi[1]]),
            2,
        )
        cmin3 = np.minimum(cmin3, c9)
        ntest3 += 1

    # =====================================================================
    # GROUP 4: Thin Cirrus
    # =====================================================================

    if vis_ok and r19 is not None:
        doref3 = thresholds.get("doref3", [0.04, 0.035, 0.03, 1.0])
        c11 = conf_test(r19, doref3[0], doref3[2], doref3[3], doref3[1], -1)
        cmin4 = np.minimum(cmin4, c11)
        ntest4 += 1

    # =====================================================================
    # Aggregate confidence
    # =====================================================================

    # Geometric mean of group minima (matching Fortran logic)
    groups_active = np.zeros(shape, dtype=np.float32)
    if ntest1 > 0:
        groups_active += 1
    if ntest2 > 0:
        groups_active += 1
    if ntest3 > 0:
        groups_active += 1
    if ntest4 > 0:
        groups_active += 1

    pre_conf = np.ones(shape, dtype=np.float32)
    if ntest1 > 0:
        pre_conf = pre_conf * cmin1
    if ntest2 > 0:
        pre_conf = pre_conf * cmin2
    if ntest3 > 0:
        pre_conf = pre_conf * cmin3
    if ntest4 > 0:
        pre_conf = pre_conf * cmin4

    with np.errstate(divide="ignore", invalid="ignore"):
        confidence = np.where(
            groups_active > 0,
            pre_conf ** (1.0 / groups_active),
            0.0,
        )
    confidence = np.clip(confidence, 0.0, 1.0).astype(np.float32)

    # Convert to 0-3 cloud mask
    cloud_mask = np.full(shape, 2, dtype=np.uint8)  # default: probably clear
    cloud_mask = np.where(confidence >= 0.95, 3, cloud_mask)  # clear
    cloud_mask = np.where(confidence >= 0.66, 2, cloud_mask)  # probably clear
    cloud_mask = np.where(confidence >= 0.33, 1, cloud_mask)  # probably cloudy
    cloud_mask = np.where(confidence < 0.33, 0, cloud_mask)   # cloudy

    return cloud_mask, confidence
