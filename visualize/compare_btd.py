"""
compare_btd.py — Compare MERSI BT/BTD distributions against MYD35 truth
========================================================================
For each scene stratum, computes MERSI BT11 and key BTD (11-12, 8-11, 7-11,
11-4) distributions separately for MYD35-cloudy and MYD35-clear pixels.
Finds optimal separation thresholds and compares with current values.

Usage:
  python compare_btd.py --data_dir /data/Data_yuq/fy3_cloud/20220803/
"""

from __future__ import annotations
import argparse
import os
import re
from pathlib import Path
from datetime import datetime, timezone

import numpy as np
import h5py

from io_mersi import (
    load_clm_hdf5, load_mersi_bt108, find_l1b_for_clm, parse_mersi_datetime,
)
from io_myd35 import load_best_myd35_for_mersi


# ─────────────────────────────────────────────────────────────────────────────
# BT/BTD calibration for all emissive channels
# ─────────────────────────────────────────────────────────────────────────────

C1 = 1.1910427e-5
C2 = 1.4387752

# EV_1KM_Emissive: idx -> (lambda_cm, cal_idx)
_EV1KM = {
    "BT8 (8.55)":  (0, 8.55e-4,  2),
    "BT11 (10.8)": (1, 10.80e-4, 3),
    "BT12 (12.0)": (2, 12.00e-4, 4),
}

_EV250 = {
    "BT4 (4.05)":  (1, 4.05e-4,  1),
}

def _cal_bt(dn, cal_avg, lam_cm):
    c0, c1, c2 = float(cal_avg[0]), float(cal_avg[1]), float(cal_avg[2])
    rad_wn = c0 + c1 * dn + c2 * dn * dn
    rad_wn = np.where(rad_wn > 0, rad_wn, np.nan)
    rad_wl = rad_wn / (lam_cm ** 2)
    bt = C2 / (lam_cm * np.log(C1 / (lam_cm ** 5 * rad_wl) + 1.0))
    return np.where((bt > 170) & (bt < 340), bt, np.nan)


def load_all_bt(l1b_path: str) -> dict[str, np.ndarray]:
    """Load and calibrate all IR brightness temperatures."""
    result = {}
    with h5py.File(l1b_path, "r") as f:
        ir_cal = f["Calibration/IR_Cal_Coeff"][:]  # (6,4,200)
        cal_avg_all = ir_cal.mean(axis=2)  # average over scans

        if "Data/EV_1KM_Emissive" in f:
            ev = f["Data/EV_1KM_Emissive"][:].astype(np.float64)
            for name, (idx, lam, cal_idx) in _EV1KM.items():
                result[name] = _cal_bt(ev[idx], cal_avg_all[cal_idx], lam)

        if "Data/EV_250_Aggr.1KM_Emissive" in f:
            ev250 = f["Data/EV_250_Aggr.1KM_Emissive"][:].astype(np.float64)
            for name, (idx, lam, cal_idx) in _EV250.items():
                result[name] = _cal_bt(ev250[idx], cal_avg_all[cal_idx], lam)

    return result


# ─────────────────────────────────────────────────────────────────────────────
# Optimal threshold finder
# ─────────────────────────────────────────────────────────────────────────────

def find_optimal_threshold(
    values: np.ndarray,
    y_true: np.ndarray,   # 0=cloudy (MYD35 cls 0-1), 1=clear (MYD35 cls 2-3)
    direction: str = "lt",  # "lt" = cloudy if value < thresh, "gt" = cloudy if value > thresh
) -> dict:
    """
    Find threshold that maximizes HSS for a given BTD test.

    Returns dict with optimal threshold, max HSS, and distribution stats.
    """
    valid = np.isfinite(values) & np.isfinite(y_true)
    vals = values[valid]
    y = y_true[valid]

    if len(vals) < 100:
        return {}

    # Distribution stats
    cloudy_vals = vals[y == 0]
    clear_vals  = vals[y == 1]

    if len(cloudy_vals) < 10 or len(clear_vals) < 10:
        return {}

    # Scan thresholds
    candidates = np.linspace(np.nanpercentile(vals, 2), np.nanpercentile(vals, 98), 100)
    best_hss = -1
    best_thr = 0
    best_pod = 0
    best_far = 0

    for thr in candidates:
        if direction == "lt":
            pred_cloud = vals < thr
        else:
            pred_cloud = vals > thr

        y_cloud = (y == 0)

        tp = (pred_cloud & y_cloud).sum()
        fp = (pred_cloud & ~y_cloud).sum()
        fn = (~pred_cloud & y_cloud).sum()
        tn = (~pred_cloud & ~y_cloud).sum()

        n = len(vals)
        pod = tp / (tp + fn + 1e-9)
        far = fp / (tp + fp + 1e-9)
        expected = ((tp + fp) * (tp + fn) + (tn + fp) * (tn + fn)) / (n + 1e-9)
        hss = (tp + tn - expected) / (n - expected + 1e-9)

        if hss > best_hss:
            best_hss = hss
            best_thr = thr
            best_pod = pod
            best_far = far

    return {
        "threshold": float(best_thr),
        "hss": float(best_hss),
        "pod": float(best_pod),
        "far": float(best_far),
        "n": len(vals),
        "n_cloudy": len(cloudy_vals),
        "n_clear": len(clear_vals),
        "cloudy_mean": float(np.mean(cloudy_vals)),
        "cloudy_std": float(np.std(cloudy_vals)),
        "clear_mean": float(np.mean(clear_vals)),
        "clear_std": float(np.std(clear_vals)),
        "separation": float(abs(np.mean(cloudy_vals) - np.mean(clear_vals)) /
                           max(np.std(vals), 0.01)),
    }


# ─────────────────────────────────────────────────────────────────────────────
# Scene-level comparison
# ─────────────────────────────────────────────────────────────────────────────

def compare_orbit(
    recal_path: str,
    myd35_dirs: list[str],
    mersi_root: str = "/data/Data_yuq/mersi",
    time_window_min: int = 15,
    min_overlap: float = 0.05,
) -> dict[str, dict]:
    """Compare BTD distributions for key scenes in one orbit."""
    recal_data = load_clm_hdf5(recal_path)
    if recal_data is None:
        return {}

    lat = recal_data["lat"]
    lon = recal_data["lon"]

    # Load GEO
    l1b_path = find_l1b_for_clm(recal_path, mersi_root)
    lsm = None
    if l1b_path:
        geo_path = l1b_path.replace("_1000M_MS.HDF", "_GEO1K_MS.HDF")
        try:
            with h5py.File(geo_path, "r") as f:
                lsm = f["Geolocation/LandSeaMask"][:].astype(np.int32)
        except Exception:
            pass
    if lsm is None:
        lsm = np.where(np.abs(lat) > 80, 1, 0).astype(np.int32)

    # Load all BT
    bt_all = {}
    if l1b_path:
        bt_all = load_all_bt(l1b_path)

    # Load MYD35
    mersi_dt = parse_mersi_datetime(recal_path)
    if mersi_dt is None:
        return {}

    myd35_data = load_best_myd35_for_mersi(
        mersi_lat=lat, mersi_lon=lon, mersi_dt=mersi_dt,
        search_dirs=myd35_dirs,
        time_window_min=time_window_min, min_overlap=min_overlap,
    )
    if myd35_data is None:
        return {}

    myd_clm = myd35_data["clm_resampled"]
    myd_cloudy = (myd_clm <= 1)  # MYD35 classes 0-1 = cloud
    myd_clear  = (myd_clm >= 2)  # MYD35 classes 2-3 = clear

    # Stratify by BT range
    bt11 = bt_all.get("BT11 (10.8)")
    if bt11 is None:
        return {}

    bt_ranges = {
        "Cold(<230K)": (bt11 < 230),
        "Cool(230-250K)": (bt11 >= 230) & (bt11 < 250),
        "Mod(250-270K)": (bt11 >= 250) & (bt11 < 270),
        "Warm(>270K)": (bt11 >= 270),
    }

    # Surface types
    surfaces = {
        "Ocean": (lsm == 0),
        "Land": (lsm == 1),
    }

    # BTD definitions and current thresholds
    btds = {
        "BTD_11_12": {
            "formula": lambda bt: bt["BT11 (10.8)"] - bt["BT12 (12.0)"],
            "direction": "lt",  # lower BTD = thinner cirrus? actually depends on scene
            "current_thr_polar_land": 3.0,  # pnl11_12hi
            "current_thr_ocean": 3.0,
            "note": "thin cirrus: lower BTD -> cloud? At cold T, sign reverses",
        },
        "BTD_8_11": {
            "formula": lambda bt: bt.get("BT8 (8.55)", bt["BT11 (10.8)"]*0+np.nan) - bt["BT11 (10.8)"],
            "direction": "lt",  # lower BTD (more water vapor) -> clear
            "current_thr_ocean": 24.0,  # pno86_73 (but for 8-7.3, not 8-11)
            "note": "water vapor: more negative = more WV = clear",
        },
        "BTD_11_4": {
            "formula": lambda bt: bt["BT11 (10.8)"] - (bt.get("BT4 (4.05)", bt["BT11 (10.8)"]*0+np.nan) - 2.0),
            "direction": "lt",  # lower = low cloud/fog
            "current_thr_ocean_night": 1.0,  # no11_4lo
            "current_thr_polar_night": -0.2,  # pn_11_4l(2)
            "note": "fog/low cloud: 4um dead at cold T",
        },
        "BT11_raw": {
            "formula": lambda bt: bt["BT11 (10.8)"],
            "direction": "lt",  # lower BT = high cloud (colder)
            "current_thr_polar": 270.0,  # pnobt11(2)
            "note": "BT threshold: colder -> cloud",
        },
    }

    orbit_tag = re.search(r'(\d{8}_\d{4})', os.path.basename(recal_path))
    orbit_tag = orbit_tag.group(1) if orbit_tag else "unknown"

    results = {}
    for bt_name, bt_range_mask in bt_ranges.items():
        for sfc_name, sfc_mask in surfaces.items():
            base = bt_range_mask & sfc_mask & (myd_clm >= 0)
            if base.sum() < 500:
                continue

            y_true = np.where(myd_cloudy, 0, 1).astype(np.float64)

            for btd_name, btd_def in btds.items():
                vals = btd_def["formula"](bt_all)
                key = f"{orbit_tag}|{bt_name}|{sfc_name}|{btd_name}"

                opt = find_optimal_threshold(vals[base], y_true[base],
                                             direction=btd_def["direction"])
                if opt:
                    opt["orbit"] = orbit_tag
                    opt["bt_range"] = bt_name
                    opt["surface"] = sfc_name
                    opt["btd"] = btd_name
                    opt["note"] = btd_def.get("note", "")
                    results[key] = opt

    return results


# ─────────────────────────────────────────────────────────────────────────────
# Reporter
# ─────────────────────────────────────────────────────────────────────────────

def print_btd_report(all_results: dict[str, dict]) -> None:
    """Print BTD comparison report grouped by scene."""
    # Group by BT range + surface
    groups: dict[str, list] = {}
    for key, r in all_results.items():
        scene = f"{r['bt_range']}_{r['surface']}"
        if scene not in groups:
            groups[scene] = []
        groups[scene].append(r)

    for scene in sorted(groups.keys()):
        entries = groups[scene]
        # Compute average N from BT11 test
        n_avg = int(np.median([e["n"] for e in entries if "BT11_raw" in e["btd"]] or [0]))

        print(f"\n{'='*80}")
        print(f"  {scene}  (N~{n_avg:,})")
        print(f"{'='*80}")
        print(f"  {'BTD':<14s} {'Cloudy_mean':>10s} {'Clear_mean':>10s} "
              f"{'Sep':>6s} {'OptThr':>8s} {'MaxHSS':>7s} {'POD':>7s} {'FAR':>7s}")
        print(f"  {'':-<14s} {'':->10s} {'':->10s} {'':->6s} {'':->8s} {'':->7s} {'':->7s} {'':->7s}")

        for e in entries:
            print(f"  {e['btd']:<14s} {e['cloudy_mean']:10.2f} {e['clear_mean']:10.2f} "
                  f"{e['separation']:5.2f}σ {e['threshold']:8.2f} {e['hss']:7.4f} "
                  f"{e['pod']*100:6.1f}% {e['far']*100:6.1f}%")


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compare MERSI BT/BTD distributions against MYD35 truth")
    parser.add_argument("--data_dir", required=True,
                        help="Directory with RECALI CLM files")
    parser.add_argument("--myd35_dir", nargs="+",
                        default=["/data/Data_yuq/aqua_modis/MYD35_L2/"],
                        help="MYD35 search directories")
    parser.add_argument("--mersi_root", default="/data/Data_yuq/mersi")
    parser.add_argument("--time_window", type=int, default=15)
    parser.add_argument("--min_overlap", type=float, default=0.05)
    args = parser.parse_args()

    data_dir = Path(args.data_dir)
    recal_files = sorted(data_dir.glob("*_RECALI.HDF"))

    all_results = {}
    for recal_path in recal_files:
        print(f"\n[COMPARE] {recal_path.name}")
        results = compare_orbit(
            str(recal_path),
            myd35_dirs=args.myd35_dir,
            mersi_root=args.mersi_root,
            time_window_min=args.time_window,
            min_overlap=args.min_overlap,
        )
        all_results.update(results)

    if not all_results:
        print("[ERROR] No comparison results.")
        exit(1)

    print_btd_report(all_results)
