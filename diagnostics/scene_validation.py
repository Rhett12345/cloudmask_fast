"""
scene_validation.py — Stratified cloud mask validation by scene type
====================================================================
Splits MERSI-MYD35 overlap pixels into strata by surface type x day/night
x latitude band x BT11 range, then computes per-stratum validation metrics.

Output: formatted table showing which scenes need the most improvement.

Usage:
  python scene_validation.py --data_dir /data/Data_yuq/fy3_cloud/20220803/

Design:
  - Reuses io_mersi / io_myd35 for data loading
  - Reuses figure_2.compute_validation_stats for metrics
  - Classifies pixels using LandSeaMask, SolarZenith, Latitude, BT11
"""

from __future__ import annotations
import argparse
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "visualize"))
import re
from pathlib import Path
from datetime import datetime, timezone

import numpy as np
import h5py

from io_mersi import (
    load_clm_hdf5, load_mersi_bt108, find_l1b_for_clm,
    parse_mersi_datetime,
)
from io_myd35 import load_best_myd35_for_mersi
from figure_2 import compute_validation_stats


# ─────────────────────────────────────────────────────────────────────────────
# Scene classification
# ─────────────────────────────────────────────────────────────────────────────

def classify_pixels(
    lat: np.ndarray,
    lsm: np.ndarray,
    sza: np.ndarray,
    bt11: np.ndarray | None,
) -> dict:
    """
    Return boolean masks for each scene stratum.

    Returns
    -------
    dict with keys like 'Ocean_Night_Polar_Cold' → boolean mask
    """
    valid = np.isfinite(lat) & np.isfinite(lsm)

    # Surface type
    ocean = valid & (lsm == 0)
    land  = valid & (lsm == 1)
    inland = valid & (lsm == 2)

    # Day/night
    sza_valid = np.isfinite(sza)
    day   = sza_valid & (sza < 85.0)
    night = sza_valid & (sza >= 85.0)

    # Latitude band
    abs_lat = np.abs(lat)
    tropical = valid & (abs_lat < 30)
    midlat   = valid & (abs_lat >= 30) & (abs_lat < 60)
    polar    = valid & (abs_lat >= 60)

    # BT11 range
    if bt11 is not None:
        bt_valid = np.isfinite(bt11)
        bt_cold   = bt_valid & (bt11 < 230)
        bt_cool   = bt_valid & (bt11 >= 230) & (bt11 < 250)
        bt_mod    = bt_valid & (bt11 >= 250) & (bt11 < 270)
        bt_warm   = bt_valid & (bt11 >= 270)
        bt_ranges = {
            "Cold(<230K)": bt_cold,
            "Cool(230-250K)": bt_cool,
            "Mod(250-270K)": bt_mod,
            "Warm(>270K)": bt_warm,
        }
    else:
        bt_ranges = {"AllBT": np.ones(lat.shape, dtype=bool)}

    strata = {}
    for sfc_name, sfc_mask in [("Ocean", ocean), ("Land", land),
                                 ("InlandWater", inland)]:
        for dn_name, dn_mask in [("Day", day), ("Night", night)]:
            for lat_name, lat_mask in [("Tropical", tropical),
                                        ("MidLat", midlat),
                                        ("Polar", polar)]:
                base = sfc_mask & dn_mask & lat_mask
                if not base.any():
                    continue
                for bt_name, bt_mask in bt_ranges.items():
                    mask = base & bt_mask
                    if mask.sum() > 100:  # minimum sample threshold
                        key = f"{sfc_name}_{dn_name}_{lat_name}_{bt_name}"
                        strata[key] = mask

    return strata


# ─────────────────────────────────────────────────────────────────────────────
# Single orbit analysis
# ─────────────────────────────────────────────────────────────────────────────

def validate_orbit(
    recal_path: str,
    onboard_path: str,
    myd35_dirs: list[str],
    mersi_root: str = "/data/Data_yuq/mersi",
    time_window_min: int = 15,
    min_overlap: float = 0.05,
) -> list[dict]:
    """
    Run stratified validation for one orbit.

    Returns list of per-stratum stat dicts, sorted by HSS (worst first).
    """
    # Load MERSI CLM
    recal_data = load_clm_hdf5(recal_path)
    onboard_data = load_clm_hdf5(onboard_path)
    if recal_data is None or onboard_data is None:
        return []

    lat = recal_data["lat"]
    lon = recal_data["lon"]
    recal_clm = recal_data["clm"]
    onboard_clm = onboard_data["clm"]

    # Load GEO ancillary
    l1b_path = find_l1b_for_clm(recal_path, mersi_root)
    lsm = None
    sza = None
    if l1b_path:
        geo_path = l1b_path.replace("_1000M_MS.HDF", "_GEO1K_MS.HDF")
        try:
            with h5py.File(geo_path, "r") as f:
                lsm = f["Geolocation/LandSeaMask"][:].astype(np.int32)
                sza = f["Geolocation/SolarZenith"][:].astype(np.float64)
        except Exception:
            pass

    if lsm is None:
        # Fallback: rough classification from latitude alone
        lsm = np.where(np.abs(lat) > 80, 1, 0).astype(np.int32)
    if sza is None:
        sza = np.full_like(lat, 85.0)  # assume night

    # Load BT11
    bt11 = None
    if l1b_path:
        bt11 = load_mersi_bt108(l1b_path)

    # Load MYD35
    mersi_dt = parse_mersi_datetime(recal_path)
    if mersi_dt is None:
        return []

    myd35_data = load_best_myd35_for_mersi(
        mersi_lat=lat, mersi_lon=lon, mersi_dt=mersi_dt,
        search_dirs=myd35_dirs,
        time_window_min=time_window_min, min_overlap=min_overlap,
    )
    if myd35_data is None:
        return []

    myd_resamp = myd35_data["clm_resampled"]

    # Classify pixels
    strata = classify_pixels(lat, lsm, sza, bt11)

    # Compute stats per stratum
    orbit_tag = re.search(r'(\d{8}_\d{4})', os.path.basename(recal_path))
    orbit_tag = orbit_tag.group(1) if orbit_tag else "unknown"

    rows = []
    for key, mask in strata.items():
        ov = mask & (recal_clm >= 0) & (myd_resamp >= 0)
        if ov.sum() < 100:
            continue

        s_recal = compute_validation_stats(
            np.where(ov, recal_clm, -1),
            np.where(ov, myd_resamp, -1),
            label="recal")
        s_onboard = compute_validation_stats(
            np.where(ov, onboard_clm, -1),
            np.where(ov, myd_resamp, -1),
            label="onboard")

        if s_recal and s_onboard:
            n_mersi = int((recal_clm[ov] >= 0).sum())
            m_cloud = 100.0 * (recal_clm[ov] <= 1).sum() / ov.sum()
            y_cloud = 100.0 * (myd_resamp[ov] <= 1).sum() / ov.sum()

            rows.append({
                "orbit": orbit_tag,
                "scene": key,
                "n_pixels": int(ov.sum()),
                "n_mersi": n_mersi,
                "mersi_cloud_pct": m_cloud,
                "myd35_cloud_pct": y_cloud,
                "cloud_bias": m_cloud - y_cloud,
                "recal_agree": s_recal["agree_pct"],
                "recal_pod": s_recal["pod"],
                "recal_far": s_recal["far"],
                "recal_hss": s_recal["hss"],
                "onboard_agree": s_onboard["agree_pct"],
                "onboard_pod": s_onboard["pod"],
                "onboard_far": s_onboard["far"],
                "onboard_hss": s_onboard["hss"],
            })

    return sorted(rows, key=lambda r: r["recal_hss"])


# ─────────────────────────────────────────────────────────────────────────────
# Reporter
# ─────────────────────────────────────────────────────────────────────────────

def print_stratified_table(rows: list[dict], title: str = "") -> None:
    """Print a formatted stratified validation table."""
    if not rows:
        print("No strata with sufficient samples.")
        return

    if title:
        print(f"\n{'='*100}")
        print(f"  {title}")
        print(f"{'='*100}")

    header = (f"{'Orbit':>10s}  {'Scene':<38s}  {'N':>7s}  "
              f"{'Agree':>7s}  {'POD':>7s}  {'FAR':>7s}  {'HSS':>7s}  "
              f"{'M_Cloud':>8s}  {'Y_Cloud':>8s}  {'Bias':>7s}")
    print(header)
    print("-" * len(header))

    for r in rows:
        print(f"{r['orbit']:>10s}  {r['scene']:<38s}  {r['n_pixels']:>7,}  "
              f"{r['recal_agree']:6.2f}%  {r['recal_pod']:6.2f}%  "
              f"{r['recal_far']:6.2f}%  {r['recal_hss']:7.4f}  "
              f"{r['mersi_cloud_pct']:7.1f}%  {r['myd35_cloud_pct']:7.1f}%  "
              f"{r['cloud_bias']:+7.1f}%")

    # Summary
    total_n = sum(r["n_pixels"] for r in rows)
    w_agree = sum(r["recal_agree"] * r["n_pixels"] for r in rows) / total_n
    w_pod = sum(r["recal_pod"] * r["n_pixels"] for r in rows) / total_n
    w_far = sum(r["recal_far"] * r["n_pixels"] for r in rows) / total_n
    w_hss = sum(r["recal_hss"] * r["n_pixels"] for r in rows) / total_n
    print("-" * len(header))
    print(f"{'WEIGHTED AVG':>10s}  {'':<38s}  {total_n:>7,}  "
          f"{w_agree:6.2f}%  {w_pod:6.2f}%  {w_far:6.2f}%  {w_hss:7.4f}")

    # Bottom-N by HSS
    print(f"\n--- Bottom 10 strata by HSS (most need improvement) ---")
    for r in rows[:10]:
        print(f"  HSS={r['recal_hss']:.4f}  {r['scene']}  "
              f"N={r['n_pixels']:,}  bias={r['cloud_bias']:+.1f}%  "
              f"({r['orbit']})")


def print_summary_by_dimension(rows: list[dict]) -> None:
    """Aggregate rows by surface type, day/night, latitude, BT range."""
    for dim in ["sfc", "dn", "lat", "bt"]:
        print(f"\n{'='*80}")
        print(f"  Aggregated by {dim}")
        print(f"{'='*80}")

        groups: dict[str, list] = {}
        for r in rows:
            parts = r["scene"].split("_")
            if dim == "sfc":
                key = parts[0]
            elif dim == "dn":
                key = parts[1]
            elif dim == "lat":
                key = parts[2]
            elif dim == "bt":
                key = "_".join(parts[3:]) if len(parts) > 3 else parts[3]
            else:
                key = "all"

            if key not in groups:
                groups[key] = []
            groups[key].append(r)

        for key in sorted(groups.keys()):
            grp = groups[key]
            tn = sum(r["n_pixels"] for r in grp)
            w = lambda field: sum(r[field] * r["n_pixels"] for r in grp) / tn
            print(f"  {key:<20s}  N={tn:>10,}  "
                  f"Agree={w('recal_agree'):5.1f}%  "
                  f"POD={w('recal_pod'):5.1f}%  "
                  f"FAR={w('recal_far'):5.1f}%  "
                  f"HSS={w('recal_hss'):.4f}  "
                  f"cloud_bias={w('cloud_bias'):+.1f}%")


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Stratified MERSI vs MYD35 cloud mask validation")
    parser.add_argument("--data_dir", required=True,
                        help="Directory with CLM pairs (e.g. .../20220803/)")
    parser.add_argument("--myd35_dir", nargs="+", default=["/data/Data_yuq/aqua_modis/MYD35_L2/"],
                        help="MYD35 search directories")
    parser.add_argument("--mersi_root", default="/data/Data_yuq/mersi")
    parser.add_argument("--time_window", type=int, default=15)
    parser.add_argument("--min_overlap", type=float, default=0.05)
    args = parser.parse_args()

    data_dir = Path(args.data_dir)
    recal_files = sorted(data_dir.glob("*_RECALI.HDF"))

    if not recal_files:
        # Try legacy naming
        recal_files = sorted(data_dir.glob("*_CLM_CLA_recal.h5"))

    all_rows = []
    for recal_path in recal_files:
        fname = recal_path.name
        # Derive onboard path
        onboard_path = recal_path.with_name(
            fname.replace("_RECALI.HDF", "_BUSINESS.HDF")
                .replace("_CLM_CLA_recal.h5", "_CLM_CLA.h5"))
        if not onboard_path.exists():
            print(f"[SKIP] No matching onboard for {fname}")
            continue

        print(f"\n[VALIDATE] {fname}")
        rows = validate_orbit(
            str(recal_path), str(onboard_path),
            myd35_dirs=args.myd35_dir,
            mersi_root=args.mersi_root,
            time_window_min=args.time_window,
            min_overlap=args.min_overlap,
        )
        all_rows.extend(rows)

    if not all_rows:
        print("[ERROR] No validation results.")
        exit(1)

    # Sort by HSS (worst first) for the main table
    all_rows.sort(key=lambda r: r["recal_hss"])

    print_stratified_table(all_rows, "Stratified Validation (RECAL vs MYD35)")
    print_summary_by_dimension(all_rows)
