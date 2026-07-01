"""
analyze_accuracy.py — 全天精度评估：MERSI CLM vs MYD35 真值，分层统计
=====================================================================
分层维度：
  - 总体 (overall)
  - 纬度带：高纬 60-90°, 中纬 30-60°, 低纬 0-30°
  - 下垫面：水体 (LandSeaMask=0), 陆地 (LandSeaMask=1)
  - 组合：纬度带 × 下垫面

只统计 MERSI 与 MYD35 重叠区域（pixel-level overlap > 50%）。
"""

from __future__ import annotations
import argparse
import sys, os, re, time, json
from pathlib import Path
from datetime import datetime, timedelta, timezone
from collections import defaultdict

import numpy as np

# Add visualize/ to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "visualize"))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "python"))

from io_mersi import load_clm_hdf5, parse_mersi_datetime
from io_myd35 import load_best_myd35_for_mersi, read_myd35, resample_to_mersi_grid
from fylat.mersi_io import read_l1b


# ─────────────────────────────────────────────────────────────────────────────
# Constants
# ─────────────────────────────────────────────────────────────────────────────

CLM_LABEL = {0: "Cloudy", 1: "Prob.Cloudy", 2: "Prob.Clear", 3: "Conf.Clear"}
CLM_DIR   = Path("/data/Data_yuq/fy3_cloud/20220803")
MYD35_DIR = "/data/Data_yuq/aqua_modis/MYD35_L2/20220803"
GEO_ROOT  = "/data/Data_yuq/mersi/20220803"
NWP_ROOT  = "/data/nwp"

LAT_BANDS = {
    "Low   0-30":   (0, 30),
    "Mid  30-60":   (30, 60),
    "High 60-90":   (60, 90),
}

SURFACE_MASKS = {
    "Water": lambda lsm: lsm == 0,
    "Land": lambda lsm: lsm == 1,
    "Coast/InlandWater": lambda lsm: (lsm >= 2) & (lsm <= 5),
    "OtherSurface": lambda lsm: lsm > 5,
}

FOCUS_STRATA = [
    "Overall",
    "Day SZA<80",
    "Mid  30-60  Land",
    "High 60-90  Land",
    "Coast/InlandWater",
    "Day MidLand",
    "Warm BT11>=270 Land",
]


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def geo_path_for_clm(clm_path: str) -> Path | None:
    """Return the matching 1-km GEO file path for a CLM product."""
    import h5py
    m = re.search(r'(\d{8})_(\d{4})', os.path.basename(clm_path))
    if not m:
        return None
    date_str, time_tag = m.group(1), m.group(2)
    root = Path(GEO_ROOT)
    if root.name != date_str:
        root = root / date_str
    geo_path = root / f"FY3D_MERSI_GBAL_L1_{date_str}_{time_tag}_GEO1K_MS.HDF"
    return geo_path if geo_path.exists() else None


def l1b_path_for_clm(clm_path: str) -> Path | None:
    """Return the matching 1-km L1B file path for a CLM product."""
    m = re.search(r'(\d{8})_(\d{4})', os.path.basename(clm_path))
    if not m:
        return None
    date_str, time_tag = m.group(1), m.group(2)
    root = Path(GEO_ROOT)
    if root.name != date_str:
        root = root / date_str
    l1b_path = root / f"FY3D_MERSI_GBAL_L1_{date_str}_{time_tag}_1000M_MS.HDF"
    return l1b_path if l1b_path.exists() else None


def _scaled_geo_dataset(h5, dataset_path: str) -> np.ndarray:
    ds = h5[dataset_path]
    raw = ds[:].astype(np.float64)
    slope = float(np.squeeze(ds.attrs.get("Slope", 1.0)))
    intercept = float(np.squeeze(ds.attrs.get("Intercept", 0.0)))
    fill_value = float(np.squeeze(ds.attrs.get("FillValue", -32767)))
    return np.where(raw != fill_value, (raw + intercept) * slope, np.nan)


def load_geo_context(clm_path: str) -> dict:
    """Load GEO context used for diagnostic stratification."""
    import h5py
    geo_path = geo_path_for_clm(clm_path)
    if geo_path is None:
        return {}
    if not geo_path.exists():
        return {}
    try:
        with h5py.File(str(geo_path), "r") as f:
            return {
                "lsm": f["Geolocation/LandSeaMask"][:].astype(np.int32),
                "sza": _scaled_geo_dataset(f, "Geolocation/SolarZenith").astype(np.float32),
            }
    except Exception:
        return {}


def load_l1b_features(clm_path: str) -> dict:
    """Load spectral features used to diagnose daytime threshold failures."""
    l1b_path = l1b_path_for_clm(clm_path)
    if l1b_path is None:
        return {}
    try:
        l1b = read_l1b(str(l1b_path))
        bt = l1b["bt_ir"]
        ref = l1b["ref_vis"]
        return {
            "BT11": bt[:, :, 4],
            "BTD11_12": bt[:, :, 4] - bt[:, :, 5],
            "BTD8_11": bt[:, :, 3] - bt[:, :, 4],
            "Ref065": ref[:, :, 2],
            "Ref086": ref[:, :, 3],
            "Ref138": ref[:, :, 18],
        }
    except Exception as exc:
        print(f"  [WARN] L1B spectral diagnostics unavailable: {exc}")
        return {}


def nearest_nwp_age_hours(mersi_dt: datetime | None) -> float | None:
    """Estimate the nearest available NWP valid-hour age for this orbit."""
    if mersi_dt is None:
        return None
    date_str = mersi_dt.strftime("%Y%m%d")
    org_dir = Path(NWP_ROOT) / date_str / "ORG"
    if not org_dir.exists():
        return None
    ages = []
    for path in org_dir.glob("gfs.t*z.pgrb2.0p25.f*"):
        m = re.search(r"gfs\.t(\d{2})z\.pgrb2\.0p25\.f(\d{3})", path.name)
        if not m:
            continue
        cycle_hour = int(m.group(1))
        lead_hour = int(m.group(2))
        cycle = datetime(
            mersi_dt.year, mersi_dt.month, mersi_dt.day,
            cycle_hour, tzinfo=timezone.utc,
        )
        valid_time = cycle + timedelta(hours=lead_hour)
        ages.append(abs((mersi_dt - valid_time).total_seconds()) / 3600.0)
    return min(ages) if ages else None


def _crop_like(arr: np.ndarray | None, shape: tuple[int, int]) -> np.ndarray | None:
    if arr is None:
        return None
    if arr.shape == shape:
        return arr
    return arr[:shape[0], :shape[1]]


def _nwp_age_label(age_hours: float | None) -> str | None:
    if age_hours is None:
        return None
    if age_hours <= 1.5:
        return "NWP age <=1.5h"
    if age_hours <= 3.0:
        return "NWP age 1.5-3h"
    if age_hours <= 6.0:
        return "NWP age 3-6h"
    return "NWP age >6h"


# ─────────────────────────────────────────────────────────────────────────────
# Stratified statistics
# ─────────────────────────────────────────────────────────────────────────────

def compute_stratified_stats(
    mersi_clm: np.ndarray,
    myd35_clm: np.ndarray,
    lat: np.ndarray,
    lsm: np.ndarray | None,
    sza: np.ndarray | None = None,
    bt11: np.ndarray | None = None,
    tag: str | None = None,
    nwp_age_hours: float | None = None,
) -> dict:
    """
    Compute validation metrics stratified by latitude band and surface type.
    Only pixels with valid MERSI CLM, valid MYD35 CLM, and valid lat/lsm are used.
    """
    # Base valid mask: both CLMs valid + lat valid
    base_mask = (mersi_clm >= 0) & (myd35_clm >= 0) & np.isfinite(lat)
    if not base_mask.any():
        return {}

    abs_lat = np.abs(lat)

    stratifications = {}

    def add(name: str, mask: np.ndarray) -> None:
        stratifications[name] = base_mask & mask

    # ── Overall ──
    add("Overall", np.ones_like(base_mask, dtype=bool))

    # ── Latitude bands ──
    for name, (lo, hi) in LAT_BANDS.items():
        add(name, (abs_lat >= lo) & (abs_lat < hi))

    # ── Surface type (if LSM available) ──
    if lsm is not None:
        surface_masks = {name: fn(lsm) for name, fn in SURFACE_MASKS.items()}
        for name, mask in surface_masks.items():
            add(name, mask)

        # ── Combinations ──
        for lat_name, (lo, hi) in LAT_BANDS.items():
            lat_mask = (abs_lat >= lo) & (abs_lat < hi)
            for surface_name, surface_mask in surface_masks.items():
                add(f"{lat_name}  {surface_name}", lat_mask & surface_mask)

    # ── Solar zenith / day-night diagnostics ──
    if sza is not None:
        add("Day SZA<80", sza < 80.0)
        add("Twilight SZA80-90", (sza >= 80.0) & (sza < 90.0))
        add("Night SZA>=90", sza >= 90.0)
        add("Day SZA0-60", sza < 60.0)
        add("Day SZA60-80", (sza >= 60.0) & (sza < 80.0))
        if lsm is not None:
            land_mask = lsm == 1
            coast_mask = (lsm >= 2) & (lsm <= 5)
            add("Day MidLand", (sza < 80.0) & (abs_lat >= 30.0) & (abs_lat < 60.0) & land_mask)
            add("Day Coast/InlandWater", (sza < 80.0) & coast_mask)

    # ── Hemisphere ──
    add("North", lat >= 0)
    add("South", lat < 0)

    # ── Thermal background diagnostics ──
    if bt11 is not None:
        cold = bt11 < 230.0
        moderate = (bt11 >= 230.0) & (bt11 < 270.0)
        warm = bt11 >= 270.0
        add("Cold BT11<230", cold)
        add("Moderate BT11 230-270", moderate)
        add("Warm BT11>=270", warm)
        if lsm is not None:
            add("Warm BT11>=270 Land", warm & (lsm == 1))
            add("Cold BT11<230 HighLand", cold & (abs_lat >= 60.0) & (lsm == 1))

    if tag:
        hhmm = tag.split("_")[-1]
        hour = int(hhmm[:2])
        add(f"Time {hhmm}", np.ones_like(base_mask, dtype=bool))
        add(f"UTC {hour // 6 * 6:02d}-{hour // 6 * 6 + 6:02d}", np.ones_like(base_mask, dtype=bool))

    nwp_label = _nwp_age_label(nwp_age_hours)
    if nwp_label:
        add(nwp_label, np.ones_like(base_mask, dtype=bool))

    # ── Compute metrics for each stratum ──
    results = {}
    for name, mask in stratifications.items():
        if mask.sum() < 100:
            continue  # skip tiny strata

        a = mersi_clm[mask]
        b = myd35_clm[mask]
        n = len(a)

        agree_pct = 100.0 * np.mean(a == b)

        a_cloud = a <= 1
        b_cloud = b <= 1
        TP = int((a_cloud & b_cloud).sum())
        FP = int((a_cloud & ~b_cloud).sum())
        FN = int((~a_cloud & b_cloud).sum())
        TN = int((~a_cloud & ~b_cloud).sum())

        pod = 100.0 * TP / (TP + FN + 1e-9)
        far = 100.0 * FP / (TP + FP + 1e-9)
        csi = 100.0 * TP / (TP + FP + FN + 1e-9)
        expected = ((TP + FP) * (TP + FN) + (TN + FP) * (TN + FN)) / (n + 1e-9)
        hss = (TP + TN - expected) / (n - expected + 1e-9)

        # Per-class agreement
        class_agree = {}
        for c in range(4):
            cm = (a == c)
            if cm.sum() > 0:
                class_agree[f"agr_c{c}"] = 100.0 * np.mean(b[cm] == c)
            else:
                class_agree[f"agr_c{c}"] = np.nan

        # Class distribution (MERSI / MYD35)
        mer_dist = {f"mer_c{c}": 100.0 * (a == c).sum() / n for c in range(4)}
        myd_dist = {f"myd_c{c}": 100.0 * (b == c).sum() / n for c in range(4)}

        results[name] = dict(
            n=int(n), agree=agree_pct, pod=pod, far=far,
            csi=csi, hss=hss,
            TP=TP, FP=FP, FN=FN, TN=TN,
            **class_agree, **mer_dist, **myd_dist,
        )

    return results


def compute_feature_stats(
    features: dict,
    masks: dict[str, np.ndarray],
    myd35_clm: np.ndarray,
) -> dict:
    """Summarize feature means/stds by MYD35 cloud/clear label."""
    if not features:
        return {}

    out = {}
    truth_cloud = myd35_clm <= 1
    truth_clear = myd35_clm >= 2
    for mask_name, mask in masks.items():
        group = {}
        for feature_name, arr in features.items():
            if arr is None:
                continue
            arr = _crop_like(arr, myd35_clm.shape)
            if arr is None:
                continue
            valid = mask & np.isfinite(arr)
            stats = {}
            for label, label_mask in [("cloud", truth_cloud), ("clear", truth_clear)]:
                values = arr[valid & label_mask].astype(np.float64)
                if values.size < 100:
                    continue
                stats[label] = {
                    "n": int(values.size),
                    "sum": float(values.sum()),
                    "sumsq": float(np.square(values).sum()),
                }
            if stats:
                group[feature_name] = stats
        if group:
            out[mask_name] = group
    return out


def accumulate_feature_stats(results_list: list[dict]) -> dict:
    acc = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: {"n": 0, "sum": 0.0, "sumsq": 0.0})))
    for result in results_list:
        for mask_name, features in result.items():
            for feature_name, labels in features.items():
                for label, stats in labels.items():
                    slot = acc[mask_name][feature_name][label]
                    slot["n"] += int(stats["n"])
                    slot["sum"] += float(stats["sum"])
                    slot["sumsq"] += float(stats["sumsq"])

    final = {}
    for mask_name, features in acc.items():
        final[mask_name] = {}
        for feature_name, labels in features.items():
            final[mask_name][feature_name] = {}
            for label, stats in labels.items():
                n = stats["n"]
                if n == 0:
                    continue
                mean = stats["sum"] / n
                var = max(stats["sumsq"] / n - mean * mean, 0.0)
                final[mask_name][feature_name][label] = {
                    "n": n,
                    "mean": mean,
                    "std": var ** 0.5,
                }
    return final


# ─────────────────────────────────────────────────────────────────────────────
# Accumulator
# ─────────────────────────────────────────────────────────────────────────────

def accumulate(results_list: list[dict]) -> dict:
    """Merge per-orbit stats into weighted aggregate."""
    acc = defaultdict(lambda: {"n": 0, "n_agree": 0})
    # For cloud detection metrics we need the full confusion matrix
    # So we accumulate TP/FP/FN/TN, and per-class agreement separately
    acc_cm = defaultdict(lambda: {"TP": 0, "FP": 0, "FN": 0, "TN": 0})
    acc_class = defaultdict(lambda: {
        "mer_c0": 0, "mer_c1": 0, "mer_c2": 0, "mer_c3": 0,
        "myd_c0": 0, "myd_c1": 0, "myd_c2": 0, "myd_c3": 0,
        "class_n": 0,          # pixel count for class dist
        "agr_c0_num": 0,       # numerator for per-class agreement
        "agr_c1_num": 0,
        "agr_c2_num": 0,
        "agr_c3_num": 0,
        "agr_c0_den": 0,       # denominator for per-class agreement
        "agr_c1_den": 0,
        "agr_c2_den": 0,
        "agr_c3_den": 0,
    })

    for r in results_list:
        for name, stats in r.items():
            n = stats["n"]
            acc[name]["n"] += n
            acc[name]["n_agree"] += int(round(stats["agree"] * n / 100.0))

            acc_cm[name]["TP"] += stats.get("TP", 0)
            acc_cm[name]["FP"] += stats.get("FP", 0)
            acc_cm[name]["FN"] += stats.get("FN", 0)
            acc_cm[name]["TN"] += stats.get("TN", 0)

            # Class distribution weighted by n
            cm_key = acc_class[name]
            cm_key["class_n"] += n
            for c in range(4):
                cm_key[f"mer_c{c}"] += round(stats.get(f"mer_c{c}", 0) * n / 100.0)
                cm_key[f"myd_c{c}"] += round(stats.get(f"myd_c{c}", 0) * n / 100.0)
                # Per-class agreement - we need to reconstruct numerator/denominator
                # agr_cX = % of MERSI class X pixels where MYD35 agrees
                # = (number of pixels where a==X and b==X) / (number of pixels where a==X)
                # So den = mer_cX count, num = den * agr_cX / 100
                denom = round(stats.get(f"mer_c{c}", 0) * n / 100.0)
                cm_key[f"agr_c{c}_den"] += denom
                agr_val = stats.get(f"agr_c{c}", np.nan)
                if not np.isnan(agr_val):
                    cm_key[f"agr_c{c}_num"] += round(denom * agr_val / 100.0)

    # Compute final metrics
    final = {}
    for name in acc:
        n = acc[name]["n"]
        n_agree = acc[name]["n_agree"]
        cm = acc_cm[name]
        TP, FP, FN, TN = cm["TP"], cm["FP"], cm["FN"], cm["TN"]

        pod = 100.0 * TP / (TP + FN + 1e-9)
        far = 100.0 * FP / (TP + FP + 1e-9)
        csi = 100.0 * TP / (TP + FP + FN + 1e-9)
        n_total = TP + FP + FN + TN
        expected = ((TP + FP) * (TP + FN) + (TN + FP) * (TN + FN)) / (n_total + 1e-9)
        hss = (TP + TN - expected) / (n_total - expected + 1e-9) if n_total > 0 else 0

        cl = acc_class[name]
        # Class distribution
        mer_dist = {c: 100.0 * cl[f"mer_c{c}"] / cl["class_n"] for c in range(4) if cl["class_n"] > 0}
        myd_dist = {c: 100.0 * cl[f"myd_c{c}"] / cl["class_n"] for c in range(4) if cl["class_n"] > 0}
        # Per-class agreement
        class_agree = {}
        for c in range(4):
            den = cl[f"agr_c{c}_den"]
            if den > 0:
                class_agree[c] = 100.0 * cl[f"agr_c{c}_num"] / den
            else:
                class_agree[c] = np.nan

        final[name] = dict(
            n=n, agree=100.0 * n_agree / n,
            pod=pod, far=far, csi=csi, hss=hss,
            mer_dist=mer_dist, myd_dist=myd_dist,
            class_agree=class_agree,
        )

    return final


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def find_orbit_pairs():
    """Discover recal/business CLM pairs in 20220803 directory."""
    pairs = []
    for business in sorted(CLM_DIR.glob("*_BUSINESS.HDF")):
        recal = business.with_name(business.name.replace("_BUSINESS.HDF", "_RECALI.HDF"))
        if not recal.exists():
            continue
        m = re.search(r'(\d{8})_(\d{4})', business.name)
        if not m:
            continue
        pairs.append((str(recal), str(business), m.group(1), m.group(2)))
    return pairs


def process_orbit(
    recal_path,
    business_path,
    date_str,
    time_str,
    spectral_diagnostics: bool,
    time_window_min: int,
    min_overlap: float,
):
    """Process a single orbit, return stratified stats dict or None."""
    tag = f"{date_str}_{time_str}"
    print(f"\n{'='*60}")
    print(f"  [{tag}]")

    # Load CLM
    recal_data   = load_clm_hdf5(recal_path)
    onboard_data = load_clm_hdf5(business_path)
    if recal_data is None or onboard_data is None:
        print(f"  [SKIP] CLM loading failed")
        return None

    lat = recal_data["lat"]
    lon = recal_data["lon"]
    recal_clm   = recal_data["clm"]
    onboard_clm = onboard_data["clm"]

    # Load GEO/L1B diagnostic context.
    geo_context = load_geo_context(business_path)
    lsm = _crop_like(geo_context.get("lsm"), lat.shape)
    sza = _crop_like(geo_context.get("sza"), lat.shape)
    features = load_l1b_features(business_path) if spectral_diagnostics else {}
    bt11 = _crop_like(features.get("BT11"), lat.shape)
    print(f"  LSM: {'available' if lsm is not None else 'missing'}")
    print(f"  SZA: {'available' if sza is not None else 'missing'}")
    print(f"  L1B spectral features: {'available' if features else 'missing/skipped'}")

    # Match MYD35
    mersi_dt = parse_mersi_datetime(business_path)
    if mersi_dt is None:
        print(f"  [SKIP] Could not parse datetime")
        return None

    myd35_data = load_best_myd35_for_mersi(
        mersi_lat=lat,
        mersi_lon=lon,
        mersi_dt=mersi_dt,
        search_dirs=MYD35_DIR,
        mersi_clm=recal_clm,
        time_window_min=time_window_min,
        min_overlap=min_overlap,
    )
    if myd35_data is None:
        print(f"  [SKIP] No matching MYD35")
        return None

    myd35_resampled = myd35_data["clm_resampled"]
    nwp_age_hours = nearest_nwp_age_hours(mersi_dt)

    # Compute stats for both calibrations
    print(f"  Computing stratified stats...")
    stats_recal = compute_stratified_stats(
        recal_clm, myd35_resampled, lat, lsm,
        sza=sza, bt11=bt11, tag=tag, nwp_age_hours=nwp_age_hours,
    )
    stats_onboard = compute_stratified_stats(
        onboard_clm, myd35_resampled, lat, lsm,
        sza=sza, bt11=bt11, tag=tag, nwp_age_hours=nwp_age_hours,
    )

    focus_masks = {}
    base = (myd35_resampled >= 0) & np.isfinite(lat)
    if lsm is not None:
        abs_lat = np.abs(lat)
        focus_masks["MidLand"] = base & (abs_lat >= 30.0) & (abs_lat < 60.0) & (lsm == 1)
        focus_masks["HighLand"] = base & (abs_lat >= 60.0) & (lsm == 1)
        focus_masks["CoastInlandWater"] = base & (lsm >= 2) & (lsm <= 5)
        if sza is not None:
            focus_masks["DayMidLand"] = focus_masks["MidLand"] & (sza < 80.0)
        if bt11 is not None:
            focus_masks["WarmLand"] = focus_masks["MidLand"] & (bt11 >= 270.0)
    focus_masks["Overall"] = base
    feature_stats = compute_feature_stats(features, focus_masks, myd35_resampled)

    return {
        "tag": tag,
        "recal": stats_recal,
        "onboard": stats_onboard,
        "feature_stats": feature_stats,
        "nwp_age_hours": nwp_age_hours,
        "dt_diff_min": myd35_data.get("dt_diff_min", 0),
        "pixel_overlap": myd35_data.get("pixel_overlap", 0),
    }


def print_summary(acc_recal, acc_onboard, n_orbits, total_pixels, feature_stats=None):
    """Pretty-print the final stratified accuracy table."""

    def print_section(title, acc, n_orbits):
        print(f"\n{'═'*110}")
        print(f"  {title}")
        print(f"  Orbits processed: {n_orbits}")
        print(f"{'─'*110}")
        header = (f"  {'Stratum':<28s} {'Pixels':>12s} {'Agree':>7s} "
                  f"{'POD':>7s} {'FAR':>7s} {'CSI':>7s} {'HSS':>7s} "
                  f"{'CldAg':>7s} {'PClAg':>7s} {'PClrAg':>7s} {'ClrAg':>7s}")
        print(header)
        print(f"  {'':-<110s}")

        # Sort order
        order = [
            "Overall",
            "Low   0-30",
            "Mid  30-60",
            "High 60-90",
            "Water",
            "Land",
            "Low   0-30  Water",
            "Low   0-30  Land",
            "Mid  30-60  Water",
            "Mid  30-60  Land",
            "High 60-90  Water",
            "High 60-90  Land",
            "Coast/InlandWater",
            "Low   0-30  Coast/InlandWater",
            "Mid  30-60  Coast/InlandWater",
            "High 60-90  Coast/InlandWater",
            "Day SZA<80",
            "Twilight SZA80-90",
            "Night SZA>=90",
            "Day SZA0-60",
            "Day SZA60-80",
            "Day MidLand",
            "Day Coast/InlandWater",
            "Cold BT11<230",
            "Moderate BT11 230-270",
            "Warm BT11>=270",
            "Warm BT11>=270 Land",
            "Cold BT11<230 HighLand",
            "NWP age <=1.5h",
            "NWP age 1.5-3h",
            "NWP age 3-6h",
            "NWP age >6h",
            "North",
            "South",
        ]
        for name in order:
            if name not in acc:
                continue
            s = acc[name]
            ca = s.get("class_agree", {})
            print(f"  {name:<28s} {s['n']:>12,d} {s['agree']:>6.2f}% "
                  f"{s['pod']:>6.2f}% {s['far']:>6.2f}% {s['csi']:>6.2f}% "
                  f"{s['hss']:>6.3f}  "
                  f"{ca.get(0, np.nan):>6.2f}% {ca.get(1, np.nan):>6.2f}% "
                  f"{ca.get(2, np.nan):>6.2f}% {ca.get(3, np.nan):>6.2f}%")

        # Print class distribution
        print(f"\n  Class distribution (%):")
        dist_header = (f"  {'Stratum':<28s} "
                       f"{'M:Cld':>7s} {'M:PCld':>7s} {'M:PClr':>7s} {'M:Clr':>7s}"
                       f"  |  "
                       f"{'T:Cld':>7s} {'T:PCld':>7s} {'T:PClr':>7s} {'T:Clr':>7s}")
        print(dist_header)
        print(f"  {'':-<95s}")
        for name in order:
            if name not in acc:
                continue
            s = acc[name]
            md = s.get("mer_dist", {})
            td = s.get("myd_dist", {})
            print(f"  {name:<28s} "
                  f"{md.get(0,0):>6.1f}% {md.get(1,0):>6.1f}% "
                  f"{md.get(2,0):>6.1f}% {md.get(3,0):>6.1f}%  |  "
                  f"{td.get(0,0):>6.1f}% {td.get(1,0):>6.1f}% "
                  f"{td.get(2,0):>6.1f}% {td.get(3,0):>6.1f}%")

        rows = [
            (s["hss"], s["agree"], s["pod"], s["far"], s["n"], name)
            for name, s in acc.items()
            if name != "Overall" and s.get("n", 0) >= 1000
        ]
        print(f"\n  Worst HSS strata:")
        for hss, agree, pod, far, n, name in sorted(rows)[:10]:
            print(f"  {name:<32s} n={n:>10,d} HSS={hss:>6.3f} "
                  f"Agree={agree:>6.2f}% POD={pod:>6.2f}% FAR={far:>6.2f}%")

        print(f"\n  Highest FAR strata:")
        for hss, agree, pod, far, n, name in sorted(rows, key=lambda x: x[3], reverse=True)[:10]:
            print(f"  {name:<32s} n={n:>10,d} HSS={hss:>6.3f} "
                  f"Agree={agree:>6.2f}% POD={pod:>6.2f}% FAR={far:>6.2f}%")

        print(f"\n  Lowest POD strata:")
        for hss, agree, pod, far, n, name in sorted(rows, key=lambda x: x[2])[:10]:
            print(f"  {name:<32s} n={n:>10,d} HSS={hss:>6.3f} "
                  f"Agree={agree:>6.2f}% POD={pod:>6.2f}% FAR={far:>6.2f}%")

    # Print recal
    print_section("RECALIBRATION  (recal)", acc_recal, n_orbits)
    # Print onboard
    print_section("ONBOARD / BUSINESS  (onboard)", acc_onboard, n_orbits)

    # Summary comparison
    print(f"\n{'═'*110}")
    print(f"  RECAL vs ONBOARD — key metrics comparison (Overall)")
    print(f"{'─'*110}")
    if "Overall" in acc_recal and "Overall" in acc_onboard:
        r = acc_recal["Overall"]
        o = acc_onboard["Overall"]
        for key in ["agree", "pod", "far", "csi", "hss"]:
            diff = r[key] - o[key]
            print(f"  Δ {key.upper():<6s}: recal {r[key]:.4f}  onboard {o[key]:.4f}  diff {diff:+.4f}")

    if feature_stats:
        print(f"\n{'═'*110}")
        print("  Spectral diagnostics by MYD35 cloud/clear label")
        print(f"{'─'*110}")
        for mask_name in ["Overall", "DayMidLand", "MidLand", "HighLand", "CoastInlandWater", "WarmLand"]:
            if mask_name not in feature_stats:
                continue
            print(f"\n  [{mask_name}]")
            for feature_name, labels in feature_stats[mask_name].items():
                cloud = labels.get("cloud")
                clear = labels.get("clear")
                if not cloud or not clear:
                    continue
                delta = cloud["mean"] - clear["mean"]
                print(
                    f"    {feature_name:<10s} "
                    f"cloud={cloud['mean']:>8.3f}±{cloud['std']:<7.3f} "
                    f"clear={clear['mean']:>8.3f}±{clear['std']:<7.3f} "
                    f"Δ={delta:>8.3f} "
                    f"n={cloud['n']:,}/{clear['n']:,}"
                )

    print(f"\n{'═'*110}")


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--date", default="20220803", help="Observation date YYYYMMDD")
    parser.add_argument("--clm-dir", default=None,
                        help="Directory containing BUSINESS/RECALI CLM products")
    parser.add_argument("--myd35-dir", default=None,
                        help="Directory containing MYD35_L2 reference files")
    parser.add_argument("--geo-root", default=None,
                        help="MERSI L1/GEO root or date directory")
    parser.add_argument("--nwp-root", default=NWP_ROOT,
                        help="NWP root used to estimate nearest forecast age")
    parser.add_argument("--output", default=None,
                        help="Output JSON path")
    parser.add_argument("--time-window-min", type=int, default=15,
                        help="MYD35 temporal search window")
    parser.add_argument("--min-overlap", type=float, default=0.50,
                        help="Minimum MYD35/MERSI overlap")
    parser.add_argument("--no-spectral-diagnostics", action="store_true",
                        help="Skip L1B spectral feature diagnostics")
    return parser.parse_args()


def main():
    global CLM_DIR, MYD35_DIR, GEO_ROOT, NWP_ROOT
    args = parse_args()
    CLM_DIR = Path(args.clm_dir or f"/data/Data_yuq/fy3_cloud/{args.date}")
    MYD35_DIR = args.myd35_dir or f"/data/Data_yuq/aqua_modis/MYD35_L2/{args.date}"
    GEO_ROOT = args.geo_root or f"/data/Data_yuq/mersi/{args.date}"
    NWP_ROOT = args.nwp_root

    pairs = find_orbit_pairs()
    print(f"Found {len(pairs)} orbit pairs for {args.date}")

    recal_results   = []
    onboard_results = []
    feature_results = []
    n_matched = 0
    t0 = time.time()

    for recal_path, business_path, date_str, time_str in pairs:
        result = process_orbit(
            recal_path, business_path, date_str, time_str,
            spectral_diagnostics=not args.no_spectral_diagnostics,
            time_window_min=args.time_window_min,
            min_overlap=args.min_overlap,
        )
        if result is None:
            continue
        n_matched += 1
        if result.get("recal"):
            recal_results.append(result["recal"])
        if result.get("onboard"):
            onboard_results.append(result["onboard"])
        if result.get("feature_stats"):
            feature_results.append(result["feature_stats"])

        print(f"  [OK] matched orbits so far: {n_matched}")

    elapsed = time.time() - t0
    print(f"\n{'='*60}")
    print(f"Done. {n_matched}/{len(pairs)} orbits matched MYD35 in {elapsed:.0f}s")

    if n_matched == 0:
        print("No orbits matched — exiting.")
        return

    # Accumulate
    acc_recal   = accumulate(recal_results)
    acc_onboard = accumulate(onboard_results)
    acc_features = accumulate_feature_stats(feature_results)

    total_pixels = acc_recal.get("Overall", {}).get("n", 0)
    print_summary(acc_recal, acc_onboard, n_matched, total_pixels, acc_features)

    # Save JSON for later use
    output = {
        "n_orbits": n_matched,
        "total_pixels": total_pixels,
        "recal": acc_recal,
        "onboard": acc_onboard,
        "feature_stats": acc_features,
    }
    output_path = args.output or f"accuracy_{args.date}.json"
    with open(output_path, "w") as f:
        json.dump(output, f, indent=2, default=lambda x: float(x) if isinstance(x, (np.integer, np.floating)) else str(x))
    print(f"\nSaved results to {output_path}")


if __name__ == "__main__":
    main()
