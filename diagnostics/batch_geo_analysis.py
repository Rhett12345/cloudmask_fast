"""
batch_geo_analysis.py — Day/night + geographic region stratified MYD35 comparison.

Usage:
  python batch_geo_analysis.py --date 20220803
"""

import argparse, os, re, sys
from pathlib import Path
import numpy as np
import h5py

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "visualize"))
from io_mersi import load_clm_hdf5, parse_mersi_datetime
from io_myd35 import load_best_myd35_for_mersi

DATA_DIR = "/data/Data_yuq/fy3_cloud/"
MYD35_DIR = "/data/Data_yuq/aqua_modis/MYD35_L2/"
MERSI_ROOT = "/data/Data_yuq/mersi/"


def classify_sza(sza):
    """Classify by solar zenith angle."""
    day = sza < 80
    twilight = (sza >= 80) & (sza < 90)
    night = sza >= 90
    return day, twilight, night


def compute_stats(mersi_clm, myd_clm, mask):
    """Compute validation stats on a boolean mask."""
    valid = mask & (mersi_clm >= 0) & (myd_clm >= 0)
    n = valid.sum()
    if n < 50:
        return None
    m_cloud = (mersi_clm[valid] <= 1)
    y_cloud = (myd_clm[valid] <= 1)
    tp = (m_cloud & y_cloud).sum()
    fp = (m_cloud & ~y_cloud).sum()
    fn = (~m_cloud & y_cloud).sum()
    tn = (~m_cloud & ~y_cloud).sum()
    nv = tp + fp + fn + tn
    pod = tp / (tp + fn + 1e-9)
    far = fp / (tp + fp + 1e-9)
    agree = (tp + tn) / nv
    expected = ((tp + fp) * (tp + fn) + (tn + fp) * (tn + fn)) / (nv + 1e-9)
    hss = (tp + tn - expected) / (nv - expected + 1e-9)
    return {"n": int(nv), "agree": agree, "pod": pod, "far": far, "hss": hss,
            "m_cloud": m_cloud.mean(), "y_cloud": y_cloud.mean(),
            "tp": int(tp), "fp": int(fp), "fn": int(fn), "tn": int(tn)}


def load_geo_ancillary(l1b_path):
    """Load SZA and LandSeaMask from GEO file.

    GEO datasets are stored as scaled integers (dtype int16/uint16).
    Real value = raw * Slope + Intercept (from HDF5 attributes).
    """
    geo_path = l1b_path.replace("_1000M_MS.HDF", "_GEO1K_MS.HDF")
    try:
        with h5py.File(geo_path, "r") as f:
            lsm = f["Geolocation/LandSeaMask"][:].astype(np.int32)
            lat = f["Geolocation/Latitude"][:].astype(np.float64)

            sz_ds = f["Geolocation/SolarZenith"]
            sz_slope = sz_ds.attrs.get("Slope", 1.0)
            sz_intercept = sz_ds.attrs.get("Intercept", 0.0)
            sz_fill = sz_ds.attrs.get("FillValue", -32767)
            sza_raw = sz_ds[:].astype(np.float64)
            sza = np.where(sza_raw != sz_fill, sza_raw * sz_slope + sz_intercept, np.nan)
        return lat, lsm, sza
    except Exception:
        return None, None, None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--date", default="20220803")
    args = parser.parse_args()

    clm_dir = os.path.join(DATA_DIR, args.date)
    recal_files = sorted(Path(clm_dir).glob("*_RECALI.HDF"))

    print(f"{'Time':>6s} {'dn':>4s} {'lat_band':>9s} {'sfc':>5s} "
          f"{'N':>7s} {'Agree':>7s} {'POD':>6s} {'FAR':>6s} {'HSS':>7s} "
          f"{'M%':>6s} {'Y%':>6s} {'Bias':>6s}")
    print("-" * 95)

    all_rows = []
    summary = {"day": [], "twilight": [], "night": [],
               "tropical": [], "midlat": [], "polar": [],
               "ocean": [], "land": [], "inland": []}

    for rpath in recal_files:
        ts = re.search(r'_(\d{4})_', rpath.name).group(1)
        data = load_clm_hdf5(str(rpath))
        if data is None:
            continue
        clm = data["clm"]

        dt = parse_mersi_datetime(str(rpath))
        if dt is None:
            continue

        # Load GEO
        l1b = os.path.join(MERSI_ROOT, args.date,
                           f"FY3D_MERSI_GBAL_L1_{args.date}_{ts}_1000M_MS.HDF")
        lat, lsm, sza = load_geo_ancillary(l1b)
        if lat is None:
            continue

        # Load MYD35
        myd = load_best_myd35_for_mersi(
            mersi_lat=data["lat"], mersi_lon=data["lon"], mersi_dt=dt,
            search_dirs=[os.path.join(MYD35_DIR, args.date)],
            time_window_min=15, min_overlap=0.20,
        )
        if myd is None:
            continue
        myd_clm = myd["clm_resampled"]

        # Classify
        day, twilight, night = classify_sza(sza)

        abs_lat = np.abs(lat)
        tropical = abs_lat < 30
        midlat = (abs_lat >= 30) & (abs_lat < 60)
        polar = abs_lat >= 60

        # FY-3D LandSeaMask: 0=ShallowOcean,1=Land,2=Coast/LakeShore,
        #   3=ShallowInlandWater,4=EphemeralWater,5=DeepInlandWater,
        #   6=ModerateOcean,7=DeepOcean
        ocean = (lsm == 0) | (lsm == 6) | (lsm == 7)
        land = (lsm == 1)
        inland = (lsm >= 2) & (lsm <= 5)

        # Stratify: 3 dn x 3 lat x 3 sfc = 27 strata
        dn_map = [("Day", day), ("Twi", twilight), ("Nit", night)]
        lat_map = [("Tropic", tropical), ("MidLat", midlat), ("Polar", polar)]
        sfc_map = [("Ocn", ocean), ("Lnd", land), ("Inl", inland)]

        for dn_name, dn_mask in dn_map:
            for lat_name, lat_mask in lat_map:
                for sfc_name, sfc_mask in sfc_map:
                    mask = dn_mask & lat_mask & sfc_mask
                    if mask.sum() < 100:
                        continue
                    s = compute_stats(clm, myd_clm, mask)
                    if s is None:
                        continue
                    bias = (s["m_cloud"] - s["y_cloud"]) * 100
                    row = (ts, dn_name, lat_name, sfc_name, s)
                    all_rows.append(row)

                    summary["day"].append(s["hss"]) if dn_name == "Day" else None
                    summary["twilight"].append(s["hss"]) if dn_name == "Twi" else None
                    summary["night"].append(s["hss"]) if dn_name == "Nit" else None
                    summary["tropical"].append(s["hss"]) if lat_name == "Tropic" else None
                    summary["midlat"].append(s["hss"]) if lat_name == "MidLat" else None
                    summary["polar"].append(s["hss"]) if lat_name == "Polar" else None
                    summary["ocean"].append(s["hss"]) if sfc_name == "Ocn" else None
                    summary["land"].append(s["hss"]) if sfc_name == "Lnd" else None
                    summary["inland"].append(s["hss"]) if sfc_name == "Inl" else None

    # Print all rows sorted by HSS
    all_rows.sort(key=lambda r: r[4]["hss"])
    for ts, dn, latb, sfc, s in all_rows:
        bias = (s["m_cloud"] - s["y_cloud"]) * 100
        print(f"{ts:>6s} {dn:>4s} {latb:>9s} {sfc:>5s} {s['n']:>7,} "
              f"{s['agree']*100:6.1f}% {s['pod']*100:5.1f}% {s['far']*100:5.1f}% "
              f"{s['hss']:7.4f} {s['m_cloud']*100:5.1f}% {s['y_cloud']*100:5.1f}% "
              f"{bias:+5.1f}%")

    # Summary by dimension
    print(f"\n{'='*60}")
    print("  Summary by dimension (mean HSS)")
    print(f"{'='*60}")
    for label, vals in [("Day", summary["day"]), ("Twilight", summary["twilight"]),
                          ("Night", summary["night"]), ("-", []),
                          ("Tropical", summary["tropical"]), ("MidLat", summary["midlat"]),
                          ("Polar", summary["polar"]), ("-", []),
                          ("Ocean", summary["ocean"]), ("Land", summary["land"]),
                          ("InlandWater", summary["inland"])]:
        if label == "-":
            print()
            continue
        vv = [v for v in vals if not np.isnan(v) and np.isfinite(v)]
        if vv:
            print(f"  {label:<15s}  N={len(vv):>5d}  "
                  f"HSS_mean={np.mean(vv):.4f}  HSS_p50={np.median(vv):.4f}  "
                  f"HSS_worst={np.min(vv):.4f}")

    # Bottom 10 worst
    print(f"\n{'='*60}")
    print("  Worst 10 strata (most need improvement)")
    print(f"{'='*60}")
    for ts, dn, latb, sfc, s in all_rows[:10]:
        bias = (s["m_cloud"] - s["y_cloud"]) * 100
        print(f"  {ts} {dn} {latb} {sfc}: HSS={s['hss']:.4f}  "
              f"POD={s['pod']*100:.0f}%  FAR={s['far']*100:.0f}%  "
              f"N={s['n']:,}  bias={bias:+.1f}%")

    # Best 10
    print(f"\n  Best 10 strata")
    print(f"{'='*60}")
    for ts, dn, latb, sfc, s in all_rows[-10:]:
        bias = (s["m_cloud"] - s["y_cloud"]) * 100
        print(f"  {ts} {dn} {latb} {sfc}: HSS={s['hss']:.4f}  "
              f"POD={s['pod']*100:.0f}%  FAR={s['far']*100:.0f}%  "
              f"N={s['n']:,}  bias={bias:+.1f}%")

    return 0


if __name__ == "__main__":
    sys.exit(main())
