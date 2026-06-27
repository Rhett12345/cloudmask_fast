"""
batch_myd35_overlap.py — Match all MERSI time slots against MYD35 and rank by consistency.

Usage:
  python batch_myd35_overlap.py --date 20220803
"""

import argparse, os, re, sys
from pathlib import Path
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "visualize"))
from io_mersi import load_clm_hdf5, parse_mersi_datetime
from io_myd35 import load_best_myd35_for_mersi


def classify_day_night(lat, sza):
    """Rough day/night classification."""
    return "Day" if np.nanmean(sza) < 85 else "Night"


def compute_stats(mersi_clm, myd_clm):
    """Compute validation stats. mersi_clm: 0-1=cloud, 2-3=clear. myd_clm: same."""
    valid = (mersi_clm >= 0) & (myd_clm >= 0)
    n = valid.sum()
    if n < 100:
        return None

    m_cloud = (mersi_clm[valid] <= 1)
    y_cloud = (myd_clm[valid] <= 1)

    tp = (m_cloud & y_cloud).sum()
    fp = (m_cloud & ~y_cloud).sum()
    fn = (~m_cloud & y_cloud).sum()
    tn = (~m_cloud & ~y_cloud).sum()

    pod = tp / (tp + fn + 1e-9)
    far = fp / (tp + fp + 1e-9)
    agree = (tp + tn) / n
    csi = tp / (tp + fp + fn + 1e-9)
    expected = ((tp + fp) * (tp + fn) + (tn + fp) * (tn + fn)) / (n + 1e-9)
    hss = (tp + tn - expected) / (n - expected + 1e-9)

    return {
        "n": int(n), "agree": agree, "pod": pod, "far": far,
        "csi": csi, "hss": hss,
        "mersi_cloud_pct": m_cloud.mean(),
        "myd35_cloud_pct": y_cloud.mean(),
    }


def main():
    parser = argparse.ArgumentParser(description="Batch MYD35 overlap analysis")
    parser.add_argument("--date", default="20220803")
    parser.add_argument("--data_dir", default="/data/Data_yuq/fy3_cloud/")
    parser.add_argument("--myd35_dir", default="/data/Data_yuq/aqua_modis/MYD35_L2/")
    parser.add_argument("--mersi_root", default="/data/Data_yuq/mersi/")
    parser.add_argument("--time_window", type=int, default=15)
    args = parser.parse_args()

    clm_dir = os.path.join(args.data_dir, args.date)
    recal_files = sorted(Path(clm_dir).glob("*_RECALI.HDF"))

    if not recal_files:
        print(f"No RECALI files found in {clm_dir}")
        return 1

    print(f"Analyzing {len(recal_files)} time slots for {args.date}")
    print()

    rows = []
    no_match = []

    for rpath in recal_files:
        fname = rpath.name
        ts_match = re.search(r'_(\d{4})_', fname)
        ts = ts_match.group(1) if ts_match else "????"

        data = load_clm_hdf5(str(rpath))
        if data is None:
            no_match.append((ts, "LOAD_FAIL"))
            continue

        lat = data["lat"]
        lon = data["lon"]
        clm = data["clm"]

        dt = parse_mersi_datetime(str(rpath))
        if dt is None:
            no_match.append((ts, "NO_DATETIME"))
            continue

        myd = load_best_myd35_for_mersi(
            mersi_lat=lat, mersi_lon=lon, mersi_dt=dt,
            search_dirs=[os.path.join(args.myd35_dir, args.date)],
            time_window_min=args.time_window, min_overlap=0.05,
        )
        if myd is None:
            no_match.append((ts, "NO_MYD35_OVERLAP"))
            continue

        myd_clm = myd["clm_resampled"]
        overlap_pct = (myd_clm >= 0).sum() / clm.size * 100
        delta_t = abs((dt - myd["dt"]).total_seconds() / 60) if myd.get("dt") else 0

        stats = compute_stats(clm, myd_clm)
        if stats is None:
            no_match.append((ts, "TOO_FEW_PIXELS"))
            continue

        # Rough day/night from lat (polar night check)
        abs_lat = np.abs(lat)
        polar = (abs_lat >= 60).mean()
        rows.append({
            "ts": ts,
            "dt_min": delta_t,
            "overlap_pct": overlap_pct,
            "granule": os.path.basename(myd.get("path", "?")),
            **stats,
        })

    # Sort by HSS descending
    rows.sort(key=lambda r: r["hss"], reverse=True)

    print(f"{'Time':>6s}  {'dt':>5s}  {'Overlap':>7s}  {'N':>8s}  "
          f"{'Agree':>7s}  {'POD':>7s}  {'FAR':>7s}  {'HSS':>8s}  "
          f"{'M_Cloud':>8s}  {'Y_Cloud':>8s}  {'Bias':>7s}  Granule")
    print("-" * 130)

    for r in rows:
        bias = r["mersi_cloud_pct"] * 100 - r["myd35_cloud_pct"] * 100
        print(f"{r['ts']:>6s}  {r['dt_min']:4.0f}m  {r['overlap_pct']:6.1f}%  {r['n']:>8,}  "
              f"{r['agree']*100:6.2f}%  {r['pod']*100:6.2f}%  {r['far']*100:6.2f}%  "
              f"{r['hss']:8.4f}  {r['mersi_cloud_pct']*100:7.1f}%  "
              f"{r['myd35_cloud_pct']*100:7.1f}%  {bias:+6.1f}%  "
              f"{r['granule'][:45]}")

    # Summary
    print(f"\n--- Summary ---")
    n = len(rows)
    avg_hss = np.mean([r["hss"] for r in rows])
    avg_agree = np.mean([r["agree"] for r in rows])
    n_bad = sum(1 for r in rows if r["hss"] < 0.3)
    n_good = sum(1 for r in rows if r["hss"] > 0.5)
    print(f"  Matched: {n}/{len(recal_files)} time slots")
    print(f"  Mean HSS: {avg_hss:.4f}")
    print(f"  Mean Agreement: {avg_agree*100:.1f}%")
    print(f"  Good (HSS>0.5): {n_good} slots")
    print(f"  Poor (HSS<0.3): {n_bad} slots")

    if no_match:
        print(f"\n  No match ({len(no_match)}):")
        for ts, reason in no_match:
            print(f"    {ts}: {reason}")

    # Top/bottom 5
    print(f"\n--- Best 5 ---")
    for r in rows[-5:]:
        print(f"  {r['ts']}: HSS={r['hss']:.4f}  POD={r['pod']*100:.1f}%  FAR={r['far']*100:.1f}%  N={r['n']:,}")
    print(f"\n--- Worst 5 ---")
    for r in rows[:5]:
        print(f"  {r['ts']}: HSS={r['hss']:.4f}  POD={r['pod']*100:.1f}%  FAR={r['far']*100:.1f}%  N={r['n']:,}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
