#!/usr/bin/env python3
"""
batch_run.py — Multi-date batch inversion for FYLAT cloud mask.

Usage:
  python batch_run.py --start 20180103 --end 20250302 --cores 8
  python batch_run.py --cores 8 --dry-run
"""

import argparse
import os
import subprocess
import sys
import tempfile
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path

PROJECT_ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, PROJECT_ROOT)

from run_fylat import (
    discover_time_slots, _build_nwp_map, discover_nwp_for_slot,
    discover_oisst_file, generate_nml, setup_calibration,
)

L1B_PATH = "/data/Data_yuq/mersi/"
NWP_PATH = "/data/nwp/"
OISST_PATH = "/data/Data_minmin/oisst/"
OUTPUT_PATH = "/data/Data_yuq/fy3_cloud/"
EXE_PATH = os.path.join(PROJECT_ROOT, "fylat_FY3_MERSI_II_PGS")


def find_valid_dates(start_date: str = None, end_date: str = None):
    """Find all dates with both MERSI L1B and NWP data."""
    mersi_dates = set()
    for d in os.listdir(L1B_PATH):
        p = os.path.join(L1B_PATH, d)
        if os.path.isdir(p) and len(d) == 8:
            try:
                datetime.strptime(d, "%Y%m%d")
                mersi_dates.add(d)
            except ValueError:
                pass

    nwp_dates = set()
    for d in os.listdir(NWP_PATH):
        p = os.path.join(NWP_PATH, d)
        if os.path.isdir(p) and len(d) == 8:
            nwp_dates.add(d)

    valid = sorted(mersi_dates & nwp_dates)

    if start_date:
        valid = [d for d in valid if d >= start_date]
    if end_date:
        valid = [d for d in valid if d <= end_date]

    return valid


def run_single_task(args):
    """Run one (date, time_slot, calibration) task."""
    (date, time_slot, calibration, nwp_grib1, nwp_grib2, oisst_file) = args
    label = f"{date}_{time_slot}_{calibration}"
    work_dir = tempfile.mkdtemp(prefix=f"fylat_{label}_")

    try:
        actual_cal = setup_calibration(calibration, work_dir, PROJECT_ROOT, date)
        nml_content = generate_nml(
            date, time_slot, calibration,
            L1B_PATH, NWP_PATH, OISST_PATH, OUTPUT_PATH,
            PROJECT_ROOT, nwp_grib1, nwp_grib2, oisst_file,
        )
        nml_path = os.path.join(work_dir, f"config_{label}.nml")
        with open(nml_path, "w") as f:
            f.write(nml_content + "\n")

        if not os.path.exists(EXE_PATH):
            return (date, time_slot, calibration, -1, "EXE_NOT_FOUND")

        t0 = time.time()
        log_path = os.path.join(work_dir, "run.log")
        with open(log_path, "w") as log_f:
            result = subprocess.run(
                [EXE_PATH, nml_path],
                cwd=work_dir,
                stdout=log_f,
                stderr=subprocess.STDOUT,
                timeout=7200,
            )
        elapsed = time.time() - t0

        if result.returncode == 0:
            # Verify output
            cal_suffix = "BUSINESS" if calibration == "business" else "RECALI"
            clm_file = os.path.join(
                OUTPUT_PATH, date,
                f"FY3D_MERSI_ORBT_L2_CLM_MLT_NUL_{date}_{time_slot}_1000M_MS_{cal_suffix}.HDF"
            )
            if os.path.exists(clm_file):
                size_mb = os.path.getsize(clm_file) / (1024 * 1024)
                return (date, time_slot, calibration, 0, f"OK {elapsed:.0f}s {size_mb:.0f}MB")
            else:
                return (date, time_slot, calibration, 1, f"NoOutput {elapsed:.0f}s")
        else:
            return (date, time_slot, calibration, result.returncode, f"FAIL rc={result.returncode} {elapsed:.0f}s")

    except subprocess.TimeoutExpired:
        return (date, time_slot, calibration, -2, "TIMEOUT")
    except Exception as e:
        return (date, time_slot, calibration, -3, f"ERROR:{e}")
    finally:
        import shutil
        shutil.rmtree(work_dir, ignore_errors=True)


def main():
    parser = argparse.ArgumentParser(description="Batch FYLAT inversion over multiple dates")
    parser.add_argument("--start", help="Start date YYYYMMDD")
    parser.add_argument("--end", help="End date YYYYMMDD")
    parser.add_argument("--cores", type=int, default=8, help="Parallel cores (default: 8)")
    parser.add_argument("--dry-run", action="store_true", help="Only list tasks, don't run")
    args = parser.parse_args()

    if not os.path.exists(EXE_PATH):
        sys.exit(f"Executable not found: {EXE_PATH}. Run build.sh first.")

    # Find valid dates
    dates = find_valid_dates(args.start, args.end)
    if not dates:
        sys.exit("No valid dates found with both MERSI and NWP data.")

    # Discover tasks
    all_tasks = []
    skipped_dates = []
    for date in dates:
        try:
            time_slots = discover_time_slots(L1B_PATH, date)
            nwp_map = _build_nwp_map(NWP_PATH, date)
            oisst_file = discover_oisst_file(OISST_PATH)
        except (FileNotFoundError, Exception) as e:
            skipped_dates.append((date, str(e)))
            continue

        for ts in time_slots:
            nwp1, nwp2 = discover_nwp_for_slot(NWP_PATH, date, ts, nwp_map)
            for cal in ["business", "recali"]:
                all_tasks.append((date, ts, cal, nwp1, nwp2, oisst_file))

    n_dates = len(dates)
    n_tasks = len(all_tasks)
    print(f"Dates: {n_dates}  |  Tasks: {n_tasks}  |  Cores: {args.cores}")
    if skipped_dates:
        print(f"Skipped {len(skipped_dates)} dates (no NWP/OISST):")
        for d, reason in skipped_dates[:10]:
            print(f"  {d}: {reason}")
        if len(skipped_dates) > 10:
            print(f"  ... and {len(skipped_dates) - 10} more")

    if args.dry_run:
        print("\n[TASKS]")
        for i, t in enumerate(all_tasks):
            print(f"  {i+1}. {t[0]}_{t[1]}_{t[2]}")
        print(f"\nTotal: {n_tasks} tasks, ~{n_tasks * 80 / args.cores / 60:.0f} min estimated")
        return 0

    print(f"Estimated: ~{n_tasks * 80 / args.cores / 60:.0f} min\n")

    # Run
    t_start = time.time()
    completed = 0
    failed = 0
    results = []

    with ProcessPoolExecutor(max_workers=args.cores) as executor:
        futures = {executor.submit(run_single_task, t): t for t in all_tasks}
        for future in as_completed(futures):
            result = future.result()
            results.append(result)
            date, ts, cal, rc, status = result
            if rc == 0:
                completed += 1
            else:
                failed += 1
            elapsed = time.time() - t_start
            rate = (completed + failed) / elapsed * 60
            eta = (n_tasks - completed - failed) / rate if rate > 0 else 0
            marker = "OK" if rc == 0 else "FAIL"
            print(f"  [{marker}] {date}_{ts}_{cal}: {status}  "
                  f"[{completed}/{failed}] {rate:.1f}/min ETA {eta:.0f}min")

    # Summary
    elapsed_total = time.time() - t_start
    print(f"\n{'='*60}")
    print(f"  BATCH DONE — {elapsed_total/60:.0f} min")
    print(f"  OK: {completed}  |  FAIL: {failed}  |  Total: {n_tasks}")
    print(f"  Rate: {n_tasks / elapsed_total * 60:.1f} tasks/min")

    # List failures
    failures = [r for r in results if r[3] != 0]
    if failures:
        print(f"\n  FAILED ({len(failures)}):")
        for date, ts, cal, rc, status in failures:
            print(f"    {date}_{ts}_{cal}: {status}")

    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
