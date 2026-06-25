#!/usr/bin/env python3
"""
run_fylat.py — FY-3D MERSI-II 云检测全链条自动化反演

一条命令完成: 发现时次 → NWP预处理 → 编译 → 双定标反演 → 验证

Usage:
  python run_fylat.py --date 20220803
  python run_fylat.py --date 20220803 --cores 4
  python run_fylat.py --date 20220803 --skip-build
  python run_fylat.py --date 20220803 --verify-only
  python run_fylat.py --date 20220803 --dry-run       # 仅生成配置不运行
"""

import argparse
import glob
import os
import shutil
import subprocess
import sys
import tempfile
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# ---------------------------------------------------------------------------
# Default paths — override via environment or command line
# ---------------------------------------------------------------------------
PROJECT_ROOT = os.path.dirname(os.path.abspath(__file__))
DEFAULT_L1B_PATH = "/data/Data_yuq/mersi/"
DEFAULT_NWP_PATH = "/data/nwp/"
DEFAULT_OISST_PATH = "/data/Data_minmin/oisst/"
DEFAULT_OUTPUT_PATH = "/data/Data_yuq/fy3_cloud/"
DEFAULT_OISST_FILE = "sst.day.mean.20200401.hdf5"


# ---------------------------------------------------------------------------
# NWP mapping utilities
# ---------------------------------------------------------------------------
def _int2str(v: int) -> str:
    return f"{v:02d}"


def find_nwp_forecast_hours(obs_hour: int) -> Tuple[int, int]:
    """Map observation hour to GFS forecast lead times (+18h offset).

    NWP cycles: [0, 3, 6, 9, 12, 15, 18, 21, 24]
    Forecast lead = cycle_hour + 18 (from batch_run.py convention).
    """
    nwp_cycles = [0, 3, 6, 9, 12, 15, 18, 21, 24]
    for i in range(len(nwp_cycles) - 1):
        if nwp_cycles[i] <= obs_hour < nwp_cycles[i + 1]:
            n1 = nwp_cycles[i] + 18
            n2 = nwp_cycles[i + 1] + 18
            return n1, n2
    return 18, 21  # fallback


def discover_nwp_files(nwp_path: str, date: str) -> Tuple[str, str]:
    """Auto-discover GFS GRIB2 files for a given date.

    Returns (grib1, grib2) paths suitable for nwp_grib_data1/2 in .nml.
    """
    org_dir = os.path.join(nwp_path, date, "ORG")
    if not os.path.isdir(org_dir):
        raise FileNotFoundError(f"NWP ORG directory not found: {org_dir}")

    gfs_files = sorted(glob.glob(os.path.join(org_dir, "gfs.t*z.pgrb2.0p25.f*")))
    # Filter out .ok marker files
    gfs_files = [f for f in gfs_files if not f.endswith(".ok")]

    if not gfs_files:
        raise FileNotFoundError(f"No GFS GRIB2 files found in {org_dir}")

    # Group by cycle (t06z, t12z, etc.)
    cycles: Dict[str, List[str]] = {}
    for f in gfs_files:
        basename = os.path.basename(f)
        # gfs.t06z.pgrb2.0p25.f018
        parts = basename.split(".")
        if len(parts) >= 4:
            cycle = parts[1]  # t06z
            if cycle not in cycles:
                cycles[cycle] = []
            cycles[cycle].append(f)

    if not cycles:
        raise FileNotFoundError(f"Cannot parse GFS cycle from files in {org_dir}")

    # Use the earliest available cycle (typically 00z or 06z)
    chosen_cycle = sorted(cycles.keys())[0]
    cycle_files = sorted(cycles[chosen_cycle])

    # Pick f018 and f021 as defaults (matching existing working configs)
    f018 = None
    f021 = None
    for f in cycle_files:
        basename = os.path.basename(f)
        if "f018" in basename:
            f018 = f
        elif "f021" in basename:
            f021 = f

    if f018 is None or f021 is None:
        # Fallback: use first two available
        if len(cycle_files) >= 2:
            f018 = cycle_files[0]
            f021 = cycle_files[1]

    if f018 is None or f021 is None:
        raise FileNotFoundError(f"Need at least 2 GFS forecast files, found {len(cycle_files)}")

    print(f"  NWP cycle: {chosen_cycle}, files: {os.path.basename(f018)}, {os.path.basename(f021)}")
    return f018, f021


# ---------------------------------------------------------------------------
# Time slot discovery
# ---------------------------------------------------------------------------
def discover_time_slots(l1b_path: str, date: str) -> List[str]:
    """Find all observation time slots (HHMM) with valid GEO + L1B pairs."""
    data_dir = os.path.join(l1b_path, date)
    if not os.path.isdir(data_dir):
        raise FileNotFoundError(f"L1B data directory not found: {data_dir}")

    geo_pattern = os.path.join(data_dir, f"FY3D_MERSI_GBAL_L1_{date}_*_GEO1K_MS.HDF")
    l1b_pattern = os.path.join(data_dir, f"FY3D_MERSI_GBAL_L1_{date}_*_1000M_MS.HDF")

    geo_files = sorted(glob.glob(geo_pattern))
    l1b_files = sorted(glob.glob(l1b_pattern))

    # Extract time slots from GEO files (HHMM)
    geo_times = set()
    for f in geo_files:
        basename = os.path.basename(f)
        # FY3D_MERSI_GBAL_L1_20220803_0740_GEO1K_MS.HDF
        # parts: [FY3D, MERSI, GBAL, L1, 20220803, 0740, GEO1K, MS.HDF]
        parts = basename.split("_")
        if len(parts) >= 6:
            geo_times.add(parts[5])

    l1b_times = set()
    for f in l1b_files:
        basename = os.path.basename(f)
        parts = basename.split("_")
        if len(parts) >= 6:
            l1b_times.add(parts[5])

    # Only keep time slots with both GEO and L1B
    valid_times = sorted(geo_times & l1b_times)
    if not valid_times:
        raise FileNotFoundError(f"No valid GEO+L1B pairs found for date {date}")

    print(f"  Discovered {len(valid_times)} time slot(s): {', '.join(valid_times)}")
    return valid_times


# ---------------------------------------------------------------------------
# OISST file discovery
# ---------------------------------------------------------------------------
def discover_oisst_file(oisst_path: str) -> str:
    """Find an OISST daily file."""
    hdf5_files = sorted(glob.glob(os.path.join(oisst_path, "sst.day.mean.*.hdf5")))
    if hdf5_files:
        return hdf5_files[0]
    # Fallback
    fallback = os.path.join(oisst_path, DEFAULT_OISST_FILE)
    if os.path.exists(fallback):
        return fallback
    raise FileNotFoundError(f"No OISST file found in {oisst_path}")


# ---------------------------------------------------------------------------
# Calibration setup
# ---------------------------------------------------------------------------
def setup_calibration(cal_name: str, work_dir: str, code_root: str, date: str) -> str:
    """Set up calibration mode files.

    Fortran reads:
      - cal_mode.txt from code_root_path/ (line 430 in io_module.f90)
      - VIS_Cal_Coeff.xcfg from CWD (line 458 in io_module.f90)

    - "business": remove any calibration override files
    - "recali": auto-map to YYYYMM based on observation date
    - "YYYYMM": explicit recalibration month

    Returns the actual calibration name used.
    """
    # cal_mode.txt lives in code_root (not CWD)
    cal_mode_path = os.path.join(code_root, "cal_mode.txt")
    xcfg_path = os.path.join(work_dir, "VIS_Cal_Coeff.xcfg")

    # Clean up any existing files
    for p in [cal_mode_path, xcfg_path]:
        if os.path.exists(p):
            os.remove(p)

    if cal_name == "business":
        return "business"

    # Map "recali" to date-based YYYYMM
    sys.path.insert(0, os.path.join(PROJECT_ROOT, "python"))
    from fylat.calibration import list_calibrations, generate_xcfg

    available = list_calibrations()

    if cal_name == "recali":
        # Auto-map to YYYYMM based on observation date
        cal_month = date[:6]  # 20220803 → 202208
        if cal_month in available:
            cal_name = cal_month
        else:
            # Try to find any available recalibration for that year
            year_cals = [c for c in available if c.startswith(date[:4])]
            if year_cals:
                cal_name = year_cals[0]
                print(f"  Calibration {date[:6]} not found, using {cal_name}")
            else:
                print(f"  No recalibration data for {date[:4]}, using 'business'")
                return "business"

    if cal_name not in available:
        print(f"  Warning: calibration '{cal_name}' not found, using 'business'")
        return "business"

    with open(cal_mode_path, "w") as f:
        f.write("recali\n")

    content = generate_xcfg(cal_name)
    with open(xcfg_path, "w") as f:
        f.write(content)

    print(f"  Calibration: {cal_name} (RECALI mode)")
    return cal_name


# ---------------------------------------------------------------------------
# Namelist generation
# ---------------------------------------------------------------------------
def generate_nml(
    date: str,
    time_slot: str,
    calibration: str,
    l1b_path: str,
    nwp_path: str,
    oisst_path: str,
    output_path: str,
    code_root: str,
    nwp_grib1: str,
    nwp_grib2: str,
    oisst_file: str,
) -> str:
    """Generate Fortran namelist content for a single scene."""
    cal_suffix = "BUSINESS" if calibration == "business" else "RECALI"
    out_dir = os.path.join(output_path, date)
    os.makedirs(out_dir, exist_ok=True)

    geo_file = os.path.join(l1b_path, date,
                            f"FY3D_MERSI_GBAL_L1_{date}_{time_slot}_GEO1K_MS.HDF")
    l1b_file = os.path.join(l1b_path, date,
                            f"FY3D_MERSI_GBAL_L1_{date}_{time_slot}_1000M_MS.HDF")

    base = f"FY3D_MERSI_ORBT_L2_{{}}_MLT_NUL_{date}_{time_slot}_1000M_MS_{cal_suffix}.HDF"
    base_5km = f"FY3D_MERSI_ORBT_L2_{{}}_MLT_NUL_{date}_{time_slot}_5000M_MS_{cal_suffix}.HDF"

    lines = ["&config"]
    lines.append(f"  fylat_sensor_id     = 21,")
    lines.append(f'  code_root_path      = "{code_root}/",')
    lines.append(f'  L1b_data_path       = "{os.path.join(l1b_path, date)}/",')
    lines.append(f'  nwp_data_path       = "{os.path.join(nwp_path, date)}/",')
    lines.append(f'  oisst_data_path     = "{oisst_path}",')
    lines.append(f'  fy3_mersi_GEO_data  = "{geo_file}",')
    lines.append(f'  fy3_mersi_L1b_data  = "{l1b_file}",')
    lines.append(f'  fy3_mersi_CLM_data  = "{os.path.join(out_dir, base.format("CLM"))}",')
    lines.append(f'  fy3_mersi_CLA_data  = "{os.path.join(out_dir, base_5km.format("CLA"))}",')
    lines.append(f'  fy3_mersi_CLP_data  = "{os.path.join(out_dir, base.format("CLP"))}",')
    lines.append(f'  fy3_mersi_CTP_data  = "{os.path.join(out_dir, base.format("CTP"))}",')
    lines.append(f'  fy3_mersi_COT_data  = "{os.path.join(out_dir, base.format("COT"))}",')
    lines.append(f'  fy3_mersi_CON_data  = "{os.path.join(out_dir, base.format("CON"))}",')
    lines.append(f'  fy3_mersi_SST_data  = "{os.path.join(out_dir, base.format("SST"))}",')
    lines.append(f'  fy3_intermediate    = "{os.path.join(out_dir, f"FY3D_MERSI_ORBT_L2_XXX_MLT_NUL_{date}_{time_slot}_INTERMED_{cal_suffix}.HDF")}",')
    lines.append(f"  fylat_nwp_opt       = 10,")
    lines.append(f"  fylat_rtm_opt       = 1,")
    lines.append(f'  nwp_grib_data1      = "{nwp_grib1}",')
    lines.append(f'  nwp_grib_data2      = "{nwp_grib2}",')
    lines.append(f'  oisst_data          = "{oisst_file}",')
    lines.append(f"  cloudmask_id        = 1,")
    lines.append(f"  cloudamount_id      = 0,")
    lines.append(f"  cloudphase_id       = 0,")
    lines.append(f"  cloudtopz_id        = 0,")
    lines.append(f"  cloudtau_day_id     = 0,")
    lines.append(f"  cloudtau_night_id   = 0,")
    lines.append(f"  cloudtypeII_id      = 0,")
    lines.append(f"  surface_sst_id      = 0,")
    lines.append(f"  write_inter_id      = 0/")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Build
# ---------------------------------------------------------------------------
def build_executable() -> None:
    """Run build.sh to compile the Fortran executable."""
    build_script = os.path.join(PROJECT_ROOT, "build.sh")
    if not os.path.exists(build_script):
        raise FileNotFoundError(f"Build script not found: {build_script}")

    print("\n" + "=" * 60)
    print("  BUILD: Compiling Fortran executable")
    print("=" * 60)

    result = subprocess.run(
        ["bash", "build.sh"],
        cwd=PROJECT_ROOT,
        capture_output=False,
    )

    if result.returncode != 0:
        print("  BUILD FAILED. Check output above for errors.")
        sys.exit(1)

    exe_path = os.path.join(PROJECT_ROOT, "fylat_FY3_MERSI_II_PGS")
    if not os.path.exists(exe_path):
        print(f"  BUILD FAILED: executable not found at {exe_path}")
        sys.exit(1)

    print("  BUILD: OK")


# ---------------------------------------------------------------------------
# Single scene execution
# ---------------------------------------------------------------------------
def run_single_scene(args: Tuple) -> Tuple[str, str, str, int, str]:
    """Run a single Fortran execution. Designed for parallel invocation.

    Args:
        args: (date, time_slot, calibration, l1b_path, nwp_path, oisst_path,
               output_path, code_root, nwp_grib1, nwp_grib2, oisst_file)

    Returns:
        (date, time_slot, calibration, returncode, stdout_log)
    """
    (date, time_slot, calibration,
     l1b_path, nwp_path, oisst_path, output_path,
     code_root, nwp_grib1, nwp_grib2, oisst_file) = args

    label = f"{date}_{time_slot}_{calibration}"

    # Create isolated working directory
    work_dir = tempfile.mkdtemp(prefix=f"fylat_{label}_")
    log_path = os.path.join(work_dir, "run.log")

    try:
        # Set up calibration mode (returns resolved calibration name)
        actual_cal = setup_calibration(calibration, work_dir, code_root, date)

        # Generate .nml config in work_dir (use original cal name for file suffix)
        nml_content = generate_nml(
            date, time_slot, calibration,
            l1b_path, nwp_path, oisst_path, output_path,
            code_root, nwp_grib1, nwp_grib2, oisst_file,
        )
        nml_path = os.path.join(work_dir, f"config_{label}.nml")
        with open(nml_path, "w") as f:
            f.write(nml_content + "\n")

        # Run Fortran
        exe_path = os.path.join(code_root, "fylat_FY3_MERSI_II_PGS")
        if not os.path.exists(exe_path):
            return (date, time_slot, calibration, -1,
                    f"Executable not found: {exe_path}")

        t0 = time.time()
        with open(log_path, "w") as log_f:
            result = subprocess.run(
                [exe_path, nml_path],
                cwd=work_dir,
                stdout=log_f,
                stderr=subprocess.STDOUT,
                timeout=7200,  # 2 hour timeout per scene
            )
        elapsed = time.time() - t0

        if result.returncode == 0:
            status = f"OK ({elapsed:.0f}s)"
        else:
            status = f"FAILED (rc={result.returncode}, {elapsed:.0f}s)"

        return (date, time_slot, calibration, result.returncode, status)

    except subprocess.TimeoutExpired:
        return (date, time_slot, calibration, -1, "TIMEOUT (>2h)")
    except Exception as e:
        return (date, time_slot, calibration, -1, f"ERROR: {e}")
    finally:
        # Clean up temp directory (keep on failure for debugging)
        pass
        # Uncomment to auto-clean: shutil.rmtree(work_dir, ignore_errors=True)


# ---------------------------------------------------------------------------
# Output verification
# ---------------------------------------------------------------------------
def verify_outputs(
    date: str,
    time_slots: List[str],
    calibrations: List[str],
    output_path: str,
) -> Dict[str, bool]:
    """Verify output files exist and have reasonable sizes."""
    print("\n" + "=" * 60)
    print("  VERIFY: Checking output files")
    print("=" * 60)

    results = {}

    for ts in time_slots:
        for cal in calibrations:
            cal_suffix = "BUSINESS" if cal == "business" else "RECALI"
            clm_file = os.path.join(
                output_path, date,
                f"FY3D_MERSI_ORBT_L2_CLM_MLT_NUL_{date}_{ts}_1000M_MS_{cal_suffix}.HDF"
            )
            label = f"{date}_{ts}_{cal}"

            if not os.path.exists(clm_file):
                print(f"  [FAIL] {label}: CLM file not found: {clm_file}")
                results[label] = False
                continue

            size_mb = os.path.getsize(clm_file) / (1024 * 1024)
            if size_mb < 1.0:
                print(f"  [WARN] {label}: CLM file suspiciously small ({size_mb:.1f} MB)")
                results[label] = False
            else:
                print(f"  [OK]   {label}: CLM {size_mb:.1f} MB")
                results[label] = True

    return results


# ---------------------------------------------------------------------------
# Main driver
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="FYLAT FY-3D MERSI-II 云检测全链条自动化反演",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python run_fylat.py --date 20220803
  python run_fylat.py --date 20220803 --cores 4
  python run_fylat.py --date 20220803 --dry-run
  python run_fylat.py --date 20220803 --skip-build --verify-only
        """,
    )
    parser.add_argument("--date", required=True, help="观测日期 YYYYMMDD")
    parser.add_argument("--cores", type=int, default=None,
                        help="并行核心数 (默认: 全部可用)")
    parser.add_argument("--l1b-path", default=DEFAULT_L1B_PATH,
                        help="L1B 数据根目录")
    parser.add_argument("--nwp-path", default=DEFAULT_NWP_PATH,
                        help="NWP 数据根目录")
    parser.add_argument("--oisst-path", default=DEFAULT_OISST_PATH,
                        help="OISST 数据目录")
    parser.add_argument("--output-path", default=DEFAULT_OUTPUT_PATH,
                        help="输出根目录")
    parser.add_argument("--calibrations", default="business,recali",
                        help="定标模式，逗号分隔 (default: business,recali)")
    parser.add_argument("--skip-nwp", action="store_true",
                        help="跳过 NWP 预处理")
    parser.add_argument("--skip-build", action="store_true",
                        help="跳过编译 (假设已编译)")
    parser.add_argument("--verify-only", action="store_true",
                        help="仅验证已有输出，不运行反演")
    parser.add_argument("--dry-run", action="store_true",
                        help="仅生成配置文件和定标设置，不运行")

    args = parser.parse_args()

    # Validate date format
    try:
        datetime.strptime(args.date, "%Y%m%d")
    except ValueError:
        sys.exit(f"Invalid date format: {args.date}. Expected YYYYMMDD.")

    code_root = PROJECT_ROOT
    calibrations = [c.strip() for c in args.calibrations.split(",")]

    # =====================================================================
    # Step 1: Discover time slots
    # =====================================================================
    print("\n" + "=" * 60)
    print(f"  STEP 1: Discover data for {args.date}")
    print("=" * 60)

    time_slots = discover_time_slots(args.l1b_path, args.date)

    # =====================================================================
    # Step 2: Discover NWP files
    # =====================================================================
    print("\n  STEP 2: Discover NWP files")
    nwp_grib1, nwp_grib2 = discover_nwp_files(args.nwp_path, args.date)

    # =====================================================================
    # Step 3: Discover OISST
    # =====================================================================
    print("\n  STEP 3: Discover OISST file")
    oisst_file = discover_oisst_file(args.oisst_path)
    print(f"  OISST: {os.path.basename(oisst_file)}")

    # =====================================================================
    # Step 4: Build
    # =====================================================================
    if not args.skip_build and not args.verify_only and not args.dry_run:
        build_executable()

    # =====================================================================
    # Step 5: Dry run (generate configs only)
    # =====================================================================
    if args.dry_run:
        print("\n" + "=" * 60)
        print("  DRY RUN: Generating configuration files")
        print("=" * 60)
        for ts in time_slots:
            for cal in calibrations:
                nml = generate_nml(
                    args.date, ts, cal,
                    args.l1b_path, args.nwp_path, args.oisst_path,
                    args.output_path, code_root,
                    nwp_grib1, nwp_grib2, oisst_file,
                )
                nml_file = os.path.join(
                    code_root, f"temp_fy3d_config_{args.date}_{ts}_{cal}.nml"
                )
                with open(nml_file, "w") as f:
                    f.write(nml + "\n")
                print(f"  Written: {nml_file}")
        print("\n  Dry run complete.")
        return 0

    # =====================================================================
    # Step 5: Verify only
    # =====================================================================
    if args.verify_only:
        results = verify_outputs(args.date, time_slots, calibrations, args.output_path)
        all_ok = all(results.values()) if results else False
        if all_ok:
            print("\n  VERIFY: ALL OUTPUTS OK")
        else:
            print("\n  VERIFY: SOME OUTPUTS MISSING OR INVALID")
        return 0 if all_ok else 1

    # =====================================================================
    # Step 6: Run all scenes in parallel
    # =====================================================================
    tasks = []
    for ts in time_slots:
        for cal in calibrations:
            tasks.append((
                args.date, ts, cal,
                args.l1b_path, args.nwp_path, args.oisst_path,
                args.output_path, code_root,
                nwp_grib1, nwp_grib2, oisst_file,
            ))

    n_total = len(tasks)
    n_cores = args.cores or min(os.cpu_count() or 1, n_total)

    print("\n" + "=" * 60)
    print(f"  STEP 6: Run {n_total} retrieval task(s) on {n_cores} core(s)")
    print("=" * 60)
    for ts in time_slots:
        for cal in calibrations:
            print(f"  Task: {args.date}_{ts}_{cal}")
    print()

    results_list: List[Tuple] = []

    if n_total == 1:
        # Sequential for single task
        results_list.append(run_single_scene(tasks[0]))
    else:
        with ProcessPoolExecutor(max_workers=n_cores) as executor:
            futures = {executor.submit(run_single_scene, t): t for t in tasks}
            for future in as_completed(futures):
                result = future.result()
                results_list.append(result)
                _, ts, cal, rc, status = result
                marker = "OK" if rc == 0 else "FAIL"
                print(f"  [{marker}] {args.date}_{ts}_{cal}: {status}")

    # =====================================================================
    # Step 7: Summary
    # =====================================================================
    print("\n" + "=" * 60)
    print("  SUMMARY")
    print("=" * 60)
    all_ok = True
    for _, ts, cal, rc, status in results_list:
        label = f"{args.date}_{ts}_{cal}"
        if rc == 0:
            print(f"  [OK]   {label}: {status}")
        else:
            print(f"  [FAIL] {label}: {status}")
            all_ok = False

    # =====================================================================
    # Step 8: Verify outputs
    # =====================================================================
    verify_ok = True
    if all_ok:
        verify_results = verify_outputs(args.date, time_slots, calibrations, args.output_path)
        verify_ok = all(verify_results.values())

    if all_ok and verify_ok:
        print("\n  ALL DONE — Everything OK")
        return 0
    else:
        print("\n  DONE WITH ISSUES — check output above")
        return 1


if __name__ == "__main__":
    sys.exit(main())
