#!/usr/bin/env python3
"""Run FYLAT cloud mask for a single scene using YAML configuration.

Usage:
  python scripts/run_single.py config/scenes/20220803_0740.yaml
  python scripts/run_single.py config/scenes/20220803_0740.yaml --calibration recali
  python scripts/run_single.py --date 20220803 --time 0740 --calibration business

The script:
  1. Loads YAML config (defaults + scene overrides)
  2. Sets up calibration mode (cal_mode.txt + VIS_Cal_Coeff.xcfg)
  3. Generates a Fortran namelist (.nml) file
  4. Runs the fylat_FY3_MERSI_II_PGS executable
"""

import argparse
import os
import sys

# Add project python/ to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "python"))

from fylat.config import load_config, run_scene


def main():
    parser = argparse.ArgumentParser(
        description="Run FYLAT cloud mask for a single scene (YAML config)"
    )
    parser.add_argument(
        "scene_yaml", nargs="?", default=None,
        help="Path to scene YAML config file"
    )
    parser.add_argument(
        "--date", default=None,
        help="Observation date YYYYMMDD (auto-selects config/scenes/<date>_<time>.yaml)"
    )
    parser.add_argument(
        "--time", default=None,
        help="Observation time HHMM"
    )
    parser.add_argument(
        "--calibration", default=None,
        help="Calibration mode (business=onboard, YYYYMM=auto-load from ../fy3d_recali/)"
    )
    parser.add_argument(
        "--list-calibrations", action="store_true",
        help="List all available calibration options and exit"
    )
    parser.add_argument(
        "--base-config", default=None,
        help="Path to base YAML config (default: config/default.yaml)"
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Generate .nml and setup calibration but don't run Fortran"
    )
    args = parser.parse_args()

    if args.list_calibrations:
        from fylat.calibration import list_calibrations
        cals = list_calibrations()
        print(f"Available calibrations ({len(cals)}):")
        for c in cals:
            print(f"  {c}")
        return 0

    # Determine scene YAML path
    if args.scene_yaml:
        scene_path = args.scene_yaml
    elif args.date and args.time:
        project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        scene_path = os.path.join(
            project_root, "config", "scenes", f"{args.date}_{args.time}.yaml"
        )
    else:
        parser.error("Either provide scene_yaml path or --date and --time")

    if not os.path.exists(scene_path):
        sys.exit(f"Scene config not found: {scene_path}")

    # Load and optionally override calibration
    cfg = load_config(scene_path=scene_path, base_path=args.base_config)
    if args.calibration:
        cfg.setdefault("scene", {})["calibration"] = args.calibration

    scene = cfg.get("scene", {})
    print(f"FYLAT Cloud Mask — {scene.get('date', '?')} {scene.get('time', '?')} "
          f"({scene.get('calibration', 'business')})")
    print(f"  Sensor ID: {cfg.get('sensor', {}).get('id', '?')}")
    print(f"  NWP opt:   {cfg.get('sensor', {}).get('nwp_opt', '?')}")
    print(f"  RTM opt:   {cfg.get('sensor', {}).get('rtm_opt', '?')}")

    if args.dry_run:
        from fylat.config import generate_namelist, setup_calibration
        code_root = cfg.get("paths", {}).get("code_root", ".")
        setup_calibration(cfg, code_root)
        nml = generate_namelist(cfg)
        print("\n--- Generated namelist ---")
        print(nml)
        print("--- End namelist ---")
        print("\nDry run complete. .nml NOT written to disk.")
        return 0

    ret = run_scene(scene_path, base_yaml_path=args.base_config,
                    calibration=args.calibration)
    if ret == 0:
        print("FYLAT completed successfully.")
    else:
        print(f"FYLAT exited with code {ret}")
    return ret


if __name__ == "__main__":
    sys.exit(main())
