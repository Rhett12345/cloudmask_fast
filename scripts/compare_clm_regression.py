#!/usr/bin/env python3
"""Compare FYLAT CLM outputs for repeatable migration regression checks.

The script compares business/recalibration products for fixed time slots using:
  - whole-file SHA256
  - Cloud_Mask equality and differing-pixel count
  - Quality_Assurance equality and differing-pixel count

Use this after M1/M2 plumbing changes and before enabling any experimental
cloud-mask backend in production.
"""

from __future__ import annotations

import argparse
import hashlib
import os
from pathlib import Path
from typing import Iterable

import h5py
import numpy as np


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def product_path(output_root: Path, date: str, time_tag: str, calibration: str) -> Path:
    cal = calibration.upper()
    return (
        output_root
        / date
        / f"FY3D_MERSI_ORBT_L2_CLM_MLT_NUL_{date}_{time_tag}_1000M_MS_{cal}.HDF"
    )


def compare_dataset(path_a: Path, path_b: Path, dataset: str) -> tuple[bool, int, tuple[int, ...]]:
    with h5py.File(path_a, "r") as fa, h5py.File(path_b, "r") as fb:
        a = fa[dataset][:]
        b = fb[dataset][:]
    if a.shape != b.shape:
        return False, -1, tuple(a.shape)
    diff = int(np.count_nonzero(a != b))
    return diff == 0, diff, tuple(a.shape)


def iter_pairs(
    output_root: Path,
    date: str,
    times: Iterable[str],
    calibrations: Iterable[str],
    baseline_root: Path | None = None,
) -> Iterable[tuple[str, str, Path, Path]]:
    cals = list(calibrations)
    if baseline_root is not None:
        for time_tag in times:
            for cal in cals:
                yield (
                    time_tag,
                    f"baseline_vs_current_{cal}",
                    product_path(baseline_root, date, time_tag, cal),
                    product_path(output_root, date, time_tag, cal),
                )
        return

    if len(cals) == 1:
        return

    baseline = cals[0]
    for time_tag in times:
        base_path = product_path(output_root, date, time_tag, baseline)
        for cal in cals[1:]:
            yield time_tag, f"{baseline}_vs_{cal}", base_path, product_path(output_root, date, time_tag, cal)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--date", default="20220803", help="Observation date YYYYMMDD")
    parser.add_argument("--times", default="0715,0720,0740",
                        help="Comma-separated HHMM time slots")
    parser.add_argument("--calibrations", default="business,recali",
                        help="Comma-separated calibration products; first is baseline")
    parser.add_argument("--output-root", default="/data/Data_yuq/fy3_cloud",
                        help="Root containing date/product HDF5 outputs")
    parser.add_argument("--baseline-root", default=None,
                        help="Optional baseline output root for same-calibration old/new comparison")
    parser.add_argument("--require-identical", action="store_true",
                        help="Return non-zero when any compared dataset or file SHA differs")
    args = parser.parse_args()

    output_root = Path(args.output_root)
    baseline_root = Path(args.baseline_root) if args.baseline_root else None
    times = [x.strip() for x in args.times.split(",") if x.strip()]
    calibrations = [x.strip().lower() for x in args.calibrations.split(",") if x.strip()]
    datasets = ["Cloud_Mask", "Quality_Assurance"]

    any_missing = False
    any_diff = False

    for time_tag, label, path_a, path_b in iter_pairs(
        output_root, args.date, times, calibrations, baseline_root=baseline_root,
    ):
        print(f"\n[{args.date}_{time_tag}] {label}")
        print(f"  A: {path_a}")
        print(f"  B: {path_b}")
        if not path_a.exists() or not path_b.exists():
            print("  MISSING")
            any_missing = True
            continue

        sha_a = sha256(path_a)
        sha_b = sha256(path_b)
        same_sha = sha_a == sha_b
        any_diff = any_diff or not same_sha
        print(f"  SHA256 equal: {same_sha}")
        if not same_sha:
            print(f"    A {sha_a}")
            print(f"    B {sha_b}")

        for dataset in datasets:
            same, diff, shape = compare_dataset(path_a, path_b, dataset)
            any_diff = any_diff or not same
            diff_label = "shape mismatch" if diff < 0 else f"{diff:,} differing values"
            print(f"  {dataset:<18s} equal={str(same):<5s} shape={shape} {diff_label}")

    if any_missing:
        return 2
    if args.require_identical and any_diff:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
