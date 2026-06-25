"""Visible band calibration coefficients for FY-3D MERSI-II.

Auto-discovers recalibration files from ../fy3d_recali/ directory.
No more hardcoded coefficients or manual .xcfg management.

Usage:
    from fylat.calibration import list_calibrations, setup_cal_mode
    print(list_calibrations())       # ['business', '20220803', '20200308', ...]
    setup_cal_mode("20220803", ".")  # writes cal_mode.txt + VIS_Cal_Coeff.xcfg
    setup_cal_mode("business", ".")  # use HDF5 built-in coefficients
"""

import os
import glob
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np

# Band indices 7-25: same for all calibration sets (thermal IR — not calibrated)
_COMMON_IR = [
    (7, -1.6386, 0.0101), (8, -1.606, 0.008869), (9, -1.147, 0.00843),
    (10, -1.856, 0.01083), (11, -1.735, 0.009302), (12, -1.739, 0.008765),
    (13, -1.693, 0.008661), (14, -1.299, 0.008634), (15, -4.5605, 0.0261),
    (16, -5.6694, 0.0242), (17, -6.784, 0.0286), (18, -6.0429, 0.0335),
    (19, 0, 0), (20, 0, 0), (21, 0, 0), (22, 0, 0), (23, 0, 0),
    (24, 0, 0), (25, 0, 0),
]

# Recalibration data directory (relative to project root)
_RECALI_DIR = os.path.join(os.path.dirname(__file__), "..", "..", "..", "fy3d_recali")


def _find_recali_dir() -> str:
    """Find the recalibration data directory."""
    path = os.path.abspath(_RECALI_DIR)
    if os.path.isdir(path):
        return path
    # Try alternate locations
    alt = "/home/liusy2020/yuq/cloudmask/fy3d_recali"
    if os.path.isdir(alt):
        return alt
    return ""


def list_calibrations() -> List[str]:
    """List all available calibration options (business + auto-discovered dates)."""
    available = ["business"]
    recali_dir = _find_recali_dir()
    if recali_dir:
        for entry in sorted(os.listdir(recali_dir)):
            entry_path = os.path.join(recali_dir, entry)
            if os.path.isdir(entry_path) and len(entry) == 6:
                # Check if it contains .h5 files
                h5_files = glob.glob(os.path.join(entry_path, "RAD_*.h5"))
                if h5_files:
                    available.append(entry)
    return available


def load_coefficients(cal_name: str) -> List[Tuple[int, float, float]]:
    """Load visible band coefficients for a given calibration.

    Args:
        cal_name: "business" for HDF5 built-in, or "YYYYMM" for recalibration.

    Returns:
        List of (band_index, offset, slope) tuples for bands 0-25.
    """
    if cal_name == "business":
        return _COMMON_IR  # No VIS override — Fortran uses HDF5 built-in

    # Try to load from recalibration directory
    recali_dir = _find_recali_dir()
    if not recali_dir:
        raise FileNotFoundError(f"Recalibration directory not found: {_RECALI_DIR}")

    # Find the specific day file within the month directory
    month_dir = os.path.join(recali_dir, cal_name)
    if not os.path.isdir(month_dir):
        raise FileNotFoundError(f"No recalibration data for {cal_name}")

    h5_files = sorted(glob.glob(os.path.join(month_dir, "RAD_*.h5")))
    if not h5_files:
        raise FileNotFoundError(f"No .h5 files in {month_dir}")

    # Use the first available file for this month
    import h5py
    with h5py.File(h5_files[0], 'r') as f:
        coeffs = f['calibration_coeff'][:]  # shape (7, 3): offset, slope, 0

    vis_coeffs = []
    for i in range(7):
        offset = float(coeffs[i, 0])
        slope = float(coeffs[i, 1])
        vis_coeffs.append((i, offset, slope))

    return vis_coeffs + _COMMON_IR


def generate_xcfg(cal_name: str) -> str:
    """Generate VIS_Cal_Coeff.xcfg content."""
    coeffs = load_coefficients(cal_name) if cal_name != "business" else _COMMON_IR
    lines = ["VIS_Cal_Coeff = ["]
    for band, offset, slope in coeffs:
        lines.append(f" {band}, {offset}, {slope},")
    if len(lines) > 2:
        lines[-2] = lines[-2].rstrip(",")
    return "\n".join(lines) + "\n"


def setup_cal_mode(cal_name: str, code_root: str) -> None:
    """Configure calibration mode for the Fortran runtime.

    - "business": use HDF5 built-in coefficients (no .xcfg needed)
    - "YYYYMM": load from ../fy3d_recali/YYYYMM/RAD_*.h5
    """
    cal_mode_path = os.path.join(code_root, "cal_mode.txt")
    xcfg_path = os.path.join(code_root, "VIS_Cal_Coeff.xcfg")

    if cal_name == "business":
        if os.path.exists(cal_mode_path):
            os.remove(cal_mode_path)
        if os.path.exists(xcfg_path):
            os.remove(xcfg_path)
    else:
        with open(cal_mode_path, "w") as f:
            f.write("recali\n")
        content = generate_xcfg(cal_name)
        with open(xcfg_path, "w") as f:
            f.write(content)
