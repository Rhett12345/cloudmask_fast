"""Visible band calibration coefficients for FY-3D MERSI-II.

Replaces the manual .xcfg file management. Coefficients are stored in Python
dicts and written to VIS_Cal_Coeff.xcfg at runtime before Fortran execution.

Usage:
    from fylat.calibration import write_cal_xcfg, CALIBRATIONS
    write_cal_xcfg("business")     # use business (onboard) coeffs
    write_cal_xcfg("recali")      # use recali (external) coeffs
    write_cal_xcfg("20200308")    # use date-specific recalibration
"""

import os
from typing import Dict, List, Tuple

# Band indices 7-25: same for all calibration sets (thermal IR — not calibrated here)
_COMMON_IR = [
    (7, -1.6386, 0.0101),
    (8, -1.606, 0.008869),
    (9, -1.147, 0.00843),
    (10, -1.856, 0.01083),
    (11, -1.735, 0.009302),
    (12, -1.739, 0.008765),
    (13, -1.693, 0.008661),
    (14, -1.299, 0.008634),
    (15, -4.5605, 0.0261),
    (16, -5.6694, 0.0242),
    (17, -6.784, 0.0286),
    (18, -6.0429, 0.0335),
    (19, 0, 0),
    (20, 0, 0),
    (21, 0, 0),
    (22, 0, 0),
    (23, 0, 0),
    (24, 0, 0),
    (25, 0, 0),
]

# Visible band coefficients (bands 0-6) for each calibration set
_VIS_COEFFS: Dict[str, List[Tuple[int, float, float]]] = {
    "business": [
        (0, -2.8897673760346585, 0.026092955565603558),
        (1, -3.4571798300485845, 0.025337091635494854),
        (2, -6.70530755496166, 0.027094761631082855),
        (3, -3.694272567991622, 0.027065205356990654),
        (4, -3.0625245866597273, 0.021895992332693938),
        (5, -3.4559566108943396, 0.025646709986914866),
        (6, -2.8142294442892704, 0.020665962472166604),
    ],
    "recali": [
        (0, -2.670518037371008, 0.02686365784076485),
        (1, -3.3857848098742944, 0.025613468737277403),
        (2, -6.601329459108014, 0.027436789171352626),
        (3, -4.03046463353921, 0.02756314504001579),
        (4, -2.7638263960285325, 0.022092171487461485),
        (5, -3.9083634596407166, 0.02581888240234926),
        (6, -3.0910706697307515, 0.021727999435385205),
    ],
    "old": [
        (0, -2.9931, 0.0249),
        (1, -3.253, 0.0254),
        (2, -6.2679, 0.0261),
        (3, -3.9423, 0.0272),
        (4, -3.2619, 0.0219),
        (5, -3.882, 0.0259),
        (6, -2.9822, 0.0213),
    ],
}

# Date-specific recalibration (same format as recali but different dates)
_DATE_COEFFS: Dict[str, List[Tuple[int, float, float]]] = {
    "20200308": [
        (0, -2.66591297485679, 0.0257526576160977),
        (1, -3.217270198338979, 0.025105374805732106),
        (2, -6.752785327333475, 0.02708525459110231),
        (3, -3.696114985920946, 0.026963024593657226),
        (4, -3.040902007535548, 0.021919875482267495),
        (5, -3.6007150796502416, 0.025762141350630992),
        (6, -3.2846206564969105, 0.02066528880040864),
    ],
    "20200715": [
        (0, -2.493780010063756, 0.025503602060310816),
        (1, -2.8710511719712137, 0.024570012691957962),
        (2, -6.070098165720142, 0.02652233234819512),
        (3, -3.7613520566291507, 0.02745375510561476),
        (4, -3.35129350622244, 0.02178564402746276),
        (5, -3.3606389047805183, 0.025427156834131867),
        (6, -2.856404403626115, 0.021348302374734136),
    ],
}

# Merge date-specific into main lookup
CALIBRATIONS = dict(_VIS_COEFFS)
CALIBRATIONS.update(_DATE_COEFFS)


def build_coeff_list(cal_name: str) -> List[Tuple[int, float, float]]:
    """Build full 26-band coefficient list for a given calibration name."""
    if cal_name not in CALIBRATIONS:
        raise KeyError(f"Unknown calibration: {cal_name}. Available: {list(CALIBRATIONS.keys())}")
    vis = CALIBRATIONS[cal_name]
    return vis + _COMMON_IR


def generate_xcfg(cal_name: str) -> str:
    """Generate the content of a VIS_Cal_Coeff.xcfg file."""
    coeffs = build_coeff_list(cal_name)
    lines = ["VIS_Cal_Coeff = ["]
    for band, offset, slope in coeffs:
        lines.append(f" {band}, {offset}, {slope},")
    lines.append("];")
    # Remove trailing comma from last data line
    if len(lines) > 2:
        lines[-2] = lines[-2].rstrip(",")
    return "\n".join(lines) + "\n"


def write_cal_xcfg(cal_name: str, output_path: str = "VIS_Cal_Coeff.xcfg") -> None:
    """Write the calibration coefficient file for Fortran to read."""
    content = generate_xcfg(cal_name)
    with open(output_path, "w") as f:
        f.write(content)


def setup_cal_mode(cal_name: str, code_root: str) -> None:
    """Configure calibration mode for the Fortran runtime.

    Writes cal_mode.txt (containing 'recali' for external cal, empty for business)
    and VIS_Cal_Coeff.xcfg with the appropriate coefficients.
    """
    cal_mode_path = os.path.join(code_root, "cal_mode.txt")
    xcfg_path = os.path.join(code_root, "VIS_Cal_Coeff.xcfg")

    if cal_name == "business":
        # Business: use HDF5 built-in coeffs, remove cal_mode.txt
        if os.path.exists(cal_mode_path):
            os.remove(cal_mode_path)
        # Still write xcfg for reference (not used by Fortran when cal_mode is absent)
    else:
        # External calibration: write cal_mode.txt and xcfg
        with open(cal_mode_path, "w") as f:
            f.write("recali\n")
        write_cal_xcfg(cal_name, xcfg_path)
