"""Planck radiance ↔ brightness temperature conversion.
Exact reimplementation of Fortran planck_module.f90 for FY3D MERSI-II.
"""

import math
import numpy as np

# Planck constants (from planck_module.f90)
C1 = 1.191042722e-16  # 2*h*c^2  [W·m²]
C2 = 1.4387752e-02     # h*c/k    [K·m]

# FY3D MERSI-II band coefficients (sensor_id=21)
# Bands 20-25: effective central wavenumbers + linear correction
FY3D_COEFFS = {
    # band: (cwn_cm1, tcs_slope, tci_intercept_K)
    20: (2643.4359, 0.9992917440, 0.50718071650),   # 3.8 um
    21: (2471.654,  0.9994814177, 0.3493280160),    # 4.0 um (7.2 um proxy)
    22: (1382.621,  0.9989956900, 0.40925130837),   # 7.2 um (8.6 proxy)
    23: (1168.182,  0.9997135336, 0.1014073981),    # 8.6 um (10.8 proxy)
    24: (933.3640,  0.9980397975, 0.57633464244),   # 10.8 um (~11 um)
    25: (836.9410,  0.9983777125, 0.4317181810),    # 12.0 um
}


def rad2bt_w_um(rad_w_um: float, band: int) -> float:
    """Convert radiance (W/m2/sr/um) to brightness temperature (K).
    Uses the Fortran bright_m function — wavelength-domain Planck."""
    if rad_w_um <= 0 or band < 20 or band > 25:
        return 0.0
    cwn, tcs, tci = FY3D_COEFFS[band]
    w_um = 1.0e4 / cwn          # effective wavelength (um)
    ws = w_um * 1.0e-6           # wavelength (m)
    ws5 = ws ** 5
    bt_raw = C2 / (ws * math.log(C1 / (1.0e6 * rad_w_um * ws5) + 1.0))
    return (bt_raw - tci) / tcs


def rad2bt_array(rad: np.ndarray, band: int) -> np.ndarray:
    """Vectorized radiance-to-BT conversion."""
    if band < 20 or band > 25:
        return np.zeros_like(rad)
    cwn, tcs, tci = FY3D_COEFFS[band]
    w_um = 1.0e4 / cwn
    ws = w_um * 1.0e-6
    ws5 = ws ** 5
    valid = rad > 0
    bt_raw = np.zeros_like(rad)
    bt_raw[valid] = C2 / (ws * np.log(C1 / (1.0e6 * rad[valid] * ws5) + 1.0))
    return (bt_raw - tci) / tcs


def bt2rad_w_um(bt: float, band: int) -> float:
    """Convert BT (K) to radiance (W/m2/sr/um). Inverse Planck."""
    if bt <= 0 or band < 20 or band > 25:
        return 0.0
    cwn = FY3D_COEFFS[band][0]
    w_um = 1.0e4 / cwn
    ws = w_um * 1.0e-6
    ws5 = ws ** 5
    bt_raw = bt * FY3D_COEFFS[band][1] + FY3D_COEFFS[band][2]
    return 1.0e-6 * C1 / (ws5 * (math.exp(C2 / (ws * bt_raw)) - 1.0))
