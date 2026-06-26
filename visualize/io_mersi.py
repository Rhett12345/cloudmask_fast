"""
io_mersi.py — FY-3D MERSI-II L1B + CLM data readers
=====================================================
Clean, modular version of the I/O functions from visualize_clm_nature.py.

Public API
----------
load_mersi_l1b(l1b_path)
    → dict with keys: rgb (H,W,3), lat, lon

load_clm_hdf5(path)
    → dict with keys: clm (H,W), lat, lon

parse_mersi_datetime(filename)
    → datetime (UTC-aware)

find_l1b_for_clm(clm_path, mersi_root)
    → path string or None
"""

from __future__ import annotations
import os
import re
from pathlib import Path
from datetime import datetime, timezone

import numpy as np
import h5py


# ─────────────────────────────────────────────────────────────────────────────
# 1.  Filename utilities
# ─────────────────────────────────────────────────────────────────────────────

def parse_mersi_datetime(filename: str) -> datetime | None:
    """
    Extract UTC datetime from FY-3D MERSI filename.
    Pattern: ..._YYYYMMDD_HHMM_...
    Returns UTC-aware datetime or None.
    """
    m = re.search(r'(\d{8})_(\d{4})', os.path.basename(filename))
    if not m:
        return None
    d, t = m.group(1), m.group(2)
    try:
        return datetime(int(d[:4]), int(d[4:6]), int(d[6:8]),
                        int(t[:2]),  int(t[2:]),
                        tzinfo=timezone.utc)
    except ValueError:
        return None


def find_l1b_for_clm(
    clm_path: str,
    mersi_root: str | Path = "/data/Data_yuq/mersi",
) -> str | None:
    """
    Locate the MERSI L1B 1-km HDF file corresponding to a CLM file.
    Searches <mersi_root>/<YYYYMMDD>/FY3D_MERSI_GBAL_L1_<date>_<time>_1000M_MS.HDF
    """
    m = re.search(r'(\d{8})_(\d{4})', os.path.basename(clm_path))
    if not m:
        return None
    date_str, time_tag = m.group(1), m.group(2)
    p = (Path(mersi_root) / date_str /
         f"FY3D_MERSI_GBAL_L1_{date_str}_{time_tag}_1000M_MS.HDF")
    return str(p) if p.exists() else None


# ─────────────────────────────────────────────────────────────────────────────
# 2.  MERSI L1B reader  →  RGB + geolocation
# ─────────────────────────────────────────────────────────────────────────────

def load_mersi_l1b(l1b_path: str) -> dict | None:
    """
    Read MERSI-II 1-km L1B HDF and produce a display-ready RGB array.

    Uses bands B1/B2/B3 (EV_250_Aggr.1KM_RefSB channels 0,1,2)
    and the VIS calibration coefficients.

    Returns
    -------
    dict with:
      'rgb'  : (H, W, 3) float32, values clipped to [0, 1]
      'lat'  : (H, W) float64
      'lon'  : (H, W) float64
      'path' : input file path
    or None on failure.
    """
    if not l1b_path or not os.path.exists(l1b_path):
        print(f"[WARN] L1B not found: {l1b_path}")
        return None

    try:
        with h5py.File(l1b_path, "r") as f:
            vis_250 = f["Data/EV_250_Aggr.1KM_RefSB"][:].astype(np.float64)
            vis_cal = f["Calibration/VIS_Cal_Coeff"][:]
            esd     = float(np.squeeze(f.attrs.get("EarthSun Distance Ratio", 1.0)))

            # Geolocation — try several common paths
            lat, lon = _read_l1b_geo(f)

        esd2   = esd ** 2
        n_bands, n_line, n_elem = vis_250.shape
        rgb = np.zeros((n_line, n_elem, 3), dtype=np.float32)

        for i in range(3):                        # R, G, B → bands 0, 1, 2
            c0, c1, c2 = vis_cal[i]
            dn   = vis_250[i]
            refl = (c0 + c1 * dn + c2 * dn * dn) * 0.01 / esd2
            rgb[:, :, i] = refl.astype(np.float32)

        rgb = _stretch_rgb(rgb)

        return {"rgb": rgb, "lat": lat, "lon": lon, "path": l1b_path}

    except Exception as e:
        print(f"[ERROR] Loading L1B {l1b_path}: {e}")
        return None


def _read_l1b_geo(f: h5py.File) -> tuple[np.ndarray, np.ndarray]:
    """Try multiple common geolocation paths inside L1B HDF."""
    # Ordered by likelihood
    paths = [
        ("Geolocation/Latitude",    "Geolocation/Longitude"),
        ("Data/Latitude",           "Data/Longitude"),
        ("Navigation/Latitude",     "Navigation/Longitude"),
    ]
    for lp, lonp in paths:
        if lp in f and lonp in f:
            lat = f[lp][:].astype(np.float64)
            lon = f[lonp][:].astype(np.float64)
            lat = np.where((lat < -90)  | (lat > 90),   np.nan, lat)
            lon = np.where((lon < -180) | (lon > 180),  np.nan, lon)
            return lat, lon
    raise KeyError("Geolocation datasets not found in L1B HDF")


def _stretch_rgb(rgb: np.ndarray, lo_pct: float = 2, hi_pct: float = 98) -> np.ndarray:
    """Per-channel percentile stretch → [0, 1]."""
    out = np.zeros_like(rgb)
    for ch in range(3):
        band  = rgb[:, :, ch]
        valid = band > 0
        if not valid.any():
            continue
        p2, p98 = np.percentile(band[valid], [lo_pct, hi_pct])
        if p98 <= p2:
            continue
        out[:, :, ch] = np.clip((band - p2) / (p98 - p2 + 1e-10), 0, 1)
    return out


# ─────────────────────────────────────────────────────────────────────────────
# 3.  MERSI CLM reader  (HDF5)
# ─────────────────────────────────────────────────────────────────────────────

def load_clm_hdf5(path: str, geo_root: str = "/data/Data_yuq/mersi") -> dict | None:
    """
    Read MERSI CLM from HDF5.  Handles three storage conventions:

      Flat:   /cm, /conf, /lat, /lon  (root-level datasets)
      Group:  Cloud_Mask_1km/  (grouped datasets)
      Raw:    /Cloud_Mask (N,H,W) raw bitmask, /Quality_Assurance (M,H,W)
              — lat/lon loaded from external GEO HDF5 file.

    Returns
    -------
    dict with keys: clm (H,W int32), lat, lon, path
    or None on failure.
    """
    if not os.path.exists(path):
        print(f"[ERROR] CLM file not found: {path}")
        return None
    try:
        with h5py.File(path, "r") as f:
            lat, lon, clm = None, None, None

            # ── Convention A: flat cm/conf/lat/lon ──
            if "cm" in f:
                lat = f["lat"][:].astype(np.float64)
                lon = f["lon"][:].astype(np.float64)
                clm = f["cm"][:].astype(np.int32)
                if np.all(clm == 0) and "conf" in f:
                    derived = _clm_from_confidence(f["conf"][:].astype(np.float64))
                    if np.any(derived >= 0):
                        clm = derived

            # ── Convention B: Cloud_Mask_1km group ──
            elif "Cloud_Mask_1km" in f:
                grp = f["Cloud_Mask_1km"]
                lat = grp["Latitude"][:].astype(np.float64)
                lon = grp["Longitude"][:].astype(np.float64)
                clm = grp["Cloud_Mask_Value"][:].astype(np.int32)
                if np.all(clm == 0) and "Cloud_Mask" in grp:
                    decoded = _decode_bitmask(grp["Cloud_Mask"][:])
                    if np.any(decoded >= 0):
                        clm = decoded
                if np.all(clm == 0) and "Confidence" in grp:
                    derived = _clm_from_confidence(grp["Confidence"][:])
                    if np.any(derived >= 0):
                        clm = derived

            # ── Convention C: raw /Cloud_Mask bitmask (N,H,W) ──
            elif "Cloud_Mask" in f:
                cm_raw = f["Cloud_Mask"][:]
                clm = _decode_bitmask(cm_raw)
                # Try to read embedded geolocation; if absent, load from GEO file
                lat, lon = _read_embedded_geo(f)
                if lat is None:
                    geo_path = find_geo_for_clm(path, geo_root)
                    if geo_path:
                        lat, lon = _read_geo_from_file(geo_path)
                    else:
                        print(f"[ERROR] No GEO file found for {path}")
                        return None

            else:
                print(f"[ERROR] Unknown HDF5 structure in {path}")
                return None

        if lat is None or lon is None or clm is None:
            print(f"[ERROR] Incomplete data in {path}")
            return None

        lat = np.where((lat < -90)  | (lat > 90),   np.nan, lat)
        lon = np.where((lon < -180) | (lon > 180),  np.nan, lon)
        clm = np.where((clm < 0)   | (clm > 3),     -1,    clm)

        return {"clm": clm, "lat": lat, "lon": lon, "path": path}

    except Exception as e:
        print(f"[ERROR] Loading CLM {path}: {e}")
        return None


def _read_embedded_geo(f: h5py.File) -> tuple[np.ndarray | None, np.ndarray | None]:
    """Try to find lat/lon inside a CLM HDF5."""
    for lp, lonp in [
        ("lat", "lon"),
        ("Latitude", "Longitude"),
        ("Geolocation/Latitude", "Geolocation/Longitude"),
    ]:
        if lp in f and lonp in f:
            return (f[lp][:].astype(np.float64),
                    f[lonp][:].astype(np.float64))
    return None, None


def _read_geo_from_file(geo_path: str) -> tuple[np.ndarray | None, np.ndarray | None]:
    """Read lat/lon from a MERSI GEO HDF5 file."""
    try:
        with h5py.File(geo_path, "r") as f:
            for lp, lonp in [
                ("Geolocation/Latitude", "Geolocation/Longitude"),
                ("Latitude", "Longitude"),
            ]:
                if lp in f and lonp in f:
                    return (f[lp][:].astype(np.float64),
                            f[lonp][:].astype(np.float64))
    except Exception as e:
        print(f"[ERROR] Reading GEO {geo_path}: {e}")
    return None, None


def find_geo_for_clm(
    clm_path: str,
    mersi_root: str | Path = "/data/Data_yuq/mersi",
) -> str | None:
    """Locate the MERSI GEO 1-km HDF file corresponding to a CLM file."""
    m = re.search(r'(\d{8})_(\d{4})', os.path.basename(clm_path))
    if not m:
        return None
    date_str, time_tag = m.group(1), m.group(2)
    p = (Path(mersi_root) / date_str /
         f"FY3D_MERSI_GBAL_L1_{date_str}_{time_tag}_GEO1K_MS.HDF")
    return str(p) if p.exists() else None

def _decode_bitmask(cm_raw: np.ndarray) -> np.ndarray:
    """Decode multi-byte CLM bitmask array → class 0-3.

    Handles two storage conventions:
      • (H, W, N)  — last axis is byte index
      • (N, H, W)  — first axis is byte index

    Fortran ibset convention (LSB-based, matches HDF5 byte storage and
    Fortran convert_cloud_mask which uses ibits(...,0,1), ibits(...,1,1), ibits(...,2,1)):
      bit 0 (value 1):   Cloud Mask Flag   (0=not determined, 1=determined)
      bits 1-2 (values 2,4): FOV Quality   (00=Cloudy, 01=Uncertain,
                                              10=Probably Clear, 11=Confident Clear)
    """
    if cm_raw.ndim == 3 and cm_raw.shape[0] < cm_raw.shape[-1]:
        byte0 = cm_raw[0].astype(np.uint8)           # (N, H, W) → byte 0
    elif cm_raw.ndim == 3:
        byte0 = cm_raw[:, :, 0].astype(np.uint8)     # (H, W, N) → byte 0
    else:
        byte0 = cm_raw.astype(np.uint8)

    determined = (byte0 & 1).astype(bool)              # bit 0 (LSB) = processed flag
    cat_bits   = (byte0 >> 1) & 3                     # bits 1-2 = cloud category

    result = np.full(byte0.shape, -1, dtype=np.int32)
    result[determined] = cat_bits[determined]          # 0=cloudy, 1=uncertain,
                                                       # 2=prob_clear, 3=conf_clear
    return result


def _clm_from_confidence(conf: np.ndarray) -> np.ndarray:
    """Bin float confidence [0,1] into 4 CLM classes."""
    clm   = np.full(conf.shape, -1, dtype=np.int32)
    valid = np.isfinite(conf)
    clm[valid & (conf < 0.33)]               = 0
    clm[valid & (conf >= 0.33) & (conf < 0.66)] = 1
    clm[valid & (conf >= 0.66) & (conf < 0.90)] = 2
    clm[valid & (conf >= 0.90)]                  = 3
    return clm


# ─────────────────────────────────────────────────────────────────────────────
# 4.  Distribution printer  (moved here from visualize_clm_nature.py)
# ─────────────────────────────────────────────────────────────────────────────

CLM_LABEL = {0: "Cloudy", 1: "Prob. Cloudy", 2: "Prob. Clear", 3: "Conf. Clear"}


def print_clm_distribution(label: str, clm: np.ndarray) -> None:
    """Print pixel category distribution to stdout."""
    valid_mask  = (clm >= 0) & (clm <= 3)
    total_valid = int(valid_mask.sum())
    total_all   = clm.size
    sep = "═" * 62
    print(sep)
    print(f"  CLM Distribution — {label}")
    print(f"  Total pixels : {total_all:>10,}  |  Valid : {total_valid:>10,}")
    print("─" * 62)
    print(f"  {'Category':<20} {'Class':>5}  {'Count':>10}  {'Ratio':>7}")
    print("─" * 62)
    for v in range(4):
        cnt = int((clm == v).sum())
        pct = 100.0 * cnt / total_valid if total_valid > 0 else 0.0
        bar = "█" * int(pct / 2)
        print(f"  {CLM_LABEL[v]:<20} {v:>5}  {cnt:>10,}  {pct:>6.2f}%  {bar}")
    invalid = int((clm < 0).sum())
    inv_pct = 100.0 * invalid / total_all if total_all > 0 else 0.0
    print("─" * 62)
    print(f"  {'Invalid/unprocessed':<20} {'–':>5}  {invalid:>10,}  {inv_pct:>6.2f}%")
    print(sep + "\n")


# ─────────────────────────────────────────────────────────────────────────────
# 5.  MERSI-II thermal IR reader — 10.8 µm brightness temperature
# ─────────────────────────────────────────────────────────────────────────────

# MERSI-II 1-km thermal emissive bands stored in EV_1KM_Emissive (H5).
# Actual FY-3D MERSI-II L1B layout — EV_1KM_Emissive has 4 bands:
#   Idx  0 : 8.55 µm  (band 24)
#   Idx  1 : 10.8 µm  (band 25)  ← standard cloud / IR window
#   Idx  2 : 12.0 µm  (band 26)
#   Idx  3 : 13.5 µm  (band 27)
# EV_250_Aggr.1KM_Emissive has 2 more (bands 20, 21/22: 3.8, 4.05 µm).
#
# Calibration: DN → radiance → BT (inverse Planck)
# IR_Cal_Coeff shape (6, 4, 200): 6 IR bands × 4 coeffs × 200 scan-lines.
# EV_1KM_Emissive band i maps to IR_Cal_Coeff band (i+2).
# Calibration formula (quadratic, per-wavenumber):
#   L_wn  = c0 + c1*DN + c2*DN²          [mW m⁻² sr⁻¹ cm⁻¹]
# Convert per-wavenumber → per-wavelength for Planck:
#   L_wl  = L_wn / λ²                    [mW m⁻² sr⁻¹ cm]
# BT (K) T = c2 / (λ · ln(c1 / (λ⁵ L_wl) + 1))

_C1 = 1.1910427e-5   # mW m⁻² sr⁻¹ cm⁴
_C2 = 1.4387752      # cm·K

# EV_1KM_Emissive band index → (central wavelength cm, IR_Cal_Coeff index)
_EV1KM_BAND_META = {
    # ev_idx: (lambda_cm, cal_idx)
    0: (8.55e-4,  2),   # band 24: 8.55 µm
    1: (10.80e-4, 3),   # band 25: 10.8 µm  ← IR window
    2: (12.00e-4, 4),   # band 26: 12.0 µm
    3: (13.50e-4, 5),   # band 27: 13.5 µm
}

_BT108_EV_IDX = 1   # index of 10.8 µm inside EV_1KM_Emissive


def load_mersi_bt108(
    l1b_path: str,
    ev_band_idx: int = _BT108_EV_IDX,
) -> np.ndarray | None:
    """
    Read MERSI-II 1-km L1B HDF and return 10.8 µm brightness temperature.

    Calibration pipeline:
      1. Read raw DN from EV_1KM_Emissive[ev_band_idx]
      2. Apply quadratic calibration: L_wn = c0 + c1*DN + c2*DN²
         (coefficients from IR_Cal_Coeff, averaged over scan-lines)
      3. Convert per-wavenumber → per-wavelength: L_wl = L_wn / λ²
      4. Convert radiance → BT via inverse Planck

    Parameters
    ----------
    l1b_path : str
        Path to FY3D_MERSI_GBAL_L1_..._1000M_MS.HDF (HDF5).
    ev_band_idx : int
        Index of the desired emissive band inside EV_1KM_Emissive.
        Default 1 = 10.8 µm.

    Returns
    -------
    bt : (H, W) float32 array of brightness temperature in Kelvin.
         Invalid / fill pixels are set to NaN.
    None on any read/calibration failure.
    """
    if not l1b_path or not os.path.exists(l1b_path):
        print(f"[WARN] L1B not found for BT108: {l1b_path}")
        return None

    try:
        with h5py.File(l1b_path, "r") as f:

            # ── Locate the emissive dataset ──────────────────────────
            emissive_paths = [
                "Data/EV_1KM_Emissive",
                "EV_1KM_Emissive",
            ]
            ev_data = None
            for ep in emissive_paths:
                if ep in f:
                    ev_data = f[ep][ev_band_idx].astype(np.float64)
                    break
            if ev_data is None:
                print(f"[WARN] EV_1KM_Emissive not found in {l1b_path}")
                return None

            # ── Calibration coefficients ────────────────────────────
            cal_paths = [
                "Calibration/IR_Cal_Coeff",
                "Calibration/IRB_Cal_Coeff",
                "IR_Cal_Coeff",
                "IRB_Cal_Coeff",
            ]
            cal = None
            for cp in cal_paths:
                if cp in f:
                    cal = f[cp][:]   # shape (6, 4, 200) or (n_bands, 3+)
                    break
            if cal is None:
                print(f"[WARN] IR_Cal_Coeff not found in {l1b_path} — "
                      "BT108 not calibrated.")
                return None

            # ── Get wavelength and calibration index for this EV band ─
            lam, cal_idx = _EV1KM_BAND_META[ev_band_idx]

            # ── Parse calibration format ─────────────────────────────
            if cal.ndim == 3:
                # Per-scan format (n_bands, n_coeff, n_scans): average over scans
                cal_band = cal[cal_idx].mean(axis=1)  # (4,) → (n_coeff,)
            else:
                # Simple format (n_bands, n_coeff)
                cal_band = cal[cal_idx]

            c0 = float(cal_band[0])
            c1 = float(cal_band[1])
            c2 = float(cal_band[2]) if len(cal_band) > 2 else 0.0

            # ── Fill-value / valid range ────────────────────────────
            fill_mask = (ev_data == 0) | (ev_data >= 65535)

        # ── DN → radiance (per-wavenumber) ────────────────────────
        rad_wn = c0 + c1 * ev_data + c2 * ev_data * ev_data   # mW m⁻² sr⁻¹ cm⁻¹
        rad_wn[fill_mask] = np.nan
        rad_wn = np.where(rad_wn > 0, rad_wn, np.nan)

        # ── Per-wavenumber → per-wavelength conversion ───────────
        lam2 = lam ** 2
        rad_wl = rad_wn / lam2                                # mW m⁻² sr⁻¹ cm

        # ── Radiance → brightness temperature (inverse Planck) ──
        lam5 = lam ** 5
        bt = _C2 / (lam * np.log(_C1 / (lam5 * rad_wl) + 1.0))

        # Sanity clip: physical BT range 170–340 K
        bt = np.where((bt > 170.0) & (bt < 340.0), bt, np.nan)

        print(f"[BT108] Loaded from {os.path.basename(l1b_path)}: "
              f"valid={np.isfinite(bt).sum()}/{bt.size}, "
              f"BT range {np.nanmin(bt):.1f}–{np.nanmax(bt):.1f} K")
        return bt.astype(np.float32)

    except Exception as e:
        print(f"[ERROR] Loading BT108 from {l1b_path}: {e}")
        return None