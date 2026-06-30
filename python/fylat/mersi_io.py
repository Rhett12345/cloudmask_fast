"""MERSI-II L1 HDF5 IO module — Python replacement for Fortran io_module.f90.

Reads FY-3D MERSI-II GEO and L1B HDF5 data with full calibration
(DN -> radiance/reflectance -> brightness temperature).

Usage:
    from fylat.mersi_io import read_geo, read_l1b, MersiL1Reader

    reader = MersiL1Reader(l1b_path, geo_path, recalibration_dir=None)
    data = reader.read_all()
    print(data['bt_108'].shape)  # (2000, 2048) BT at 10.8um
"""

import os
import glob
import importlib
import sys
from typing import Dict, Optional

import h5py
import numpy as np

# --- Physical constants ------------------------------------------------------
C1 = 1.191042e-5   # mW / (m^2 sr cm^-4)
C2 = 1.4387752      # K cm

# --- FY-3D MERSI-II channel parameters (from platform_module.f90) -----------
# Central wavelengths (um) and wavenumbers (cm^-1) for IR bands 20-25
IR_WAVELENGTH = {
    20: 3.79599, 21: 4.04587, 22: 7.23264,
    23: 8.56031, 24: 10.7139, 25: 11.94827,
}
IR_WAVENUMBER = {
    20: 2643.4359, 21: 2471.654, 22: 1382.621,
    23: 1168.182, 24: 933.364, 25: 836.941,
}

# SRF correction coefficients (tci, tcs) for FY-3D (from planck_module.f90)
SRF_TCI = {20: 0.5072, 21: 0.3493, 22: 0.4093, 23: 0.1014, 24: 0.5763, 25: 0.4317}
SRF_TCS = {20: 0.9992917440, 21: 0.9994814177, 22: 0.9989956900,
           23: 0.9997135336, 24: 0.9980397975, 25: 0.9983777125}

# IR BT valid range for clipping (K)
BT_VALID_RANGE = {
    20: (180.0, 360.0), 21: (180.0, 360.0), 22: (180.0, 340.0),
    23: (180.0, 340.0), 24: (160.0, 340.0), 25: (160.0, 340.0),
}

_CPP_BACKEND = None
_CPP_BACKEND_LOADED = False


def _project_root() -> str:
    return os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def _load_cpp_backend():
    """Load the optional C++ HDF5 backend built by CMake.

    The source tree workflow usually leaves the extension in build directories
    such as build_migration/.  We discover those locations automatically so the
    Python IO layer can use the migrated C++ reader without an install step.
    """
    global _CPP_BACKEND, _CPP_BACKEND_LOADED
    if _CPP_BACKEND_LOADED:
        return _CPP_BACKEND

    _CPP_BACKEND_LOADED = True
    try:
        _CPP_BACKEND = importlib.import_module("fylat_py")
        return _CPP_BACKEND
    except ImportError:
        pass

    root = _project_root()
    candidates = sorted(glob.glob(os.path.join(root, "build*", "fylat_py*.so")))
    candidates += sorted(glob.glob(os.path.join(root, "build*", "python", "fylat_py*.so")))
    for candidate in candidates:
        module_dir = os.path.dirname(candidate)
        if module_dir not in sys.path:
            sys.path.insert(0, module_dir)
        try:
            _CPP_BACKEND = importlib.import_module("fylat_py")
            return _CPP_BACKEND
        except ImportError:
            continue

    _CPP_BACKEND = None
    return None


def io_backend_name() -> str:
    """Return the active low-level dataset IO backend name."""
    mode = os.environ.get("FYLAT_IO_BACKEND", "auto").lower()
    if mode == "h5py":
        return "h5py"
    return "cpp" if _load_cpp_backend() is not None else "h5py"


def _read_dataset(file_path: str, dataset_path: str, dtype=np.float64) -> np.ndarray:
    """Read a numeric HDF5 dataset through C++ when available.

    Attributes still use h5py because the current C++ phase-1 backend is focused
    on moving large dataset payload reads first.
    """
    mode = os.environ.get("FYLAT_IO_BACKEND", "auto").lower()
    if mode not in {"auto", "cpp", "h5py"}:
        raise ValueError(f"Unsupported FYLAT_IO_BACKEND: {mode}")

    if mode != "h5py":
        backend = _load_cpp_backend()
        if backend is not None:
            try:
                return np.asarray(backend.read_float32(file_path, dataset_path), dtype=dtype)
            except Exception:
                if mode == "cpp":
                    raise

    with h5py.File(file_path, "r") as f:
        return f[dataset_path][:].astype(dtype)


def _read_cpp_payload(file_path: str, method_name: str) -> Optional[Dict[str, np.ndarray]]:
    """Read a business-level C++ IO payload when the backend provides it."""
    mode = os.environ.get("FYLAT_IO_BACKEND", "auto").lower()
    if mode == "h5py":
        return None

    backend = _load_cpp_backend()
    if backend is None or not hasattr(backend, method_name):
        if mode == "cpp":
            raise RuntimeError(f"C++ IO backend does not provide {method_name}()")
        return None

    try:
        payload = getattr(backend, method_name)(file_path)
        return {key: np.asarray(value) for key, value in payload.items()}
    except Exception:
        if mode == "cpp":
            raise
    return None


def _payload_or_dataset(
    payload: Optional[Dict[str, np.ndarray]],
    key: str,
    file_path: str,
    dataset_path: str,
    dtype=np.float64,
) -> np.ndarray:
    if payload is not None and key in payload:
        return np.asarray(payload[key], dtype=dtype)
    return _read_dataset(file_path, dataset_path, dtype)


def _planck_rad2tbb(rad: np.ndarray, band: int) -> np.ndarray:
    """Convert Planck radiance (mW/(m^2 sr cm^-1)) to brightness temperature (K).

    Uses the wavenumber-domain inverse Planck formula with SRF correction.

    Args:
        rad: Radiance array in mW/(m^2 sr cm^-1).
        band: Instrument band number (20-25).

    Returns:
        Brightness temperature array in K.
    """
    cwn = IR_WAVENUMBER[band]
    tci = SRF_TCI[band]
    tcs = SRF_TCS[band]

    # brite_m (Fortran equivalent with pre-computed constants):
    # Fortran: C2_SI * (100*cwn) / ln(C1_SI * (100*cwn)^3 / (1e-5*R) + 1)
    # Simplifies to: C2 * cwn / ln(C1 * cwn^3 / R + 1)
    # where C1 = C1_SI * 1e11, C2 = C2_SI * 100
    rad_clipped = np.maximum(rad, 1e-10)
    bt_planck = C2 * cwn / np.log(C1 * cwn**3 / rad_clipped + 1.0)

    # SRF correction: BT_final = (BT_planck - tci) / tcs
    bt = (bt_planck - tci) / tcs

    # Safety clip
    bt_min, bt_max = BT_VALID_RANGE[band]
    bt = np.clip(bt, bt_min, bt_max)
    return bt.astype(np.float32)


def read_geo(geo_path: str) -> Dict[str, np.ndarray]:
    """Read FY-3D MERSI-II GEO HDF5 file.

    Returns dict with keys:
        lat, lon, sza, saa, vza, vaa, dem, lsm, rel_azimuth
    """
    with h5py.File(geo_path, "r") as f:
        geo = {}
        payload = _read_cpp_payload(geo_path, "read_mersi_geo")

        lat = _payload_or_dataset(payload, "lat", geo_path, "Geolocation/Latitude", np.float64)
        lon = _payload_or_dataset(payload, "lon", geo_path, "Geolocation/Longitude", np.float64)
        geo["lat"] = lat
        geo["lon"] = lon

        # Solar/Sensor angles: stored as scaled int16, real = (DN + Intercept) * Slope
        for key, hdf_key in [
            ("sza", "SolarZenith"), ("saa", "SolarAzimuth"),
            ("vza", "SensorZenith"), ("vaa", "SensorAzimuth"),
        ]:
            ds = f[f"Geolocation/{hdf_key}"]
            slope = ds.attrs.get("Slope", 1.0)
            intercept = ds.attrs.get("Intercept", 0.0)
            fill_value = ds.attrs.get("FillValue", -32767)
            raw = _payload_or_dataset(
                payload, f"{key}_raw", geo_path, f"Geolocation/{hdf_key}", np.float64)
            data = np.where(raw != fill_value, (raw + intercept) * slope, np.nan)
            geo[key] = data.astype(np.float32)

        # DEM
        geo["dem"] = _payload_or_dataset(payload, "dem", geo_path, "Geolocation/DEM", np.float32)

        # LandSeaMask (uint8 -> int32)
        geo["lsm"] = _payload_or_dataset(
            payload, "lsm", geo_path, "Geolocation/LandSeaMask", np.int32)

        # Compute relative azimuth: for FY-3D, rel_az = |180 - wrapped_diff|
        saa = geo["saa"]
        vaa = geo["vaa"]
        diff = np.abs(saa - vaa)
        mdphi = np.where(diff > 180.0, 360.0 - diff, diff)
        geo["rel_azimuth"] = np.abs(180.0 - mdphi).astype(np.float32)

        geo["nlines"] = lat.shape[0]
        geo["npixels"] = lat.shape[1]

    return geo


def read_l1b(
    l1b_path: str,
    recal_xcfg_path: Optional[str] = None,
) -> Dict[str, np.ndarray]:
    """Read FY-3D MERSI-II L1B HDF5 file with full calibration.

    Args:
        l1b_path: Path to L1B HDF5 file.
        recal_xcfg_path: Optional path to VIS_Cal_Coeff.xcfg for recalibration.
                          If None, uses built-in HDF5 calibration coefficients.

    Returns dict with keys:
        ref_vis: (nlines, npixels, 19) reflectance for bands 1-19
        rad_ir: (nlines, npixels, 6) radiance for bands 20-25 (mW/(m^2 sr cm^-1))
        bt_ir: (nlines, npixels, 6) brightness temperature for bands 20-25 (K)
        vis_cal_coef: (3, 19) VIS calibration coefficients
        ir_cal_coef: (3, 6) IR calibration coefficients
    """
    with h5py.File(l1b_path, "r") as f:
        payload = _read_cpp_payload(l1b_path, "read_mersi_l1_payload")

        # --- Read calibration coefficients ---
        # VIS_Cal_Coeff: (19, 3) for bands 1-19 (bands first, coeffs second)
        vis_cal_raw = _payload_or_dataset(
            payload, "vis_cal_coeff", l1b_path, "Calibration/VIS_Cal_Coeff", np.float64)
        if vis_cal_raw.shape[0] == 19:  # (19, 3) -> (3, 19)
            vis_cal = vis_cal_raw.T
        else:
            vis_cal = vis_cal_raw

        # IR_Cal_Coeff: (6, 4, 200) for bands 20-25 (bands, coeffs, scans)
        ir_cal_full = _payload_or_dataset(
            payload, "ir_cal_coeff", l1b_path, "Calibration/IR_Cal_Coeff", np.float64)
        ir_cal = np.zeros((6, 3), dtype=np.float64)
        for b in range(6):
            ir_cal[b, 0:3] = ir_cal_full[b, 0:3, 100]  # scan line 100

        # Apply recalibration override if provided
        if recal_xcfg_path and os.path.exists(recal_xcfg_path):
            vis_cal = _load_recal_xcfg(recal_xcfg_path, vis_cal)

        # --- Read DN data ---
        # VIS bands 1-4: EV_250_Aggr.1KM_RefSB (4, nlines, npixels)
        vis_250 = _payload_or_dataset(
            payload, "ev_250_refsb", l1b_path, "Data/EV_250_Aggr.1KM_RefSB", np.float64)
        # VIS bands 5-19: EV_1KM_RefSB (15, nlines, npixels)
        vis_1km = _payload_or_dataset(
            payload, "ev_1km_refsb", l1b_path, "Data/EV_1KM_RefSB", np.float64)

        # IR bands 20-23: EV_1KM_Emissive (4, nlines, npixels)
        ir_1km = _payload_or_dataset(
            payload, "ev_1km_emissive", l1b_path, "Data/EV_1KM_Emissive", np.float64)
        # IR bands 24-25: EV_250_Aggr.1KM_Emissive (2, nlines, npixels)
        ir_250 = _payload_or_dataset(
            payload, "ev_250_emissive", l1b_path, "Data/EV_250_Aggr.1KM_Emissive", np.float64)

        # Read Slope/Intercept for IR bands (FY-3D applies scaling before DN->rad)
        ir_1km_slope = f["Data/EV_1KM_Emissive"].attrs.get("Slope", np.ones(4))
        ir_1km_intercept = f["Data/EV_1KM_Emissive"].attrs.get("Intercept", np.zeros(4))
        ir_250_slope = f["Data/EV_250_Aggr.1KM_Emissive"].attrs.get("Slope", np.ones(2))
        ir_250_intercept = f["Data/EV_250_Aggr.1KM_Emissive"].attrs.get("Intercept", np.zeros(2))

    nlines, npixels = vis_1km.shape[1], vis_1km.shape[2]

    # --- DN to reflectance for VIS bands ---
    # Formula for FY-3D: ref = (c0 + c1*DN + c2*DN^2) * 0.01 / cos(SZA)
    # Without SZA correction here (cos(SZA) applied separately with GEO data)
    ref_vis = np.zeros((nlines, npixels, 19), dtype=np.float32)
    for b in range(19):
        c0, c1, c2 = vis_cal[0, b], vis_cal[1, b], vis_cal[2, b]
        if b < 4:
            dn = vis_250[b, :, :]
        else:
            dn = vis_1km[b - 4, :, :]
        ref = (c0 + c1 * dn + c2 * dn**2) * 0.01
        ref_vis[:, :, b] = ref.astype(np.float32)

    # --- DN to radiance for IR bands ---
    rad_ir = np.zeros((nlines, npixels, 6), dtype=np.float32)
    bt_ir = np.zeros((nlines, npixels, 6), dtype=np.float32)

    for b_local in range(6):
        band = b_local + 20

        if b_local < 4:
            dn = ir_1km[b_local, :, :]
            slope = ir_1km_slope[b_local]
            intercept = ir_1km_intercept[b_local]
        else:
            dn = ir_250[b_local - 4, :, :]
            slope = ir_250_slope[b_local - 4]
            intercept = ir_250_intercept[b_local - 4]

        # Apply Slope/Intercept scaling (FY-3D specific)
        dn_scaled = (dn + intercept) * slope

        # For FY-3D (sensor_id=21), radiance = DN (passthrough, already mW/(m^2 sr cm^-1))
        rad = dn_scaled.astype(np.float64)

        # Convert to brightness temperature
        bt = _planck_rad2tbb(rad, band)

        rad_ir[:, :, b_local] = rad.astype(np.float32)
        bt_ir[:, :, b_local] = bt.astype(np.float32)

    return {
        "ref_vis": ref_vis,
        "rad_ir": rad_ir,
        "bt_ir": bt_ir,
        "vis_cal_coef": vis_cal,
        "ir_cal_coef": ir_cal,
    }


def _load_recal_xcfg(
    xcfg_path: str,
    default_vis_cal: np.ndarray,
) -> np.ndarray:
    """Load VIS recalibration coefficients from .xcfg file.

    Format:
        header_line
        IC, coef0, coef1  (19 bands, 1 per line)

    Returns updated vis_cal array (3, 19).
    """
    vis_cal = default_vis_cal.copy()
    with open(xcfg_path) as f:
        lines = f.readlines()

    band_idx = 0
    for line in lines:
        line = line.strip()
        if not line or line.startswith("#") or line.startswith("Intercept"):
            continue
        parts = line.replace(",", " ").split()
        if len(parts) >= 3:
            try:
                ic = float(parts[0])
                coef0 = float(parts[1])
                coef1 = float(parts[2])
                if band_idx < 19:
                    vis_cal[0, band_idx] = coef0  # c0
                    vis_cal[1, band_idx] = coef1  # c1
                    vis_cal[2, band_idx] = 0.0     # c2 not in .xcfg
                    band_idx += 1
            except ValueError:
                continue

    return vis_cal


class MersiL1Reader:
    """High-level reader for FY-3D MERSI-II L1 data.

    Combines GEO and L1B reading with SZA correction for VIS reflectance.

    Usage:
        reader = MersiL1Reader(l1b_path, geo_path)
        data = reader.read_all()
        # data['bt_108'] = BT at 10.8um (ch24), shape (2000, 2048)
        # data['ref_138'] = reflectance at 1.38um (ch19), shape (2000, 2048)
    """

    def __init__(
        self,
        l1b_path: str,
        geo_path: str,
        recal_xcfg_path: Optional[str] = None,
    ):
        self.l1b_path = l1b_path
        self.geo_path = geo_path
        self.recal_xcfg_path = recal_xcfg_path

    def read_all(self) -> Dict[str, np.ndarray]:
        """Read and calibrate all data.

        Returns dict with named channel arrays for convenience:
            ref_{band} (1-19): VIS reflectance
            bt_{band} (20-25): IR brightness temperature (K)
            rad_{band} (20-25): IR radiance
            lat, lon, sza, vza, lsm, dem: geolocation
        """
        geo = read_geo(self.geo_path)
        l1b = read_l1b(self.l1b_path, self.recal_xcfg_path)

        # Apply SZA correction to VIS reflectance
        sza_rad = np.deg2rad(geo["sza"])
        cos_sza = np.cos(sza_rad)
        cos_sza = np.maximum(cos_sza, 0.01)  # avoid division by zero

        result = {}
        for key in ["lat", "lon", "sza", "vza", "lsm", "dem", "rel_azimuth"]:
            result[key] = geo[key]

        for b in range(19):
            ref = l1b["ref_vis"][:, :, b] / cos_sza
            ref = np.clip(ref, 0.0, 1.5)
            result[f"ref_{b + 1}"] = ref.astype(np.float32)

        for b_local in range(6):
            band = b_local + 20
            result[f"rad_{band}"] = l1b["rad_ir"][:, :, b_local]
            result[f"bt_{band}"] = l1b["bt_ir"][:, :, b_local]

        result["nlines"] = geo["nlines"]
        result["npixels"] = geo["npixels"]

        return result
