"""Ancillary data loading — Python replacement for Fortran get_ancil_data_module.f90.
Loads ecosystem (IGBP), NISE snow/ice, and OISST directly from HDF4/HDF5 files.
No Fortran dependency.
"""

import math
import os
from pathlib import Path
from typing import Tuple

import numpy as np
from pyhdf.SD import SD, SDC


def _project_root() -> str:
    return str(Path(__file__).parent.parent.parent)


# =========================================================================
# 1. IGBP Ecosystem (HDF4, 43200x21600, ~1km resolution)
# =========================================================================

class EcosystemReader:
    """Reads IGBP ecosystem classification from HDF4 file."""

    def __init__(self, filepath: str = None):
        if filepath is None:
            filepath = os.path.join(_project_root(), "coeff", "fylat_ecosystem.hdf")
        self._filepath = filepath
        self._data = None
        self.nlat = 21600
        self.nlon = 43200
        self.dlat = 180.0 / 21600
        self.dlon = 360.0 / 43200
        self.first_lat = 89.99583
        self.first_lon = -179.99583
        self.factor = -1  # because first_lat > 0

    def _ensure_loaded(self):
        if self._data is not None:
            return
        sd = SD(self._filepath, SDC.READ)
        sds = sd.select("ecosystem_type")
        self._data = sds.get().astype(np.int8)  # (nlat, nlon)
        sd.end()

    def get_type(self, lon: float, lat: float) -> int:
        """Return IGBP ecosystem type (0-17) for a given (lon, lat)."""
        self._ensure_loaded()
        # Exact formula from get_ancil_data_module.f90
        my = int((lat - self.first_lat + 0.5 * self.dlat * self.factor) /
                 (self.dlat * self.factor) + 1)
        mx = int((lon - self.first_lon + 0.5 * self.dlon) /
                 self.dlon + 1)
        my = max(0, min(self.nlat - 1, my - 1))  # 1-based → 0-based
        mx = max(0, min(self.nlon - 1, mx - 1))
        return int(self._data[my, mx])

    def get_types(self, lons: np.ndarray, lats: np.ndarray) -> np.ndarray:
        """Vectorized lookup for arrays of (lon, lat)."""
        self._ensure_loaded()
        # Broadcast-compatible indexing
        my = ((lats - self.first_lat + 0.5 * self.dlat * self.factor) /
              (self.dlat * self.factor) + 1).astype(int) - 1
        mx = ((lons - self.first_lon + 0.5 * self.dlon) /
              self.dlon + 1).astype(int) - 1
        my = np.clip(my, 0, self.nlat - 1)
        mx = np.clip(mx, 0, self.nlon - 1)
        return self._data[my, mx].astype(np.int8)


# =========================================================================
# 2. NISE Snow/Ice Mask (HDF4, 721x721 EASE-Grid polar projection)
# =========================================================================

class NiseReader:
    """Reads NISE near-real-time snow/ice extent from HDF4 monthly files."""

    # EASE-Grid constants
    RE_KM = 6371.228
    CELL_KM = 25.067525
    COS_PHI1 = 0.866025403
    COLS = 721
    ROWS = 721
    SCALE = 1
    RG = SCALE * RE_KM / CELL_KM
    R0 = (COLS - 1) / 2.0 * SCALE
    S0 = (ROWS - 1) / 2.0 * SCALE

    def __init__(self, month: int, filepath: str = None):
        if filepath is None:
            fname = f"NISE_SSMIF13_EASEGRID_M{month:02d}.HDF"
            filepath = os.path.join(_project_root(), "coeff", "sfc_snow_ice", fname)
        self._filepath = filepath
        self._smk_n = None
        self._smk_s = None

    def _ensure_loaded(self):
        if self._smk_n is not None:
            return
        sd = SD(self._filepath, SDC.READ)
        self._smk_n = sd.select("NL_NISE_Extent").get().astype(np.int8)  # (lat, lon) = (721, 721)
        self._smk_s = sd.select("SL_NISE_Extent").get().astype(np.int8)
        sd.end()

    def _ease_grid(self, lon: float, lat: float) -> Tuple[int, int]:
        """Convert (lon, lat) to EASE-Grid (col, row) 0-based index."""
        phi = math.radians(lat)
        lam = math.radians(lon)
        pi4 = math.pi / 4.0

        if lat >= 0:
            rho = 2.0 * self.RG * math.sin(pi4 - phi / 2.0)
            r = self.R0 + rho * math.sin(lam)
            s = self.S0 + rho * math.cos(lam)
            grid = self._smk_n
        else:
            rho = 2.0 * self.RG * math.cos(pi4 - phi / 2.0)
            r = self.R0 + rho * math.sin(lam)
            s = self.S0 - rho * math.cos(lam)
            grid = self._smk_s

        row = max(0, min(self.ROWS - 1, int(round(s)) + 0))  # 0-based
        col = max(0, min(self.COLS - 1, int(round(r)) + 0))
        return grid, col, row  # col, row for Fortran-style (lon,lat) indexing

    def get(self, lon: float, lat: float) -> int:
        """Return NISE ice concentration (0-100) for a given (lon, lat)."""
        self._ensure_loaded()
        grid, col, row = self._ease_grid(lon, lat)
        return int(grid[col, row])  # Fortran: smk_n(col,row) = (lon,lat) order

    def get_array(self, lons: np.ndarray, lats: np.ndarray) -> np.ndarray:
        """Vectorized lookup — returns concentrations for arrays of (lon, lat)."""
        self._ensure_loaded()
        result = np.zeros_like(lons, dtype=np.int8)
        phi = np.radians(lats)
        lam = np.radians(lons)
        pi4 = math.pi / 4.0

        # North hemisphere
        nmask = lats >= 0
        if nmask.any():
            rho = 2.0 * self.RG * np.sin(pi4 - phi[nmask] / 2.0)
            r = np.clip((self.R0 + rho * np.sin(lam[nmask])).round().astype(int), 0, self.COLS - 1)
            s = np.clip((self.S0 + rho * np.cos(lam[nmask])).round().astype(int), 0, self.ROWS - 1)
            result[nmask] = self._smk_n[r, s]  # Fortran: smk_n(col,row) = (lon,lat)

        # South hemisphere
        smask = lats < 0
        if smask.any():
            rho = 2.0 * self.RG * np.cos(pi4 - phi[smask] / 2.0)
            r = np.clip((self.R0 + rho * np.sin(lam[smask])).round().astype(int), 0, self.COLS - 1)
            s = np.clip((self.S0 - rho * np.cos(lam[smask])).round().astype(int), 0, self.ROWS - 1)
            result[smask] = self._smk_s[r, s]  # Fortran: smk_s(col,row) = (lon,lat)

        return result


# =========================================================================
# 3. OISST (HDF5, 1440x720, 0.25-degree daily SST)
# =========================================================================

class OisstReader:
    """Reads OISST daily SST from HDF5 file."""

    def __init__(self, filepath: str):
        import h5py
        self._filepath = filepath
        self.nlat = 720
        self.nlon = 1440
        self.dlat = 0.25
        self.dlon = 0.25
        self.first_lat = -89.75
        self.first_lon = -179.75

    def _load(self):
        import h5py
        with h5py.File(self._filepath, 'r') as f:
            sst0 = f['sst'][:]  # (720, 1440) = (lat, lon), 0-360 lon convention
        # Rotate lon from [0,360) to [-180,180) — transpose to (lon, lat)
        half = self.nlon // 2
        sst0_t = sst0.T  # (1440, 720) = (lon, lat)
        sst1 = np.zeros_like(sst0_t)
        sst1[:half, :] = sst0_t[half:, :]
        sst1[half:, :] = sst0_t[:half, :]
        # Convert to Kelvin, replace fill values
        sst1 = np.where(sst1 > -20, sst1 + 273.15, -999.0)
        return sst1  # (lon, lat) = (1440, 720)

    def get(self, lon: float, lat: float) -> float:
        sst = self._load()
        my = int((lat - self.first_lat + 0.5 * self.dlat) / self.dlat)  # 0-based
        mx = int((lon - self.first_lon + 0.5 * self.dlon) / self.dlon)
        my = max(0, min(self.nlat - 1, my))
        mx = max(0, min(self.nlon - 1, mx))
        return float(sst[mx, my])

    def get_array(self, lons: np.ndarray, lats: np.ndarray) -> np.ndarray:
        sst = self._load()
        my = ((lats - self.first_lat + 0.5 * self.dlat) / self.dlat).astype(int)
        mx = ((lons - self.first_lon + 0.5 * self.dlon) / self.dlon).astype(int)
        my = np.clip(my, 0, self.nlat - 1)
        mx = np.clip(mx, 0, self.nlon - 1)
        return sst[mx, my]


# =========================================================================
# Convenience: load all ancillary data for a scene
# =========================================================================

def load_ancillary(scene_date: str = "20220803",
                   oisst_path: str = "/data/Data_minmin/oisst/sst.day.mean.20200401.hdf5",
                   project_root: str = None) -> dict:
    """Load all ancillary data needed for cloud mask processing.

    Args:
        scene_date: YYYYMMDD format date string
        oisst_path: Path to OISST HDF5 file
        project_root: Project root directory (default: auto-detect)

    Returns:
        Dict with 'eco', 'nise', 'oisst' readers
    """
    if project_root is None:
        project_root = _project_root()

    month = int(scene_date[4:6])
    eco = EcosystemReader(os.path.join(project_root, "coeff", "fylat_ecosystem.hdf"))
    nise = NiseReader(month)
    oisst = OisstReader(oisst_path)

    return {'eco': eco, 'nise': nise, 'oisst': oisst}
