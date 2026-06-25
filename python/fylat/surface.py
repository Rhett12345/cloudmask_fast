"""Surface emissivity and albedo loaders — replaces Fortran get_ancil_data_module.f90
emissivity and albedo reading functions. Reads HDF4 files directly.
"""

import os
from pathlib import Path
from typing import Tuple

import numpy as np
from pyhdf.SD import SD, SDC


def _project_root() -> str:
    return str(Path(__file__).parent.parent.parent)


# =========================================================================
# Surface Emissivity (HDF4, 7200x3600, 0.05-degree, monthly)
# =========================================================================

# Mapping: month → day-of-year for emissivity files
EMISS_MONTH_DAYS = [1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335]

# Which SDS index maps to which IR channel
# emiss7→3.8um, emiss10→7.3um, emiss11→8.6um, emiss14→11um, emiss15→12um
EMISS_CHANNELS = {
    '38': 7, '40': 7,   # both 3.8 and 4.0 use emiss7
    '73': 10,
    '86': 11,
    '11': 14,
    '12': 15,
}

# Fallback water surface emissivity (from Fortran)
EMISS_WATER = {
    '38': 0.978, '40': 0.978,
    '73': 0.993, '86': 0.988,
    '11': 0.993, '12': 0.988,
}


class EmissivityReader:
    """Reads global surface emissivity from HDF4 monthly files."""

    def __init__(self, month: int, filepath: str = None):
        if filepath is None:
            jday = EMISS_MONTH_DAYS[month - 1]
            fname = f"global_emiss_intABI_2005{jday:03d}.hdf"
            filepath = os.path.join(_project_root(), "coeff", "sfc_emiss", fname)
        self._filepath = filepath
        self._data = {}
        self.nlon = 7200
        self.nlat = 3600
        self.dlon = 0.05
        self.dlat = 0.05
        self.first_lon = -179.975
        self.first_lat = 89.975

    def _load_channel(self, sds_index: int) -> np.ndarray:
        if sds_index in self._data:
            return self._data[sds_index]
        sd = SD(self._filepath, SDC.READ)
        sds_name = f"emiss{sds_index}"
        sds = sd.select(sds_name)
        raw = sds.get().astype(np.float32)  # int16 → float32
        # Apply scale/offset
        scale = sds.attributes().get('scale_factor', 1.0)
        offset = sds.attributes().get('add_offset', 0.0)
        emiss = raw * scale + offset
        sd.end()
        self._data[sds_index] = emiss
        return emiss

    def get(self, lon: float, lat: float, channel: str) -> float:
        """Get emissivity for a given (lon, lat) and channel name ('38','73','86','11','12')."""
        sds_idx = EMISS_CHANNELS.get(channel)
        if sds_idx is None:
            return EMISS_WATER.get(channel, 0.98)
        emiss = self._load_channel(sds_idx)
        nx = int((lon - self.first_lon) / self.dlon)
        ny = int((self.first_lat - lat) / self.dlat)
        nx = max(0, min(self.nlon - 1, nx))
        ny = max(0, min(self.nlat - 1, ny))
        val = float(emiss[ny, nx])
        if val <= 0.0 or val > 1.0:
            return EMISS_WATER.get(channel, 0.98)
        return val


# =========================================================================
# Surface Albedo (HDF4, 5400x2700, ~0.067-degree, every 16 days)
# =========================================================================

ALB_DAYS = [1, 17, 33, 49, 65, 81, 97, 113, 129, 145, 161, 177,
            193, 209, 225, 241, 257, 273, 289, 305, 321, 337, 353]

ALB_CHANNELS = {
    '066': (2, 'Albedo_Map_0.659', '0.659'),   # 0.659 um
    '086': (3, 'Albedo_Map_0.858', '0.858'),   # 0.858 um
    '124': (4, 'Albedo_Map_1.24', '1.24'),     # 1.24 um
    '164': (5, 'Albedo_Map_1.64', '1.64'),     # 1.64 um
    '213': (6, 'Albedo_Map_2.13', '2.13'),     # 2.13 um
}


class AlbedoReader:
    """Reads global white-sky albedo from HDF4 16-day files."""

    def __init__(self, day_of_year: int, filepath: str = None):
        # Find closest 16-day period
        closest = min(ALB_DAYS, key=lambda d: abs(d - day_of_year))
        self._jday = closest
        self._base_dir = os.path.join(_project_root(), "coeff", "sfc_albedo")
        self._data = {}
        self.nlon = 5400
        self.nlat = 2700
        self.dlon = 0.06666
        self.dlat = 0.06666
        self.first_lon = -179.96
        self.first_lat = 89.96

    def _load(self, wavelength: str) -> np.ndarray:
        if wavelength in self._data:
            return self._data[wavelength]
        idx, sds_name, wl_str = ALB_CHANNELS[wavelength]
        fname = f"AlbMap.WS.c004.v2.0.00-04.{self._jday:03d}.{wl_str}_x4.hdf"
        fpath = os.path.join(self._base_dir, fname)
        sd = SD(fpath, SDC.READ)
        sds = sd.select(sds_name)
        raw = sds.get().astype(np.float32)
        albedo = raw / 10.0  # hardcoded scale factor from Fortran
        sd.end()
        self._data[wavelength] = albedo
        return albedo

    def get(self, lon: float, lat: float, wavelength: str) -> float:
        """Get albedo for given (lon, lat) and wavelength ('066','086','124','164','213')."""
        if wavelength not in ALB_CHANNELS:
            return -999.0
        alb = self._load(wavelength)
        nx = int((lon - self.first_lon) / self.dlon)
        ny = int((self.first_lat - lat) / self.dlat)
        nx = max(0, min(self.nlon - 1, nx))
        ny = max(0, min(self.nlat - 1, ny))
        val = float(alb[ny, nx])
        if val > 100.0 or val <= 0.0:
            return -999.0
        return val
