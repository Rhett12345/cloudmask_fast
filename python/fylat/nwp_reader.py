"""NWP GRIB2 reader — replaces wgrib/ shell scripts with Python-native pipeline.

Reads GFS 0.25-degree GRIB2 files and converts to the flat binary format
that Fortran's read_NWP_arrays_0p25_41Layers() expects.

Usage:
    from fylat.nwp_reader import grib2_to_binary
    grib2_to_binary(grib_path, bin_path)
"""

import os
import subprocess
import sys
from typing import List

# --- Configuration -----------------------------------------------------------
# wgrib2 binary location (treated as tool dependency, like ifort/HDF5)
_WGRIB2 = "/opt/software/grib2/wgrib2/wgrib2"

# 41 pressure levels expected by Fortran (hPa)
_PRE_LEVELS = [
    0.01, 0.02, 0.04, 0.07, 0.1, 0.2, 0.4, 0.7,
    1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 15.0, 20.0, 30.0, 40.0,
    50.0, 70.0, 100.0, 150.0, 200.0, 250.0, 300.0, 350.0,
    400.0, 450.0, 500.0, 550.0, 600.0, 650.0, 700.0, 750.0,
    800.0, 850.0, 900.0, 925.0, 950.0, 975.0, 1000.0,
]

# CLWMR pressure levels (22 levels, 50 hPa and above)
_CLWMR_LEVELS = [
    50.0, 100.0, 150.0, 200.0, 250.0, 300.0, 350.0, 400.0,
    450.0, 500.0, 550.0, 600.0, 650.0, 700.0, 750.0,
    800.0, 850.0, 900.0, 925.0, 950.0, 975.0, 1000.0,
]


def _build_namelist() -> List[str]:
    """Build the 281-field NameList matching the Fortran binary layout.

    Field order (nwp_opt=10, GFS 0.25-degree 41-layer):
      1-13:   surface / special fields
      14-54:  TMP on 41 pressure levels
      55-95:  HGT on 41 pressure levels
      96-136: O3MR on 41 pressure levels
      137-177: RH on 41 pressure levels
      178-199: CLWMR on 22 pressure levels
      200-240: UGRD on 41 pressure levels
      241-281: VGRD on 41 pressure levels
    """
    names = [
        # 1-13: surface and special fields
        "PRES:surface",
        "PRMSL:mean sea level",
        "TMP:surface",
        "HGT:surface",
        "PRES:surface",               # albedo placeholder (dup of field 1)
        "TMP:0.995 sigma level",
        "RH:0.995 sigma level",
        "UGRD:0.995 sigma level",
        "VGRD:0.995 sigma level",
        "PWAT",
        "WEASD:surface",
        "TOZNE",
        "TMP:tropopause",
    ]

    # TMP on 41 pressure levels
    for p in _PRE_LEVELS:
        names.append(f"TMP:{p:g} mb")

    # HGT on 41 pressure levels
    for p in _PRE_LEVELS:
        names.append(f"HGT:{p:g} mb")

    # O3MR on 41 pressure levels
    for p in _PRE_LEVELS:
        names.append(f"O3MR:{p:g} mb")

    # RH on 41 pressure levels
    for p in _PRE_LEVELS:
        names.append(f"RH:{p:g} mb")

    # CLWMR on 22 pressure levels
    for p in _CLWMR_LEVELS:
        names.append(f"CLWMR:{p:g} mb")

    # UGRD on 41 pressure levels
    for p in _PRE_LEVELS:
        names.append(f"UGRD:{p:g} mb")

    # VGRD on 41 pressure levels
    for p in _PRE_LEVELS:
        names.append(f"VGRD:{p:g} mb")

    # Extra fields appended by shell script (10 m wind)
    names.append("UGRD:10 m above ground")
    names.append("VGRD:10 m above ground")

    return names


def grib2_to_binary(grib_path: str, bin_path: str) -> None:
    """Convert GFS GRIB2 file to flat float32 binary for Fortran consumption.

    Uses wgrib2 for GRIB2 decoding.  Replaces the shell scripts in wgrib/.

    Args:
        grib_path: Path to GFS GRIB2 file.
        bin_path: Output path for flat binary file.

    Raises:
        FileNotFoundError: If wgrib2 binary or GRIB2 file is missing.
        RuntimeError: If wgrib2 conversion fails.
    """
    if not os.path.exists(_WGRIB2):
        raise FileNotFoundError(f"wgrib2 not found at {_WGRIB2}")
    if not os.path.exists(grib_path):
        raise FileNotFoundError(f"GRIB2 file not found: {grib_path}")

    name_list = _build_namelist()

    # Remove any existing output file
    if os.path.exists(bin_path):
        os.remove(bin_path)

    # Build wgrib2 inventory once
    inv = subprocess.run(
        [_WGRIB2, "-set", "local_table", "1", "-s", grib_path],
        capture_output=True, text=True, timeout=120,
    )
    if inv.returncode != 0:
        raise RuntimeError(f"wgrib2 inventory failed: {inv.stderr}")
    inventory_lines = inv.stdout.split("\n")

    # Match each NameList entry against inventory and extract to binary
    matched_count = 0
    for name in name_list:
        # Find matching inventory lines
        matching = [line for line in inventory_lines if name in line]
        if not matching:
            continue

        # Build index file from matching line numbers
        # wgrib2 -i reads inventory records from stdin
        idx_input = "\n".join(matching)
        proc = subprocess.run(
            [_WGRIB2, "-i", grib_path, "-no_header", "-append", "-bin", bin_path],
            input=idx_input, capture_output=True, text=True, timeout=60,
        )
        if proc.returncode != 0:
            raise RuntimeError(
                f"wgrib2 extraction failed for '{name}': {proc.stderr}"
            )
        matched_count += 1

    # Verify output
    if not os.path.exists(bin_path):
        raise RuntimeError(f"Binary output not created: {bin_path}")

    expected_bytes = 1440 * 721 * 4  # one field
    actual_size = os.path.getsize(bin_path)
    n_fields = actual_size // expected_bytes
    print(f"  [NWP] {os.path.basename(grib_path)} -> "
          f"{os.path.basename(bin_path)}  "
          f"({n_fields} fields, {actual_size/1024/1024:.0f} MB)")


def generate_nwp_binary(nwp_grib1: str, nwp_grib2: str, nwp_path: str) -> str:
    """Generate the combined NWP binary from two GRIB2 files.

    Fortran expects a single .bin file containing both forecast times
    concatenated (with 'uv' suffix for the wind fields).

    Args:
        nwp_grib1: Path to first GRIB2 file (earlier forecast).
        nwp_grib2: Path to second GRIB2 file (later forecast).
        nwp_path: NWP data directory (where .bin file will be written).

    Returns:
        Path to the generated .bin file.
    """
    # Extract time info from GRIB filenames for bin naming
    # gfs.t06z.pgrb2.0p25.f018 -> extract '06z' and forecast hours
    import re

    def _extract_info(fp):
        basename = os.path.basename(fp)
        m = re.search(r't(\d{2})z.*\.f(\d{3})', basename)
        if m:
            return m.group(1) + 'z', m.group(2)
        return '00z', '000'

    t1_cycle, t1_lead = _extract_info(nwp_grib1)
    t2_cycle, t2_lead = _extract_info(nwp_grib2)

    # Name format: gfs0p25_41L_YYYYMMDD_HH_00_HH_00_uv
    # Match what Fortran's convert_grib_to_binary expects
    bin_name = f"gfs0p25_41L_{t1_cycle}_{t1_lead}_{t2_lead}_uv"
    bin_path = os.path.join(nwp_path, bin_name)

    if os.path.exists(bin_path):
        print(f"  [NWP] Binary exists, skipping: {bin_path}")
        return bin_path

    # Generate binary from first GRIB file
    grib2_to_binary(nwp_grib1, bin_path)

    # Append second GRIB file's data
    # wgrib2 -append adds to existing file
    name_list = _build_namelist()
    inv = subprocess.run(
        [_WGRIB2, "-set", "local_table", "1", "-s", nwp_grib2],
        capture_output=True, text=True, timeout=120,
    )
    if inv.returncode != 0:
        raise RuntimeError(f"wgrib2 inventory failed for {nwp_grib2}")
    inventory_lines = inv.stdout.split("\n")

    for name in name_list:
        matching = [line for line in inventory_lines if name in line]
        if not matching:
            continue
        idx_input = "\n".join(matching)
        subprocess.run(
            [_WGRIB2, "-i", nwp_grib2, "-no_header", "-append", "-bin", bin_path],
            input=idx_input, capture_output=True, text=True, timeout=60,
        )

    actual_size = os.path.getsize(bin_path)
    expected_bytes = 1440 * 721 * 4
    n_fields = actual_size // expected_bytes
    print(f"  [NWP] Combined binary: {os.path.basename(bin_path)} "
          f"({n_fields} fields, {actual_size/1024/1024:.0f} MB)")

    return bin_path
