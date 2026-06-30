"""NWP GRIB2 reader and Fortran-compatible binary writer.

The production Fortran core still reads a flat float32 binary file, but the
GRIB2 decoding now runs in Python through ECMWF ecCodes instead of the legacy
``wgrib/`` shell scripts.  ``cfgrib`` is exposed for higher-level Python users
that need xarray datasets for diagnostics, ML features, or future interpolation
work.
"""

from __future__ import annotations

from dataclasses import dataclass
import os
import re
import shutil
import subprocess
import tempfile
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

# --- Configuration -----------------------------------------------------------
_WGRIB2 = "/opt/software/grib2/wgrib2/wgrib2"
_FIELD_BYTES = 1440 * 721 * 4

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


@dataclass(frozen=True)
class FieldSpec:
    """One Fortran binary field in output order."""

    label: str
    short_names: Tuple[str, ...]
    type_of_level: Optional[str] = None
    level: Optional[float] = None
    name_contains: Tuple[str, ...] = ()


def _build_namelist() -> List[str]:
    """Build the legacy field label list used by the old wgrib scripts."""
    return [spec.label for spec in _build_field_specs()]


def _build_field_specs() -> List[FieldSpec]:
    """Build GRIB field specs in the exact order Fortran expects.

    The layout follows read_NWP_arrays_0p25_41Layers() for GFS 0.25-degree
    41-layer inputs.  The fifth field is the legacy albedo placeholder and is a
    duplicate of surface pressure, matching the old shell-script behavior.
    """
    def pressure_level(p: float) -> Tuple[str, float]:
        if p < 1.0:
            return "isobaricInPa", p * 100.0
        return "isobaricInhPa", p

    specs: List[FieldSpec] = [
        FieldSpec("PRES:surface", ("pres", "sp"), "surface"),
        FieldSpec("PRMSL:mean sea level", ("prmsl", "msl"), "meanSea"),
        FieldSpec("TMP:surface", ("t", "tmp"), "surface", name_contains=("Temperature",)),
        FieldSpec("HGT:surface", ("gh", "hgt", "orog"), "surface", name_contains=("Geopotential", "height")),
        FieldSpec("PRES:surface", ("pres", "sp"), "surface"),
        # ecCodes reports the legacy "0.995 sigma level" records as sigma
        # level 1 for these GFS files.
        FieldSpec("TMP:0.995 sigma level", ("t", "tmp"), "sigma", 1),
        FieldSpec("RH:0.995 sigma level", ("r", "rh"), "sigma", 1),
        FieldSpec("UGRD:0.995 sigma level", ("u", "ugrd"), "sigma", 1),
        FieldSpec("VGRD:0.995 sigma level", ("v", "vgrd"), "sigma", 1),
        FieldSpec("PWAT", ("pwat", "tcw"), name_contains=("Precipitable water",)),
        FieldSpec("WEASD:surface", ("sdwe", "weasd"), "surface"),
        FieldSpec("TOZNE", ("tozne", "tco3"), name_contains=("Total ozone",)),
        FieldSpec("TMP:tropopause", ("t", "tmp"), "tropopause", name_contains=("Temperature",)),
    ]

    for p in _PRE_LEVELS:
        level_type, level = pressure_level(p)
        specs.append(FieldSpec(f"TMP:{p:g} mb", ("t", "tmp"), level_type, level))
    for p in _PRE_LEVELS:
        level_type, level = pressure_level(p)
        specs.append(FieldSpec(f"HGT:{p:g} mb", ("gh", "hgt"), level_type, level))
    for p in _PRE_LEVELS:
        level_type, level = pressure_level(p)
        specs.append(FieldSpec(f"O3MR:{p:g} mb", ("o3mr",), level_type, level))
    for p in _PRE_LEVELS:
        level_type, level = pressure_level(p)
        specs.append(FieldSpec(f"RH:{p:g} mb", ("r", "rh"), level_type, level))
    for p in _CLWMR_LEVELS:
        specs.append(FieldSpec(f"CLWMR:{p:g} mb", ("clwmr", "clwat"), "isobaricInhPa", p))
    for p in _PRE_LEVELS:
        level_type, level = pressure_level(p)
        specs.append(FieldSpec(f"UGRD:{p:g} mb", ("u", "ugrd"), level_type, level))
    for p in _PRE_LEVELS:
        level_type, level = pressure_level(p)
        specs.append(FieldSpec(f"VGRD:{p:g} mb", ("v", "vgrd"), level_type, level))

    specs.append(FieldSpec("UGRD:10 m above ground", ("10u", "u", "ugrd"), "heightAboveGround", 10))
    specs.append(FieldSpec("VGRD:10 m above ground", ("10v", "v", "vgrd"), "heightAboveGround", 10))
    return specs


def read_gfs_grib2(path: str, filter_by_keys: Optional[Dict[str, object]] = None):
    """Open a GFS GRIB2 file as cfgrib/xarray datasets.

    ``cfgrib.open_datasets`` is used instead of a single ``open_dataset`` call
    because GFS files contain many level types and variable groups.
    """
    import cfgrib

    backend_kwargs = {"indexpath": ""}
    if filter_by_keys:
        backend_kwargs["filter_by_keys"] = filter_by_keys
        return cfgrib.open_dataset(path, backend_kwargs=backend_kwargs)
    return cfgrib.open_datasets(path, backend_kwargs=backend_kwargs)


def _level_matches(actual: object, expected: Optional[float]) -> bool:
    if expected is None:
        return True
    try:
        return abs(float(actual) - float(expected)) < 1.0e-6
    except (TypeError, ValueError):
        return False


def _text_contains_any(text: str, needles: Iterable[str]) -> bool:
    folded = text.lower()
    return any(needle.lower() in folded for needle in needles)


def _matches(meta: Dict[str, object], spec: FieldSpec) -> bool:
    short_name = str(meta.get("shortName", "")).lower()
    param = str(meta.get("paramId", ""))
    name = str(meta.get("name", ""))
    type_of_level = str(meta.get("typeOfLevel", ""))

    if short_name not in spec.short_names and param not in spec.short_names:
        if not spec.name_contains or not _text_contains_any(name, spec.name_contains):
            return False
    if spec.type_of_level and type_of_level != spec.type_of_level:
        return False
    return _level_matches(meta.get("level"), spec.level)


def _get_key(gid, key: str):
    from eccodes import CodesInternalError, codes_get

    try:
        return codes_get(gid, key)
    except CodesInternalError:
        return None


def _message_meta(gid) -> Dict[str, object]:
    return {
        "shortName": _get_key(gid, "shortName"),
        "name": _get_key(gid, "name"),
        "paramId": _get_key(gid, "paramId"),
        "typeOfLevel": _get_key(gid, "typeOfLevel"),
        "level": _get_key(gid, "level"),
    }


def _write_values(path: str, values: np.ndarray) -> None:
    values = np.asarray(values, dtype="<f4")
    with open(path, "wb") as f:
        values.tofile(f)


def _grib2_to_binary_eccodes(grib_path: str, bin_path: str) -> Tuple[int, List[str]]:
    """Convert one GRIB2 file to a Fortran-compatible binary with ecCodes."""
    from eccodes import codes_get_values, codes_grib_new_from_file, codes_release

    specs = _build_field_specs()
    matched_paths: Dict[int, str] = {}
    missing: List[str] = []
    out_dir = os.path.dirname(os.path.abspath(bin_path)) or "."
    os.makedirs(out_dir, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="fylat_nwp_", dir=out_dir) as tmp:
        with open(grib_path, "rb") as f:
            while True:
                gid = codes_grib_new_from_file(f)
                if gid is None:
                    break
                try:
                    meta = _message_meta(gid)
                    matching_indexes = [
                        idx for idx, spec in enumerate(specs)
                        if idx not in matched_paths and _matches(meta, spec)
                    ]
                    if matching_indexes:
                        values = np.asarray(codes_get_values(gid), dtype="<f4")
                        for idx in matching_indexes:
                            field_path = os.path.join(tmp, f"field_{idx:03d}.bin")
                            _write_values(field_path, values)
                            matched_paths[idx] = field_path
                finally:
                    codes_release(gid)

        if os.path.exists(bin_path):
            os.remove(bin_path)
        with open(bin_path, "wb") as out:
            for idx, spec in enumerate(specs):
                field_path = matched_paths.get(idx)
                if field_path is None:
                    missing.append(spec.label)
                    continue
                with open(field_path, "rb") as src:
                    shutil.copyfileobj(src, out)

    return len(matched_paths), missing


def _grib2_to_binary_wgrib2(grib_path: str, bin_path: str) -> None:
    """Compatibility backend for explicit emergency use."""
    if not os.path.exists(_WGRIB2):
        raise FileNotFoundError(f"wgrib2 not found at {_WGRIB2}")

    name_list = _build_namelist()
    if os.path.exists(bin_path):
        os.remove(bin_path)

    inv = subprocess.run(
        [_WGRIB2, "-set", "local_table", "1", "-s", grib_path],
        capture_output=True, text=True, timeout=120,
    )
    if inv.returncode != 0:
        raise RuntimeError(f"wgrib2 inventory failed: {inv.stderr}")
    inventory_lines = inv.stdout.split("\n")

    for name in name_list:
        matching = [line for line in inventory_lines if name in line]
        if not matching:
            continue
        proc = subprocess.run(
            [_WGRIB2, "-i", grib_path, "-no_header", "-append", "-bin", bin_path],
            input="\n".join(matching), capture_output=True, text=True, timeout=60,
        )
        if proc.returncode != 0:
            raise RuntimeError(f"wgrib2 extraction failed for '{name}': {proc.stderr}")


def grib2_to_binary(grib_path: str, bin_path: str, backend: Optional[str] = None) -> None:
    """Convert a GFS GRIB2 file to flat float32 binary for Fortran.

    Backend selection:
      - ``eccodes`` (default): Python ecCodes decoder, no legacy shell scripts.
      - ``wgrib2``: explicit compatibility backend via ``FYLAT_NWP_BACKEND=wgrib2``.
    """
    if not os.path.exists(grib_path):
        raise FileNotFoundError(f"GRIB2 file not found: {grib_path}")

    backend = (backend or os.environ.get("FYLAT_NWP_BACKEND", "eccodes")).lower()
    if backend == "wgrib2":
        _grib2_to_binary_wgrib2(grib_path, bin_path)
        matched = os.path.getsize(bin_path) // _FIELD_BYTES
        missing: Sequence[str] = ()
    elif backend == "eccodes":
        matched, missing = _grib2_to_binary_eccodes(grib_path, bin_path)
    else:
        raise ValueError(f"Unsupported NWP backend: {backend}")

    if not os.path.exists(bin_path):
        raise RuntimeError(f"Binary output not created: {bin_path}")

    actual_size = os.path.getsize(bin_path)
    n_fields = actual_size // _FIELD_BYTES if _FIELD_BYTES else 0
    print(f"  [NWP:{backend}] {os.path.basename(grib_path)} -> "
          f"{os.path.basename(bin_path)}  "
          f"({n_fields} fields, {actual_size/1024/1024:.0f} MB)")
    if missing:
        preview = ", ".join(missing[:8])
        suffix = "..." if len(missing) > 8 else ""
        print(f"  [NWP:{backend}] Warning: {len(missing)} field(s) missing: "
              f"{preview}{suffix}")
    if matched == 0:
        raise RuntimeError(f"No expected NWP fields matched in {grib_path}")


def fortran_nwp_binary_name(date: str, valid_hour: int) -> str:
    """Return the legacy single-time binary name expected by Fortran."""
    return f"gfs0p25_41L_{date}_{valid_hour:02d}_00"


def generate_fortran_nwp_binary(grib_path: str, nwp_path: str,
                                date: str, valid_hour: int) -> str:
    """Generate one legacy Fortran-compatible NWP binary for a valid hour."""
    bin_path = os.path.join(nwp_path, fortran_nwp_binary_name(date, valid_hour))
    if os.path.exists(bin_path):
        print(f"  [NWP] Binary exists, skipping: {bin_path}")
        return bin_path
    grib2_to_binary(grib_path, bin_path)
    return bin_path
