"""YAML configuration loader and namelist (.nml) generator for FYLAT.

Replaces the manual Fortran namelist editing workflow with structured YAML configs.
The Fortran executable still reads .nml files — we generate them from YAML.
"""

import os
import shutil
from copy import deepcopy
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Optional

import yaml


# ---------------------------------------------------------------------------
# Default paths (machine-specific — override via YAML or env vars)
# ---------------------------------------------------------------------------
DEFAULT_CODE_ROOT = os.path.dirname(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
)


def _deep_merge(base: dict, override: dict) -> dict:
    """Recursively merge override into base. Lists are replaced, not merged."""
    result = deepcopy(base)
    for key, value in override.items():
        if key in result and isinstance(result[key], dict) and isinstance(value, dict):
            result[key] = _deep_merge(result[key], value)
        else:
            result[key] = deepcopy(value)
    return result


def load_config(
    scene_path: Optional[str] = None,
    base_path: Optional[str] = None,
) -> Dict[str, Any]:
    """Load the merged YAML configuration.

    Resolution order (later overrides earlier):
      1. config/default.yaml  (shipped defaults)
      2. config/scenes/<scene>.yaml  (scene-specific overrides)

    If scene_path is given directly, only that file is loaded (no merge).
    """
    project_root = Path(DEFAULT_CODE_ROOT)

    if base_path:
        base_file = Path(base_path)
    else:
        base_file = project_root / "config" / "default.yaml"

    if not base_file.exists():
        raise FileNotFoundError(f"Base config not found: {base_file}")

    with open(base_file) as f:
        cfg = yaml.safe_load(f)

    if scene_path:
        scene_file = Path(scene_path)
        if scene_file.exists():
            with open(scene_file) as f:
                scene_cfg = yaml.safe_load(f)
            cfg = _deep_merge(cfg, scene_cfg)

    return cfg


# ---------------------------------------------------------------------------
# Path derivation
# ---------------------------------------------------------------------------
def _derive_paths(cfg: Dict[str, Any]) -> Dict[str, str]:
    """Auto-derive I/O file paths from scene date/time/calibration."""
    scene = cfg.get("scene", {})
    date = scene.get("date", "")
    time = scene.get("time", "")
    cal = scene.get("calibration", "business")
    cal_suffix = "BUSINESS" if cal == "business" else "RECALI"

    paths = cfg.get("paths", {})
    l1b_data = paths.get("l1b_data", "")
    nwp_data = paths.get("nwp_data", "")
    oisst_data = paths.get("oisst_data", "")
    output = paths.get("output", "")

    input_cfg = cfg.get("input", {})
    output_cfg = cfg.get("output", {})

    # L1B file
    l1b_file = input_cfg.get("l1b_file", "")
    if not l1b_file:
        l1b_file = os.path.join(
            l1b_data, date,
            f"FY3D_MERSI_GBAL_L1_{date}_{time}_1000M_MS.HDF"
        )

    # GEO file
    geo_file = input_cfg.get("geo_file", "")
    if not geo_file:
        geo_file = os.path.join(
            l1b_data, date,
            f"FY3D_MERSI_GBAL_L1_{date}_{time}_GEO1K_MS.HDF"
        )

    # NWP grib files
    nwp_grib1 = input_cfg.get("nwp_grib1", "")
    nwp_grib2 = input_cfg.get("nwp_grib2", "")

    # OISST file
    oisst_file = input_cfg.get("oisst_file", "")

    # Output files
    out_dir = os.path.join(output, date)
    base_out = f"FY3D_MERSI_ORBT_L2_{{}}_MLT_NUL_{date}_{time}_1000M_MS_{cal_suffix}.HDF"
    base_out_5km = f"FY3D_MERSI_ORBT_L2_{{}}_MLT_NUL_{date}_{time}_5000M_MS_{cal_suffix}.HDF"
    base_out_intermed = f"FY3D_MERSI_ORBT_L2_XXX_MLT_NUL_{date}_{time}_INTERMED_{cal_suffix}.HDF"

    def _out(key, tmpl, use_5km=False):
        val = output_cfg.get(key, "")
        if not val:
            t = base_out_5km if use_5km else base_out
            val = os.path.join(out_dir, t.format(key.upper().replace("_FILE", "")))
        return val

    clm_file = _out("clm_file", base_out)
    cla_file = _out("cla_file", base_out, use_5km=True)
    clp_file = _out("clp_file", base_out)
    ctp_file = _out("ctp_file", base_out)
    cot_file = _out("cot_file", base_out)
    con_file = _out("con_file", base_out)
    sst_file = _out("sst_file", base_out)
    intermediate_file = output_cfg.get("intermediate_file", "")
    if not intermediate_file:
        intermediate_file = os.path.join(out_dir, base_out_intermed)

    return {
        "code_root_path": paths.get("code_root", DEFAULT_CODE_ROOT),
        "L1b_data_path": os.path.join(l1b_data, date, ""),
        "nwp_data_path": os.path.join(nwp_data, date, ""),
        "oisst_data_path": oisst_data,
        "fy3_mersi_GEO_data": geo_file,
        "fy3_mersi_L1b_data": l1b_file,
        "fy3_mersi_CLM_data": clm_file,
        "fy3_mersi_CLA_data": cla_file,
        "fy3_mersi_CLP_data": clp_file,
        "fy3_mersi_CTP_data": ctp_file,
        "fy3_mersi_COT_data": cot_file,
        "fy3_mersi_CON_data": con_file,
        "fy3_mersi_SST_data": sst_file,
        "fy3_intermediate": intermediate_file,
        "nwp_grib_data1": nwp_grib1,
        "nwp_grib_data2": nwp_grib2,
        "oisst_data": oisst_file,
    }


# ---------------------------------------------------------------------------
# Namelist generation
# ---------------------------------------------------------------------------
def generate_legacy_namelist(cfg: Dict[str, Any]) -> str:
    """Generate a legacy Fortran namelist (.nml) string from YAML config."""
    sensor = cfg.get("sensor", {})
    paths = cfg.get("paths", {})
    algo = cfg.get("algorithms", {})

    pp = _derive_paths(cfg)

    def _bool_to_int(val):
        return 1 if val else 0

    lines = ["&config"]
    lines.append(f"  fylat_sensor_id     = {sensor.get('id', 21)},")
    lines.append(f'  code_root_path      = "{pp["code_root_path"]}",')
    lines.append(f'  L1b_data_path       = "{pp["L1b_data_path"]}",')
    lines.append(f'  nwp_data_path       = "{pp["nwp_data_path"]}",')
    lines.append(f'  oisst_data_path     = "{pp["oisst_data_path"]}",')
    lines.append(f'  fy3_mersi_GEO_data  = "{pp["fy3_mersi_GEO_data"]}",')
    lines.append(f'  fy3_mersi_L1b_data  = "{pp["fy3_mersi_L1b_data"]}",')
    lines.append(f'  fy3_mersi_CLM_data  = "{pp["fy3_mersi_CLM_data"]}",')
    lines.append(f'  fy3_mersi_CLA_data  = "{pp["fy3_mersi_CLA_data"]}",')
    lines.append(f'  fy3_mersi_CLP_data  = "{pp["fy3_mersi_CLP_data"]}",')
    lines.append(f'  fy3_mersi_CTP_data  = "{pp["fy3_mersi_CTP_data"]}",')
    lines.append(f'  fy3_mersi_COT_data  = "{pp["fy3_mersi_COT_data"]}",')
    lines.append(f'  fy3_mersi_CON_data  = "{pp["fy3_mersi_CON_data"]}",')
    lines.append(f'  fy3_mersi_SST_data  = "{pp["fy3_mersi_SST_data"]}",')
    lines.append(f'  fy3_intermediate    = "{pp["fy3_intermediate"]}",')
    lines.append(f"  fylat_nwp_opt       = {sensor.get('nwp_opt', 10)},")
    lines.append(f"  fylat_rtm_opt       = {sensor.get('rtm_opt', 1)},")
    lines.append(f'  nwp_grib_data1      = "{pp["nwp_grib_data1"]}",')
    lines.append(f'  nwp_grib_data2      = "{pp["nwp_grib_data2"]}",')
    lines.append(f'  oisst_data          = "{pp["oisst_data"]}",')
    lines.append(f"  cloudmask_id        = {_bool_to_int(algo.get('cloudmask', True))},")
    lines.append(f"  cloudamount_id      = {_bool_to_int(algo.get('cloudamount', False))},")
    lines.append(f"  cloudphase_id       = {_bool_to_int(algo.get('cloudphase', False))},")
    lines.append(f"  cloudtopz_id        = {_bool_to_int(algo.get('cloudtopz', False))},")
    lines.append(f"  cloudtau_day_id     = {_bool_to_int(algo.get('cloudtau_day', False))},")
    lines.append(f"  cloudtau_night_id   = {_bool_to_int(algo.get('cloudtau_night', False))},")
    lines.append(f"  cloudtypeII_id      = {_bool_to_int(algo.get('cloudtypeII', False))},")
    lines.append(f"  surface_sst_id      = {_bool_to_int(algo.get('surface_sst', False))},")
    lines.append(f"  write_inter_id      = {_bool_to_int(algo.get('write_intermediate', False))}/")
    return "\n".join(lines)


def generate_namelist(cfg: Dict[str, Any]) -> str:
    """Compatibility alias for the legacy Fortran backend only."""
    return generate_legacy_namelist(cfg)


def write_namelist(cfg: Dict[str, Any], output_path: str) -> str:
    """Write the generated .nml file and return its path."""
    nml_content = generate_legacy_namelist(cfg)
    with open(output_path, "w") as f:
        f.write(nml_content + "\n")
    return output_path


def ensure_legacy_output_dirs(cfg: Dict[str, Any]) -> None:
    """Create output directories required by the legacy Fortran writer."""
    pp = _derive_paths(cfg)
    output_keys = [
        "fy3_mersi_CLM_data",
        "fy3_mersi_CLA_data",
        "fy3_mersi_CLP_data",
        "fy3_mersi_CTP_data",
        "fy3_mersi_COT_data",
        "fy3_mersi_CON_data",
        "fy3_mersi_SST_data",
        "fy3_intermediate",
    ]
    for key in output_keys:
        path = pp.get(key, "")
        if path:
            Path(path).expanduser().parent.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Calibration mode
# ---------------------------------------------------------------------------
def setup_calibration(cfg: Dict[str, Any], code_root: str) -> None:
    """Configure calibration mode for the Fortran runtime.

    Uses Python calibration module — no more manual .xcfg file management.
    The calibration field in YAML can be:
      - "business": use HDF5 built-in coefficients
      - "recali": auto-derive YYYYMM from scene date
      - "YYYYMM": use month-specific recalibration (e.g., "202208")
      - "YYYYMMDD": truncated to YYYYMM for lookup
    """
    from fylat.calibration import setup_cal_mode, list_calibrations

    cal = cfg.get("scene", {}).get("calibration", "business")

    # Normalize aliases to YYYYMM keys
    if cal == "recali":
        date = cfg.get("scene", {}).get("date", "")
        cal = date[:6] if len(date) >= 6 else cal
    elif len(cal) == 8 and cal.isdigit():
        cal = cal[:6]

    available = list_calibrations()
    if cal not in available:
        print(f"  Warning: unknown calibration '{cal}', falling back to 'business'")
        print(f"  Available: {available[:10]}...")
        cal = "business"
    setup_cal_mode(cal, code_root)


# ---------------------------------------------------------------------------
# Convenience: run a single scene
# ---------------------------------------------------------------------------
def run_scene(scene_yaml_path: str, base_yaml_path: Optional[str] = None,
              calibration: Optional[str] = None,
              output_path: Optional[str] = None) -> int:
    """Full workflow for a single scene: load config → write .nml → run Fortran.

    Returns the exit code of the Fortran executable.
    """
    cfg = load_config(scene_path=scene_yaml_path, base_path=base_yaml_path)
    if calibration:
        cfg.setdefault("scene", {})["calibration"] = calibration
    if output_path:
        cfg.setdefault("paths", {})["output"] = output_path
    code_root = cfg.get("paths", {}).get("code_root", DEFAULT_CODE_ROOT)

    # Setup calibration
    setup_calibration(cfg, code_root)

    # Determine scene info for temp file naming
    scene = cfg.get("scene", {})
    date = scene.get("date", "unknown")
    time = scene.get("time", "0000")
    nml_path = os.path.join(code_root, f"temp_fy3d_config_{date}_{time}.nml")

    # Write namelist
    ensure_legacy_output_dirs(cfg)
    write_namelist(cfg, nml_path)
    print(f"  Config written to: {nml_path}")

    # Run Fortran
    import subprocess

    exe_path = os.path.join(code_root, "fylat_FY3_MERSI_II_PGS")
    if not os.path.exists(exe_path):
        raise FileNotFoundError(f"Executable not found: {exe_path}. Run build.sh first.")

    print(f"  Running: {exe_path} {nml_path}")
    result = subprocess.run(
        [exe_path, nml_path],
        cwd=code_root,
    )
    return result.returncode
