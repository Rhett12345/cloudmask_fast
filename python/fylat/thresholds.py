"""FY-3D MERSI-II cloud mask threshold loader.

Reads the Fortran-format threshold file and provides Python dicts.
"""

import re
from pathlib import Path
from typing import Any, Dict, List, Optional


def _parse_float_list(line: str) -> List[float]:
    """Parse comma-separated float values from a threshold file line."""
    values = []
    line = re.sub(r'!.*$', '', line)
    for token in line.replace(',', ' ').split():
        try:
            values.append(float(token))
        except ValueError:
            pass
    return values


def load_thresholds(
    filepath: Optional[str] = None,
) -> Dict[str, Any]:
    """Load cloud mask thresholds from the Fortran threshold file.

    Returns a flat dict mapping parameter names to lists of float values.
    """
    if filepath is None:
        project_root = Path(__file__).parent.parent.parent
        filepath = project_root / "coeff" / "fylat_thresholds.mersi.aqua.v8"

    with open(filepath) as f:
        lines = f.readlines()

    thresholds: Dict[str, List[float]] = {}

    for line in lines:
        stripped = line.strip()
        if not stripped:
            continue
        # Skip pure comment lines (group headers, RCS, version)
        if stripped.startswith('!') or stripped.startswith('rcs_id') or stripped.startswith('thresholds_file_ver'):
            continue

        # Parse "name : val1, val2, ..."
        if ':' in stripped:
            name_part, _, val_part = stripped.partition(':')
            name = name_part.strip()
            values = _parse_float_list(val_part)
            if values and name:
                thresholds[name] = values

    return thresholds


def get_scene_thresholds(scene: str,
                         filepath: Optional[str] = None) -> Dict[str, List[float]]:
    """Get thresholds needed for a specific scene type.

    Args:
        scene: Scene type name (e.g., 'ocean_day', 'polar_nite_snow').
        filepath: Optional path to threshold file.

    Returns:
        Dict of parameter_name → [values] for that scene.
    """
    all_thr = load_thresholds(filepath)

    # Define which thresholds each scene needs
    scene_keys = {
        'ocean_day': [
            'dobt11', 'do11_12hi', 'do11_4lo', 'doref2', 'doref3',
            'dovrathi', 'dovratlo', 'dotci',
            'pfmft_11maxthre', 'pfmft_btd_min', 'pfmft_ocean', 'nfmft_ocean', 'nfmft_maxthre',
            'snglntv', 'snglntvcl', 'snglntvch', 'sg_tbdfl', 'sg_tbdfh', 'snglrat',
            'snglnt0', 'snglnt10', 'snglnt20', 'snglnt_bounds',
            'dovar11',
        ],
        'ocean_nite': [
            'nobt11', 'no11_12hi', 'no11_4lo', 'no86_73', 'no_11var',
            'pfmft_11maxthre', 'pfmft_btd_min', 'pfmft_ocean', 'nfmft_ocean', 'nfmft_maxthre',
            'dovar11',
        ],
        'land_day': [
            'dl11_12hi', 'dl11_4lo', 'dlref1', 'dlref3', 'dlvrat', 'dltci',
            'pfmft_11maxthre', 'pfmft_btd_min', 'pfmft_land', 'pfmft_cold',
            'nfmft_land', 'nfmft_maxthre',
        ],
        'land_nite': [
            'nl4_12hi', 'nl7_11s', 'nl11_12hi',
            'nl_11_4l', 'nl_11_4h', 'nl_11_4m', 'bt_diff_bounds',
            'pfmft_11maxthre', 'pfmft_btd_min', 'pfmft_land', 'pfmft_cold',
            'nfmft_land', 'nfmft_maxthre',
        ],
        'day_snow': [
            'ds11_12hi', 'ds11_12adj', 'ds4_11', 'ds4_11hel', 'dsref3', 'dstci',
            'pfmft_11maxthre', 'pfmft_btd_min', 'pfmft_snow', 'nfmft_snow', 'nfmft_maxthre',
        ],
        'nite_snow': [
            'ns11_12hi', 'ns11_12adj', 'ns11_4lo', 'ns4_12hi',
            'pfmft_11maxthre', 'pfmft_btd_min', 'pfmft_snow', 'nfmft_snow', 'nfmft_maxthre',
        ],
        'polar_day_land': [
            'pdl11_12hi', 'pdl11_4lo', 'pdlref1', 'pdlvrat', 'pdlref3', 'pdltci',
            'pfmft_11maxthre', 'pfmft_btd_min', 'pfmft_land', 'pfmft_cold', 'nfmft_land', 'nfmft_maxthre',
        ],
        'polar_day_ocean': [
            'pdobt11', 'pdo11_12hi', 'pdo11_4lo', 'pdoref2', 'pdoref3',
            'pdovrathi', 'pdovratlo', 'pdotci',
            'pfmft_11maxthre', 'pfmft_btd_min', 'pfmft_ocean', 'nfmft_ocean', 'nfmft_maxthre',
            'snglntv', 'snglntvcl', 'snglntvch', 'sg_tbdfl', 'sg_tbdfh', 'snglrat',
            'snglnt0', 'snglnt10', 'snglnt20', 'snglnt_bounds',
        ],
        'polar_day_snow': [
            'dps11_12hi', 'dps11_12adj', 'dps4_11l', 'dps4_11h', 'dps4_11m1',
            'bt_11_bnds3', 'dpsref3', 'dpstci',
            'pfmft_11maxthre', 'pfmft_btd_min', 'pfmft_snow', 'nfmft_snow', 'nfmft_maxthre',
        ],
        'polar_nite_land': [
            'pnl11_12hi', 'pn_11_4l', 'pn_11_4h', 'pn_11_4m1',
            'pn_7_11l', 'pn_7_11h', 'pn_7_11m1', 'pn_7_11m2', 'pn_7_11m3',
            'pn_4_12l', 'pn_4_12h', 'pn_4_12m1',
            'bt_11_bounds', 'bt_11_bnds2',
            'pfmft_11maxthre', 'pfmft_btd_min', 'pfmft_land', 'pfmft_cold', 'nfmft_land', 'nfmft_maxthre',
        ],
        'polar_nite_ocean': [
            'pnobt11', 'pno11_12hi', 'pno11_4lo', 'pno86_73', 'pno_11var',
            'pfmft_11maxthre', 'pfmft_btd_min', 'pfmft_ocean', 'nfmft_ocean', 'nfmft_maxthre',
            'dovar11',
        ],
        'polar_nite_snow': [
            'pns11_12hi', 'pn11_12adj',
            'pn_11_4l', 'pn_11_4h', 'pn_11_4m1',
            'pn_7_11l', 'pn_7_11h', 'pn_7_11m1', 'pn_7_11m2', 'pn_7_11m3',
            'pn_7_11lw', 'pn_7_11hw', 'pn_7_11m1w', 'pn_7_11m2w', 'pn_7_11m3w',
            'pn_4_12l', 'pn_4_12h', 'pn_4_12m1',
            'bt_11_bounds', 'bt_11_bnds2',
            'pfmft_11maxthre', 'pfmft_btd_min', 'pfmft_snow', 'nfmft_snow', 'nfmft_maxthre',
        ],
        'antarctic_day': [
            'ant4_11l', 'ant4_11h', 'ant4_11m1', 'bt_11_bnds4',
            'pfmft_11maxthre', 'pfmft_btd_min', 'pfmft_snow', 'nfmft_snow', 'nfmft_maxthre',
        ],
        'day_desert': [
            'lds11_12hi', 'lds11_4hi', 'lds11_4lo', 'ldsref2', 'ldsref3', 'ldstci',
            'pfmft_11maxthre', 'pfmft_btd_min', 'pfmft_land', 'nfmft_desert', 'nfmft_maxthre',
        ],
        'polar_day_desert': [
            'pds11_12hi', 'pds11_4hi', 'pds11_4lo', 'pdsref2', 'pdsref3', 'pdstci',
            'pfmft_11maxthre', 'pfmft_btd_min', 'pfmft_snow', 'nfmft_desert', 'nfmft_maxthre',
        ],
        'shadows': ['shadnir', 'shavrat', 'shad124'],
        'noncld_obs': ['nc21', 'ncrat', 'ncvrat', 'ncsig', 'nc11_12'],
        'land_restoral': ['ldsbt11', 'ldsbt11bd', 'ldsr5_4_thr', 'ldr5_4_thr', 'ld20m22', 'ld22m31', 'lnbt11'],
        'coast': ['swc_ndvi'],
        'snow_mask': ['sm_bt11', 'sm_ndsi', 'sm_ref2', 'sm_ref3', 'sm85_11', 'sm37_11', 'sm37_11hel', 'sm_mnir'],
    }

    result = {}
    for key in scene_keys.get(scene, []):
        if key in all_thr:
            result[key] = all_thr[key]
    return result
