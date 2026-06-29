"""Cloud mask utility functions — Python equivalents of Fortran algorithm primitives.

Provides the confidence S-curve functions, spectral tests, and threshold
loading used by the daytime decision tree modules.
"""

import os
from typing import Dict, List, Optional

import numpy as np
import yaml


# ---------------------------------------------------------------------------
# Confidence functions (equivalent to Fortran conf_test / conf_test_2val)
# ---------------------------------------------------------------------------

def conf_test(
    val: np.ndarray,
    locut: float,
    hicut: float,
    power: float,
    midpt: float,
    nmval: int,
) -> np.ndarray:
    """S-curve confidence mapping (Fortran conf_test).

    Maps a continuous value to confidence [0, 1] using an S-curve:
      - val outside [locut, hicut]: confidence = 0 or 1
      - val inside: sigmoid-like interpolation

    Args:
        val: Input values array.
        locut: Low cutoff (confidence=0 below this).
        hicut: High cutoff (confidence=1 above this for nmval=1).
        power: S-curve exponent (1.0 = linear).
        midpt: Midpoint where confidence = 0.5.
        nmval: 1 = monotonically increasing, -1 = decreasing.
               For nmval=1: higher val = higher confidence (e.g. BT test).
               For nmval=-1: higher val = lower confidence (e.g. reflectance test).

    Returns:
        Confidence array (0-1), same shape as val.
    """
    cld = np.abs(hicut - locut) / 2.0
    avg = (hicut + locut) / 2.0
    cld = np.maximum(cld, 1e-9)

    if nmval == 1:
        sig = (val - avg) / cld
    else:
        sig = (avg - val) / cld

    sf = (sig + np.sign(sig) * np.abs(sig) ** power) / 2.0
    conf = 0.5 + sf

    sig_mid = (midpt - avg) / cld
    sf_mid = (sig_mid + np.sign(sig_mid) * np.abs(sig_mid) ** power) / 2.0

    conf = conf - sf_mid
    conf = np.clip(conf, 0.0, 1.0)

    return conf.astype(np.float32)


def conf_test_2val(
    val: np.ndarray,
    locuta: np.ndarray,
    hicuta: np.ndarray,
    power: float,
    midpta: np.ndarray,
    nmval: int,
) -> np.ndarray:
    """Two-sided S-curve confidence (e.g. ratio test).

    For a band ratio test, values too low OR too high indicate cloud.
    """
    # Test lower bound: value below midpta[0] = cloud
    c1 = conf_test(val, locuta[0], hicuta[0], power, midpta[0], -1)
    # Test upper bound: value above midpta[1] = cloud
    c2 = conf_test(val, locuta[1], hicuta[1], power, midpta[1], 1)
    return np.minimum(c1, c2).astype(np.float32)


# ---------------------------------------------------------------------------
# Spectral tests
# ---------------------------------------------------------------------------

def trispc(btd_11_12: np.ndarray) -> np.ndarray:
    """8-11 um clear-sky BTD regression from 11-12 um BTD.

    Fortran trispc function: estimates expected clear-sky 8-11um BTD
    from 11-12um BTD using HIRS regression coefficients.

    trispc(X) = 2.7681 - 3.729*X + 1.054*X^2 - 0.102*X^3
    where X = BT(11um) - BT(12um) in Kelvin.
    """
    x = np.asarray(btd_11_12, dtype=np.float64)
    result = 2.7681 - 3.729 * x + 1.054 * x**2 - 0.102 * x**3
    return result.astype(np.float32)


# ---------------------------------------------------------------------------
# Threshold loading
# ---------------------------------------------------------------------------

def load_thresholds_yaml(
    yaml_path: Optional[str] = None,
) -> Dict:
    """Load cloud mask thresholds from YAML file.

    Args:
        yaml_path: Path to thresholds_mersi_ii.yaml.
                   Defaults to coeff/thresholds_mersi_ii.yaml relative to project root.

    Returns:
        Dict with scene_name -> {'description': ..., 'thresholds': {...}}.
    """
    if yaml_path is None:
        project_root = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
        yaml_path = os.path.join(project_root, "coeff", "thresholds_mersi_ii.yaml")

    with open(yaml_path) as f:
        data = yaml.safe_load(f)

    return data["scenes"]


def get_thresholds(scene: str, yaml_path: Optional[str] = None) -> Dict[str, List[float]]:
    """Get thresholds for a specific scene type.

    Args:
        scene: Scene type name (e.g. 'ocean_day', 'land_nite').
        yaml_path: Optional path to YAML threshold file.

    Returns:
        Dict of parameter_name -> [values].
    """
    scenes = load_thresholds_yaml(yaml_path)
    result = {}
    if scene in scenes:
        result.update(scenes[scene].get("thresholds", {}))
    # Merge shared thresholds
    if "shared" in scenes:
        result.update(scenes["shared"].get("thresholds", {}))
    return result
