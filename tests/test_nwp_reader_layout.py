import os
import sys


PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PROJECT_ROOT, "python"))

from fylat import nwp_reader  # noqa: E402


def main() -> None:
    labels = nwp_reader._build_namelist()
    specs = nwp_reader._build_field_specs()

    assert len(labels) == len(specs)
    assert len(labels) == 283
    assert labels[:5] == [
        "PRES:surface",
        "PRMSL:mean sea level",
        "TMP:surface",
        "HGT:surface",
        "PRES:surface",
    ]
    assert labels[-2:] == [
        "UGRD:10 m above ground",
        "VGRD:10 m above ground",
    ]
    assert nwp_reader.fortran_nwp_binary_name("20220803", 6) == "gfs0p25_41L_20220803_06_00"
    assert nwp_reader.fortran_nwp_binary_name("20220803", 24) == "gfs0p25_41L_20220803_24_00"


if __name__ == "__main__":
    main()
