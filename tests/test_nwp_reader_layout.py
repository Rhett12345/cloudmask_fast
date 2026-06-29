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


if __name__ == "__main__":
    main()
