import os
import tempfile

from fylat import nwp_reader


def _write_one_field(path: str) -> None:
    with open(path, "wb") as f:
        f.write(b"\0" * nwp_reader._FIELD_BYTES)


def main() -> None:
    original_eccodes = nwp_reader._grib2_to_binary_eccodes
    original_wgrib2 = nwp_reader._grib2_to_binary_wgrib2
    original_env = os.environ.get("FYLAT_NWP_BACKEND")
    calls = []

    def fake_eccodes(grib_path, bin_path):
        calls.append(("eccodes", os.path.basename(grib_path), os.path.basename(bin_path)))
        _write_one_field(bin_path)
        return 1, []

    def fake_wgrib2(grib_path, bin_path):
        calls.append(("wgrib2", os.path.basename(grib_path), os.path.basename(bin_path)))
        _write_one_field(bin_path)

    try:
        nwp_reader._grib2_to_binary_eccodes = fake_eccodes
        nwp_reader._grib2_to_binary_wgrib2 = fake_wgrib2

        with tempfile.TemporaryDirectory() as tmp:
            grib = os.path.join(tmp, "gfs.t06z.pgrb2.0p25.f024")
            open(grib, "wb").close()

            os.environ.pop("FYLAT_NWP_BACKEND", None)
            nwp_reader.grib2_to_binary(grib, os.path.join(tmp, "default.bin"))
            assert calls[-1][0] == "eccodes"

            os.environ["FYLAT_NWP_BACKEND"] = "wgrib2"
            nwp_reader.grib2_to_binary(grib, os.path.join(tmp, "wgrib2.bin"))
            assert calls[-1][0] == "wgrib2"

            try:
                nwp_reader.grib2_to_binary(grib, os.path.join(tmp, "bad.bin"), backend="bad")
            except ValueError as exc:
                assert "Unsupported NWP backend" in str(exc)
            else:
                raise AssertionError("unsupported backend did not raise")
    finally:
        nwp_reader._grib2_to_binary_eccodes = original_eccodes
        nwp_reader._grib2_to_binary_wgrib2 = original_wgrib2
        if original_env is None:
            os.environ.pop("FYLAT_NWP_BACKEND", None)
        else:
            os.environ["FYLAT_NWP_BACKEND"] = original_env


if __name__ == "__main__":
    main()
