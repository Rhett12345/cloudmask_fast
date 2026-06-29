import os
import tempfile

import h5py
import numpy as np

from fylat.mersi_io import io_backend_name, read_geo


def _angle_dataset(group, name, data, slope=0.01, intercept=0.0, fill=-32767):
    ds = group.create_dataset(name, data=np.asarray(data, dtype=np.int16))
    ds.attrs["Slope"] = slope
    ds.attrs["Intercept"] = intercept
    ds.attrs["FillValue"] = fill


def main() -> None:
    assert io_backend_name() == "cpp"

    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "geo.h5")
        with h5py.File(path, "w") as h5:
            geo = h5.create_group("Geolocation")
            geo.create_dataset("Latitude", data=np.array([[10.0, 11.0]], dtype=np.float32))
            geo.create_dataset("Longitude", data=np.array([[100.0, 101.0]], dtype=np.float32))
            _angle_dataset(geo, "SolarZenith", [[3000, 4000]])
            _angle_dataset(geo, "SolarAzimuth", [[1000, 1100]])
            _angle_dataset(geo, "SensorZenith", [[2000, 2100]])
            _angle_dataset(geo, "SensorAzimuth", [[500, 600]])
            geo.create_dataset("DEM", data=np.array([[1.0, 2.0]], dtype=np.float32))
            geo.create_dataset("LandSeaMask", data=np.array([[1, 2]], dtype=np.uint8))

        data = read_geo(path)
        np.testing.assert_allclose(data["lat"], [[10.0, 11.0]])
        np.testing.assert_allclose(data["sza"], [[30.0, 40.0]])
        np.testing.assert_array_equal(data["lsm"], [[1, 2]])
        assert data["nlines"] == 1
        assert data["npixels"] == 2


if __name__ == "__main__":
    main()
