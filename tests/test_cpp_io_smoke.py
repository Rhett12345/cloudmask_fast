import os
import tempfile

import h5py
import numpy as np

import fylat_py


def main() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "io_smoke.h5")
        with h5py.File(path, "w") as h5:
            h5.create_dataset("Data/Float32", data=np.arange(12, dtype=np.float32).reshape(3, 4))
            h5.create_dataset("Data/Mask", data=np.arange(6, dtype=np.uint8).reshape(2, 3))
            geo = h5.create_group("Geolocation")
            geo.create_dataset("Latitude", data=np.array([[10.0, 11.0]], dtype=np.float32))
            geo.create_dataset("Longitude", data=np.array([[100.0, 101.0]], dtype=np.float32))
            geo.create_dataset("SolarZenith", data=np.array([[3000, 4000]], dtype=np.int16))
            geo.create_dataset("SolarAzimuth", data=np.array([[1000, 1100]], dtype=np.int16))
            geo.create_dataset("SensorZenith", data=np.array([[2000, 2100]], dtype=np.int16))
            geo.create_dataset("SensorAzimuth", data=np.array([[500, 600]], dtype=np.int16))
            geo.create_dataset("DEM", data=np.array([[1.0, 2.0]], dtype=np.float32))
            geo.create_dataset("LandSeaMask", data=np.array([[1, 2]], dtype=np.uint8))
            cal = h5.create_group("Calibration")
            cal.create_dataset("VIS_Cal_Coeff", data=np.ones((19, 3), dtype=np.float32))
            cal.create_dataset("IR_Cal_Coeff", data=np.ones((6, 4, 2), dtype=np.float32))
            data = h5["Data"]
            data.create_dataset("EV_250_Aggr.1KM_RefSB", data=np.ones((4, 1, 2), dtype=np.float32))
            data.create_dataset("EV_1KM_RefSB", data=np.ones((15, 1, 2), dtype=np.float32))
            data.create_dataset("EV_1KM_Emissive", data=np.ones((4, 1, 2), dtype=np.float32))
            data.create_dataset("EV_250_Aggr.1KM_Emissive", data=np.ones((2, 1, 2), dtype=np.float32))

        assert fylat_py.dataset_shape(path, "Data/Float32") == [3, 4]
        np.testing.assert_allclose(
            fylat_py.read_float32(path, "Data/Float32"),
            np.arange(12, dtype=np.float32).reshape(3, 4),
        )
        np.testing.assert_array_equal(
            fylat_py.read_uint8(path, "Data/Mask"),
            np.arange(6, dtype=np.uint8).reshape(2, 3),
        )

        out = np.array([[3, 2], [1, 0]], dtype=np.uint8)
        fylat_py.write_uint8(path, "Output/Cloud_Mask", out)
        with h5py.File(path, "r") as h5:
            np.testing.assert_array_equal(h5["Output/Cloud_Mask"][:], out)

        geo_payload = fylat_py.read_mersi_geo(path)
        np.testing.assert_allclose(geo_payload["lat"], [[10.0, 11.0]])
        np.testing.assert_array_equal(geo_payload["lsm"], [[1, 2]])

        l1_payload = fylat_py.read_mersi_l1_payload(path)
        assert l1_payload["vis_cal_coeff"].shape == (19, 3)
        assert l1_payload["ev_1km_refsb"].shape == (15, 1, 2)


if __name__ == "__main__":
    main()
