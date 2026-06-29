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


if __name__ == "__main__":
    main()
