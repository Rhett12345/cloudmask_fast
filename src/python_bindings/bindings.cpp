#include "io/hdf5_io.hpp"

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <cstdint>
#include <map>
#include <string>
#include <vector>

namespace py = pybind11;

namespace {

std::vector<py::ssize_t> to_py_shape(const std::vector<hsize_t>& shape) {
    std::vector<py::ssize_t> result;
    result.reserve(shape.size());
    for (hsize_t dim : shape) {
        result.push_back(static_cast<py::ssize_t>(dim));
    }
    return result;
}

template <typename T>
py::array_t<T> vector_to_array(std::vector<T>&& values,
                               const std::vector<hsize_t>& h5_shape) {
    auto shape = to_py_shape(h5_shape);
    auto* heap_values = new std::vector<T>(std::move(values));
    py::capsule owner(heap_values, [](void* ptr) {
        delete reinterpret_cast<std::vector<T>*>(ptr);
    });
    return py::array_t<T>(shape, heap_values->data(), owner);
}

py::dict float_dataset_map_to_dict(
    std::map<std::string, fylat::io::Float32Dataset>&& datasets) {
    py::dict result;
    for (auto& [name, dataset] : datasets) {
        result[py::str(name)] =
            vector_to_array<float>(std::move(dataset.values), dataset.shape);
    }
    return result;
}

py::dict uint8_dataset_map_to_dict(
    std::map<std::string, fylat::io::Uint8Dataset>&& datasets) {
    py::dict result;
    for (auto& [name, dataset] : datasets) {
        result[py::str(name)] =
            vector_to_array<std::uint8_t>(std::move(dataset.values), dataset.shape);
    }
    return result;
}

}  // namespace

PYBIND11_MODULE(fylat_py, m) {
    m.doc() = "FYLAT C++ HDF5 IO bindings";

    m.def("dataset_shape", &fylat::io::dataset_shape,
          py::arg("file_path"), py::arg("dataset_path"));

    m.def("read_float32", [](const std::string& file_path,
                             const std::string& dataset_path) {
        std::vector<hsize_t> shape;
        auto values = fylat::io::read_float32_dataset(file_path, dataset_path, &shape);
        return vector_to_array<float>(std::move(values), shape);
    }, py::arg("file_path"), py::arg("dataset_path"));

    m.def("read_uint8", [](const std::string& file_path,
                           const std::string& dataset_path) {
        std::vector<hsize_t> shape;
        auto values = fylat::io::read_uint8_dataset(file_path, dataset_path, &shape);
        return vector_to_array<std::uint8_t>(std::move(values), shape);
    }, py::arg("file_path"), py::arg("dataset_path"));

    m.def("write_uint8", [](const std::string& file_path,
                            const std::string& dataset_path,
                            py::array_t<std::uint8_t, py::array::c_style | py::array::forcecast> array,
                            bool overwrite) {
        py::buffer_info info = array.request();
        std::vector<hsize_t> shape;
        shape.reserve(info.shape.size());
        for (py::ssize_t dim : info.shape) {
            shape.push_back(static_cast<hsize_t>(dim));
        }

        auto* ptr = static_cast<std::uint8_t*>(info.ptr);
        std::vector<std::uint8_t> values(ptr, ptr + array.size());
        fylat::io::write_uint8_dataset(file_path, dataset_path, shape, values,
                                       overwrite);
    }, py::arg("file_path"), py::arg("dataset_path"), py::arg("array"),
       py::arg("overwrite") = true);

    m.def("read_mersi_geo", [](const std::string& file_path) {
        py::dict result = float_dataset_map_to_dict(
            fylat::io::read_mersi_geo_float_fields(file_path));
        py::dict byte_fields = uint8_dataset_map_to_dict(
            fylat::io::read_mersi_geo_uint8_fields(file_path));
        for (auto item : byte_fields) {
            result[item.first] = item.second;
        }
        return result;
    }, py::arg("file_path"));

    m.def("read_mersi_l1_payload", [](const std::string& file_path) {
        return float_dataset_map_to_dict(fylat::io::read_mersi_l1_payload(file_path));
    }, py::arg("file_path"));
}
