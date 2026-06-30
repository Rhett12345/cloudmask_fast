#include "io/hdf5_io.hpp"

#include <H5Cpp.h>

#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>

namespace fylat::io {
namespace {

std::size_t element_count(const std::vector<hsize_t>& shape) {
    if (shape.empty()) {
        return 1;
    }
    return std::accumulate(shape.begin(), shape.end(), static_cast<std::size_t>(1),
                           [](std::size_t acc, hsize_t dim) {
                               return acc * static_cast<std::size_t>(dim);
                           });
}

void create_groups_for_path(H5::H5File& file, const std::string& dataset_path) {
    std::size_t pos = dataset_path.find('/', dataset_path[0] == '/' ? 1 : 0);
    while (pos != std::string::npos) {
        std::string group_path = dataset_path.substr(0, pos);
        if (!group_path.empty()) {
            try {
                H5::Group group = file.openGroup(group_path);
            } catch (const H5::Exception&) {
                file.createGroup(group_path);
            }
        }
        pos = dataset_path.find('/', pos + 1);
    }
}

template <typename T>
std::vector<T> read_dataset(const std::string& file_path,
                            const std::string& dataset_path,
                            const H5::PredType& mem_type,
                            std::vector<hsize_t>* shape_out) {
    H5::Exception::dontPrint();
    H5::H5File file(file_path, H5F_ACC_RDONLY);
    H5::DataSet dataset = file.openDataSet(dataset_path);
    H5::DataSpace space = dataset.getSpace();

    int rank = space.getSimpleExtentNdims();
    std::vector<hsize_t> shape(static_cast<std::size_t>(rank));
    if (rank > 0) {
        space.getSimpleExtentDims(shape.data(), nullptr);
    }

    std::vector<T> values(element_count(shape));
    dataset.read(values.data(), mem_type);

    if (shape_out != nullptr) {
        *shape_out = shape;
    }
    return values;
}

Float32Dataset read_float32_named(const std::string& file_path,
                                  const std::string& dataset_path) {
    Float32Dataset result;
    result.values = read_float32_dataset(file_path, dataset_path, &result.shape);
    return result;
}

Uint8Dataset read_uint8_named(const std::string& file_path,
                              const std::string& dataset_path) {
    Uint8Dataset result;
    result.values = read_uint8_dataset(file_path, dataset_path, &result.shape);
    return result;
}

}  // namespace

std::vector<hsize_t> dataset_shape(const std::string& file_path,
                                   const std::string& dataset_path) {
    H5::Exception::dontPrint();
    H5::H5File file(file_path, H5F_ACC_RDONLY);
    H5::DataSet dataset = file.openDataSet(dataset_path);
    H5::DataSpace space = dataset.getSpace();

    int rank = space.getSimpleExtentNdims();
    std::vector<hsize_t> shape(static_cast<std::size_t>(rank));
    if (rank > 0) {
        space.getSimpleExtentDims(shape.data(), nullptr);
    }
    return shape;
}

std::vector<float> read_float32_dataset(const std::string& file_path,
                                        const std::string& dataset_path,
                                        std::vector<hsize_t>* shape_out) {
    return read_dataset<float>(file_path, dataset_path, H5::PredType::NATIVE_FLOAT,
                               shape_out);
}

std::vector<std::uint8_t> read_uint8_dataset(const std::string& file_path,
                                             const std::string& dataset_path,
                                             std::vector<hsize_t>* shape_out) {
    return read_dataset<std::uint8_t>(file_path, dataset_path,
                                      H5::PredType::NATIVE_UINT8, shape_out);
}

void write_uint8_dataset(const std::string& file_path,
                         const std::string& dataset_path,
                         const std::vector<hsize_t>& shape,
                         const std::vector<std::uint8_t>& values,
                         bool overwrite) {
    H5::Exception::dontPrint();
    if (shape.empty()) {
        throw std::invalid_argument("write_uint8_dataset requires a non-empty shape");
    }
    if (element_count(shape) != values.size()) {
        throw std::invalid_argument("shape does not match number of values");
    }

    H5::H5File file(file_path, H5F_ACC_RDWR);
    create_groups_for_path(file, dataset_path);

    if (H5Lexists(file.getId(), dataset_path.c_str(), H5P_DEFAULT) > 0) {
        if (!overwrite) {
            throw std::runtime_error("dataset already exists: " + dataset_path);
        }
        H5Ldelete(file.getId(), dataset_path.c_str(), H5P_DEFAULT);
    }

    H5::DataSpace space(static_cast<int>(shape.size()), shape.data());
    H5::DataSet dataset = file.createDataSet(dataset_path, H5::PredType::NATIVE_UINT8,
                                             space);
    dataset.write(values.data(), H5::PredType::NATIVE_UINT8);
}

std::map<std::string, Float32Dataset>
read_mersi_geo_float_fields(const std::string& file_path) {
    const std::pair<const char*, const char*> fields[] = {
        {"lat", "Geolocation/Latitude"},
        {"lon", "Geolocation/Longitude"},
        {"sza_raw", "Geolocation/SolarZenith"},
        {"saa_raw", "Geolocation/SolarAzimuth"},
        {"vza_raw", "Geolocation/SensorZenith"},
        {"vaa_raw", "Geolocation/SensorAzimuth"},
        {"dem", "Geolocation/DEM"},
    };

    std::map<std::string, Float32Dataset> result;
    for (const auto& [name, path] : fields) {
        result.emplace(name, read_float32_named(file_path, path));
    }
    return result;
}

std::map<std::string, Uint8Dataset>
read_mersi_geo_uint8_fields(const std::string& file_path) {
    std::map<std::string, Uint8Dataset> result;
    result.emplace("lsm", read_uint8_named(file_path, "Geolocation/LandSeaMask"));
    return result;
}

std::map<std::string, Float32Dataset>
read_mersi_l1_payload(const std::string& file_path) {
    const std::pair<const char*, const char*> fields[] = {
        {"vis_cal_coeff", "Calibration/VIS_Cal_Coeff"},
        {"ir_cal_coeff", "Calibration/IR_Cal_Coeff"},
        {"ev_250_refsb", "Data/EV_250_Aggr.1KM_RefSB"},
        {"ev_1km_refsb", "Data/EV_1KM_RefSB"},
        {"ev_1km_emissive", "Data/EV_1KM_Emissive"},
        {"ev_250_emissive", "Data/EV_250_Aggr.1KM_Emissive"},
    };

    std::map<std::string, Float32Dataset> result;
    for (const auto& [name, path] : fields) {
        result.emplace(name, read_float32_named(file_path, path));
    }
    return result;
}

}  // namespace fylat::io
