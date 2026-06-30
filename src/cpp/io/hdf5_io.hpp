#pragma once

#include <H5public.h>

#include <cstdint>
#include <map>
#include <string>
#include <vector>

namespace fylat::io {

struct Float32Dataset {
    std::vector<hsize_t> shape;
    std::vector<float> values;
};

struct Uint8Dataset {
    std::vector<hsize_t> shape;
    std::vector<std::uint8_t> values;
};

std::vector<hsize_t> dataset_shape(const std::string& file_path,
                                   const std::string& dataset_path);

std::vector<float> read_float32_dataset(const std::string& file_path,
                                        const std::string& dataset_path,
                                        std::vector<hsize_t>* shape_out = nullptr);

std::vector<std::uint8_t> read_uint8_dataset(const std::string& file_path,
                                             const std::string& dataset_path,
                                             std::vector<hsize_t>* shape_out = nullptr);

void write_uint8_dataset(const std::string& file_path,
                         const std::string& dataset_path,
                         const std::vector<hsize_t>& shape,
                         const std::vector<std::uint8_t>& values,
                         bool overwrite = true);

std::map<std::string, Float32Dataset>
read_mersi_geo_float_fields(const std::string& file_path);

std::map<std::string, Uint8Dataset>
read_mersi_geo_uint8_fields(const std::string& file_path);

std::map<std::string, Float32Dataset>
read_mersi_l1_payload(const std::string& file_path);

}  // namespace fylat::io
