#pragma once

#include <H5public.h>

#include <cstdint>
#include <string>
#include <vector>

namespace fylat::io {

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

}  // namespace fylat::io
