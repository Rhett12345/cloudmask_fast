#pragma once
// Common type aliases for the FYLAT C++ library

#include <array>
#include <cstdint>

namespace fylat {

// Bit array types matching Fortran dimensions
using CloudMaskBits = std::array<uint8_t, 6>;   // testbits(6)  — 48-bit cloud mask
using QABits        = std::array<uint8_t, 10>;  // qa_bits(10)  — 80-bit QA flags

// 3×3 spatial context arrays (row-major: [line][elem])
template <typename T>
using Context3x3 = std::array<std::array<T, 3>, 3>;

}  // namespace fylat
