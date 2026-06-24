#include <cstdint>
#include <pybind11/numpy.h>

namespace py = pybind11;

namespace fylat {

// ---------------------------------------------------------------------------
// Bit operations on raw byte arrays (called from pybind11 wrappers)
// ---------------------------------------------------------------------------

inline void set_bit_raw(uint8_t* bits, int bit_num) {
    int ibyte = bit_num / 8;
    int ibit  = bit_num - ibyte * 8;
    bits[ibyte] |= (1u << ibit);
}

inline void clear_bit_raw(uint8_t* bits, int bit_num) {
    int ibyte = bit_num / 8;
    int ibit  = bit_num - ibyte * 8;
    bits[ibyte] &= ~(1u << ibit);
}

inline int check_bits_raw(const uint8_t* bits, int bit_num) {
    int ibyte = bit_num / 8;
    int ibit  = bit_num - ibyte * 8;
    return (bits[ibyte] & (1u << ibit)) ? 1 : 0;
}

// ---------------------------------------------------------------------------
// pybind11-callable wrappers (zero-copy numpy)
// ---------------------------------------------------------------------------

void set_bit(py::array_t<uint8_t> testbits, int bit_num) {
    auto buf = testbits.mutable_unchecked<1>();
    set_bit_raw(buf.mutable_data(0), bit_num);
}

void clear_bit(py::array_t<uint8_t> testbits, int bit_num) {
    auto buf = testbits.mutable_unchecked<1>();
    clear_bit_raw(buf.mutable_data(0), bit_num);
}

int check_bits(py::array_t<uint8_t> testbits, int bit_num) {
    auto buf = testbits.unchecked<1>();
    return check_bits_raw(buf.data(0), bit_num);
}

void set_qa_bit(py::array_t<uint8_t> qa_bits, int bit_num) {
    auto buf = qa_bits.mutable_unchecked<1>();
    set_bit_raw(buf.mutable_data(0), bit_num);
}

int check_qa_bits(py::array_t<uint8_t> qa_bits, int bit_num) {
    auto buf = qa_bits.unchecked<1>();
    return check_bits_raw(buf.data(0), bit_num);
}

}  // namespace fylat
