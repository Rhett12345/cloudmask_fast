#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "fylat/constants.hpp"
#include "fylat/global_params.hpp"

namespace py = pybind11;

namespace fylat {
void set_bit(py::array_t<uint8_t> testbits, int bit_num);
void clear_bit(py::array_t<uint8_t> testbits, int bit_num);
int check_bits(py::array_t<uint8_t> testbits, int bit_num);
void set_qa_bit(py::array_t<uint8_t> qa_bits, int bit_num);
int check_qa_bits(py::array_t<uint8_t> qa_bits, int bit_num);
}  // namespace fylat

PYBIND11_MODULE(fylat_core, m) {
    m.doc() = "FYLAT cloud mask C++ compute kernel";

    // Constants
    m.attr("PI") = fylat::PI;
    m.attr("DTOR") = fylat::DTOR;
    m.attr("RTOD") = fylat::RTOD;
    m.attr("BAD_DATA") = fylat::BAD_DATA;
    m.attr("C_1") = fylat::C_1;
    m.attr("C_2") = fylat::C_2;
    m.attr("CM_BYTE_DIM") = fylat::CM_BYTE_DIM;
    m.attr("CM_QA_DIM") = fylat::CM_QA_DIM;
    m.attr("INBAND") = fylat::INBAND;
    m.attr("IR_BAND") = fylat::IR_BAND;

    // Bit operations
    m.def("set_bit", &fylat::set_bit,
          py::arg("testbits"), py::arg("bit_num"),
          "Set a bit in the 48-bit cloud mask array (modified in-place)");
    m.def("clear_bit", &fylat::clear_bit,
          py::arg("testbits"), py::arg("bit_num"),
          "Clear a bit in the 48-bit cloud mask array (modified in-place)");
    m.def("check_bits", &fylat::check_bits,
          py::arg("testbits"), py::arg("bit_num"),
          "Check if bit is set (returns 0 or 1)");
    m.def("set_qa_bit", &fylat::set_qa_bit,
          py::arg("qa_bits"), py::arg("bit_num"),
          "Set a bit in the 80-bit QA array (modified in-place)");
    m.def("check_qa_bits", &fylat::check_qa_bits,
          py::arg("qa_bits"), py::arg("bit_num"),
          "Check if QA bit is set (returns 0 or 1)");
}
