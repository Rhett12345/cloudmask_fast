#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "fylat/cloudmask_core.hpp"
#include "fylat/cloudmask_spatial.hpp"
#include "fylat/constants.hpp"
#include "fylat/global_params.hpp"
#include "fylat/numerical.hpp"

namespace py = pybind11;

// ---------------------------------------------------------------------------
// Forward declarations of numpy-array wrappers (in bitops.cpp)
// ---------------------------------------------------------------------------
namespace fylat {
void set_bit(py::array_t<uint8_t> testbits, int bit_num);
void clear_bit(py::array_t<uint8_t> testbits, int bit_num);
int check_bits(py::array_t<uint8_t> testbits, int bit_num);
void set_qa_bit(py::array_t<uint8_t> qa_bits, int bit_num);
int check_qa_bits(py::array_t<uint8_t> qa_bits, int bit_num);
}

// ---------------------------------------------------------------------------
// Wrappers that take numpy arrays and delegate to raw-bit C++ functions
// ---------------------------------------------------------------------------
static void wrap_set_confdnc(float confdnc, py::array_t<uint8_t> testbits) {
    auto buf = testbits.mutable_unchecked<1>();
    fylat::CloudMaskBits bits;
    std::copy(buf.data(0), buf.data(0) + 6, bits.begin());
    fylat::set_confdnc(confdnc, bits);
    std::copy(bits.begin(), bits.end(), buf.mutable_data(0));
}

static void wrap_set_quality_A(int nmtests, int nbands, int lsf,
                                py::array_t<uint8_t> qa_bits) {
    auto buf = qa_bits.mutable_unchecked<1>();
    fylat::QABits bits;
    std::copy(buf.data(0), buf.data(0) + 10, bits.begin());
    fylat::set_quality_A(nmtests, nbands, lsf, bits);
    std::copy(bits.begin(), bits.end(), buf.mutable_data(0));
}

static void wrap_set_unused_bits(py::array_t<uint8_t> testbits) {
    auto buf = testbits.mutable_unchecked<1>();
    fylat::CloudMaskBits bits;
    std::copy(buf.data(0), buf.data(0) + 6, bits.begin());
    fylat::set_unused_bits(bits);
    std::copy(bits.begin(), bits.end(), buf.mutable_data(0));
}

static void wrap_proc_path(bool water, bool land, bool day, bool ice, bool snow,
                           bool snglnt, bool coast, bool desert,
                           bool smoke, bool shadow,
                           py::array_t<uint8_t> testbits) {
    auto buf = testbits.mutable_unchecked<1>();
    fylat::CloudMaskBits bits;
    std::copy(buf.data(0), buf.data(0) + 6, bits.begin());
    fylat::proc_path(water, land, day, ice, snow, snglnt, coast, desert,
                     smoke, shadow, bits);
    std::copy(bits.begin(), bits.end(), buf.mutable_data(0));
}

// numpy-array wrappers for functions taking float arrays
static float wrap_conf_test_2val(float val,
                                  py::array_t<float> locut,
                                  py::array_t<float> hicut,
                                  float power,
                                  py::array_t<float> midpt,
                                  int nmval) {
    auto lo = locut.unchecked<1>();
    auto hi = hicut.unchecked<1>();
    auto md = midpt.unchecked<1>();
    float lo_arr[2] = {lo(0), lo(1)};
    float hi_arr[2] = {hi(0), hi(1)};
    float md_arr[2] = {md(0), md(1)};
    return fylat::conf_test_2val(val, lo_arr, hi_arr, power, md_arr, nmval);
}

// ===========================================================================
// Module definition
// ===========================================================================
PYBIND11_MODULE(fylat_core, m) {
    m.doc() = "FYLAT cloud mask C++ compute kernel";

    // -- Constants --
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

    // -- Bit operations (zero-copy numpy) --
    m.def("set_bit", &fylat::set_bit, py::arg("testbits"), py::arg("bit_num"));
    m.def("clear_bit", &fylat::clear_bit, py::arg("testbits"), py::arg("bit_num"));
    m.def("check_bits", &fylat::check_bits, py::arg("testbits"), py::arg("bit_num"));
    m.def("set_qa_bit", &fylat::set_qa_bit, py::arg("qa_bits"), py::arg("bit_num"));
    m.def("check_qa_bits", &fylat::check_qa_bits, py::arg("qa_bits"), py::arg("bit_num"));

    // -- Confidence / spectral functions (pure float computation) --
    m.def("conf_test", &fylat::conf_test,
          py::arg("val"), py::arg("locut"), py::arg("hicut"),
          py::arg("power"), py::arg("midpt"), py::arg("nmval"),
          "S-function confidence from single threshold");
    m.def("conf_test_2val", &wrap_conf_test_2val,
          py::arg("val"), py::arg("locut"), py::arg("hicut"),
          py::arg("power"), py::arg("midpt"), py::arg("nmval"),
          "S-function confidence from two thresholds (numpy arrays for locut/hicut/midpt)");
    m.def("trispc", &fylat::trispc, py::arg("tdf1"),
          "8-11um clear threshold from 11-12um BTDIF (cubic regression)");
    m.def("tview", &fylat::tview,
          py::arg("key"), py::arg("xmu"), py::arg("bt11"),
          "2D Lagrange interpolation for thin cirrus threshold");

    // -- Bit-array-based state functions --
    m.def("set_confdnc", &wrap_set_confdnc,
          py::arg("confdnc"), py::arg("testbits"),
          "Set confidence bit flags from confidence value");
    m.def("set_quality_A", &wrap_set_quality_A,
          py::arg("nmtests"), py::arg("nbands"), py::arg("lsf"),
          py::arg("qa_bits"),
          "Set QA bit flags for test/band counts");
    m.def("set_unused_bits", &wrap_set_unused_bits, py::arg("testbits"),
          "Mark unused test bits as set");
    m.def("proc_path", &wrap_proc_path,
          py::arg("water"), py::arg("land"), py::arg("day"),
          py::arg("ice"), py::arg("snow"), py::arg("snglnt"),
          py::arg("coast"), py::arg("desert"), py::arg("smoke"),
          py::arg("shadow"), py::arg("testbits"),
          "Determine processing path and set bit flags");

    // -- Spatial / auxiliary functions (selected testable subset) --
    m.def("spatial_var",
          [](py::array_t<float> diff, py::array_t<float> dovar11) {
              auto d = diff.unchecked<2>();
              auto dv = dovar11.unchecked<1>();
              float d_arr[2][8];
              for (int b = 0; b < 2; ++b)
                  for (int i = 0; i < 8; ++i)
                      d_arr[b][i] = d(b, i);
              float dovar[1] = {dv(0)};
              int ipt, result;
              fylat::spatial_var(d_arr, dovar, ipt, result);
              return std::make_pair(ipt, result);
          },
          py::arg("diff"), py::arg("dovar11"),
          "Spatial uniformity test — returns (ipt, result)");

    m.def("shadows",
          [](py::array_t<float> pxldat, bool visusd,
             float shadnir0, float shadnir1, float shavrat0,
             py::array_t<uint8_t> qa_bits_arr) {
              auto pd = pxldat.unchecked<1>();
              auto qb = qa_bits_arr.mutable_unchecked<1>();
              float px[25];
              for (int i = 0; i < 25; ++i) px[i] = pd(i);
              float sn[2] = {shadnir0, shadnir1};
              float sv[1] = {shavrat0};
              fylat::QABits qa;
              std::copy(qb.data(0), qb.data(0) + 10, qa.begin());
              bool shadow;
              fylat::shadows(px, visusd, sn, sv, shadow, qa);
              std::copy(qa.begin(), qa.end(), qb.mutable_data(0));
              return shadow;
          },
          py::arg("pxldat"), py::arg("visusd"),
          py::arg("shadnir0"), py::arg("shadnir1"), py::arg("shavrat0"),
          py::arg("qa_bits"),
          "Cloud shadow detection");

    // -- Numerical / geometry functions --
    m.def("compute_earth2sun", &fylat::compute_earth2sun, py::arg("julday"),
          "Earth-Sun distance in AU for given day of year");
    m.def("leap_year", &fylat::leap_year, py::arg("year"),
          "Check if year is a leap year (returns 0 or 1)");
    m.def("compute_daynum", &fylat::compute_daynum,
          py::arg("month"), py::arg("day"), py::arg("ileap"),
          "Compute day number within a year");
    m.def("julian", &fylat::julian,
          py::arg("year"), py::arg("month"), py::arg("day"),
          py::arg("hour"), py::arg("minute"),
          "Compute Julian day from date components");
}
