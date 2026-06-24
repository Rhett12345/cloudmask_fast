#pragma once
// Physical and mathematical constants — translated from constant.f90

#include <cstdint>

namespace fylat {

// ---------------------------------------------------------------------------
// Status codes
// ---------------------------------------------------------------------------
constexpr int8_t SUCCESS = 0;
constexpr int8_t WARNING = 1;
constexpr int8_t ERROR   = 2;
constexpr int8_t FAILURE = 3;
constexpr int8_t SUCCEED = 0;
constexpr int8_t FAIL    = -1;

// ---------------------------------------------------------------------------
// Missing / sentinel values
// ---------------------------------------------------------------------------
constexpr double      MISSING_VALUE_REAL8 = -999.0;
constexpr float       MISSING_VALUE_REAL4 = -999.0f;
constexpr int8_t      MISSING_VALUE_INT1  = -128;
constexpr int16_t     MISSING_VALUE_INT2  = -32768;
constexpr int32_t     MISSING_VALUE_INT4  = -999;

// ---------------------------------------------------------------------------
// Mathematical constants
// ---------------------------------------------------------------------------
constexpr double R_EARTH      = 6378206.4;
constexpr double PI           = 3.141592653589793238462643;
constexpr double DTOR         = PI / 180.0;
constexpr double RTOD         = 180.0 / PI;
constexpr double SPEED_OF_LIGHT      = 2.99792458e+08;
constexpr double PLANCK_CONSTANT      = 6.62606896e-34;
constexpr double GRAVITATIONAL_CONSTANT = 6.67428e-11;
constexpr double STANDARD_ATMOSPHERE  = 101325.0;
constexpr double STANDARD_TEMPERATURE = 273.15;
constexpr double STANDARD_GRAVITY     = 9.80665;
constexpr double AVOGADRO_CONSTANT    = 6.02214179e+23;
constexpr double MOLAR_GAS_CONSTANT   = 8.314472;
constexpr double BOLTZMANN_CONSTANT   = 1.3806504e-23;
constexpr double STEFAN_BOLTZMANN     = 5.670400e-08;

// Planck constants (W·m²)
constexpr double C_1 = 1.191042722e-16;
constexpr double C_2 = 1.4387752e-02;

}  // namespace fylat
