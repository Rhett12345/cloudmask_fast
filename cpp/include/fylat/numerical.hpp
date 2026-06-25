#pragma once
// Geometry and date utility functions — translated from numerical.f90
// and frontend_module.f90
#include <cmath>
#include <cstdint>
#include "fylat/constants.hpp"

namespace fylat {

// =========================================================================
// Date/time utilities
// =========================================================================

inline int leap_year(int year) {
    return ((year % 4 == 0 && year % 100 != 0) || year % 400 == 0) ? 1 : 0;
}

inline int compute_daynum(int month, int day, int ileap) {
    constexpr int cum_days[12] = {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334};
    return cum_days[month - 1] + day + (month > 2 ? ileap : 0);
}

inline double julian(int year, int month, int day, int hour, int minute) {
    int iy = year, im = month;
    if (month <= 2) { iy = year - 1; im = month + 12; }
    double jd = static_cast<double>(
        static_cast<int>(365.25 * (iy + 4716.0)) +
        static_cast<int>(30.6001 * (im + 1.0)) + 2 -
        static_cast<int>(iy / 100.0) +
        static_cast<int>(static_cast<int>(iy / 100.0) / 4.0) +
        day - 1524);
    jd += (hour + minute / 60.0) / 24.0;
    return jd;
}

// =========================================================================
// Geometry: Earth-Sun distance, zenith angles, scattering angle
// =========================================================================

inline float compute_earth2sun(int julday) {
    return static_cast<float>(1.0 - 0.016729 * std::cos(0.9856 * (julday - 4.0) * DTOR));
}

inline void compute_cos_zenith_angles(const float* satzen, const float* solzen,
                                       float* cos_satzen, float* cos_solzen,
                                       int n) {
    for (int i = 0; i < n; ++i) {
        cos_satzen[i] = std::cos(satzen[i] * static_cast<float>(DTOR));
        cos_solzen[i] = std::cos(solzen[i] * static_cast<float>(DTOR));
    }
}

inline float compute_scat_zen(float cos_sol, float cos_sat,
                              float sin_sol, float sin_sat, float cos_relaz) {
    float val = -(cos_sol * cos_sat - sin_sol * sin_sat * cos_relaz);
    if (val > 1.0f) val = 0.0f;
    return std::acos(val) / static_cast<float>(DTOR);
}

inline void compute_scattering_angles(const float* cos_satzen, const float* cos_solzen,
                                       const float* satzen, const float* solzen,
                                       const float* relaz, float* scatzen, int n) {
    for (int i = 0; i < n; ++i) {
        float sin_sol = std::sin(solzen[i] * static_cast<float>(DTOR));
        float sin_sat = std::sin(satzen[i] * static_cast<float>(DTOR));
        float cos_rel = std::cos(relaz[i] * static_cast<float>(DTOR));
        scatzen[i] = compute_scat_zen(cos_solzen[i], cos_satzen[i], sin_sol, sin_sat, cos_rel);
    }
}

// =========================================================================
// Sun glint angle (for ocean processing)
// =========================================================================

inline float compute_glint_angle(float cos_sol, float cos_sat, float scat_zen_deg) {
    float scat_rad = scat_zen_deg * static_cast<float>(DTOR);
    float glint = std::acos(cos_sol * cos_sat + std::sin(std::acos(cos_sol)) *
                            std::sin(std::acos(cos_sat)) * std::cos(scat_rad));
    return glint / static_cast<float>(DTOR);
}

}  // namespace fylat
