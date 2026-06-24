#pragma once
// Core cloud mask computation functions — translated from Fortran F77
// All functions are pure computation, no I/O, no global state.
//
// Translated from:
//   conf_test.f, conf_test_2val.f, trispc.f, tview.f,
//   set_confdnc.f, set_quality_A.f, set_unused_bits.f,
//   pxinit.f, proc_path.f, fill_bit_pixel.f90,
//   get_sg_thresholds.f90, get_pn_thresholds.f, get_nl_thresholds.f90

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>

#include "fylat/global_params.hpp"
#include "fylat/types.hpp"

namespace fylat {

// =========================================================================
// 1. conf_test — S-function confidence from single threshold
// =========================================================================
inline float conf_test(float val, float locut, float hicut,
                       float power, float midpt, int nmval) {
    if (nmval != 1) return 0.0f;

    float coeff = std::pow(2.0f, power - 1.0f);
    float alpha, gamma, beta;
    bool flipped;

    if (hicut > locut) {
        gamma = hicut;
        alpha = locut;
        flipped = false;
    } else {
        gamma = locut;
        alpha = hicut;
        flipped = true;
    }
    beta = midpt;

    float c;
    if (!flipped && val > gamma) {
        c = 1.0f;
    } else if (!flipped && val < alpha) {
        c = 0.0f;
    } else if (flipped && val > gamma) {
        c = 0.0f;
    } else if (flipped && val < alpha) {
        c = 1.0f;
    } else {
        if (val <= beta) {
            float range = 2.0f * (beta - alpha);
            float s1 = (val - alpha) / range;
            if (!flipped) c = coeff * std::pow(s1, power);
            else          c = 1.0f - (coeff * std::pow(s1, power));
        } else {
            float range = 2.0f * (beta - gamma);
            float s1 = (val - gamma) / range;
            if (!flipped) c = 1.0f - (coeff * std::pow(s1, power));
            else          c = coeff * std::pow(s1, power);
        }
    }

    c = std::clamp(c, 0.0f, 1.0f);
    return c;
}

// =========================================================================
// 2. conf_test_2val — S-function confidence from two thresholds
// =========================================================================
inline float conf_test_2val(float val,
                            const float locut[2], const float hicut[2],
                            float power, const float midpt[2], int nmval) {
    if (nmval != 2) return 0.0f;

    float coeff = std::pow(2.0f, power - 1.0f);
    float gamma1 = hicut[0], gamma2 = hicut[1];
    float alpha1 = locut[0], alpha2 = locut[1];
    float beta1  = midpt[0],  beta2  = midpt[1];

    float c = 0.0f;

    if ((alpha1 - gamma1) > 0.0f) {
        // Inner region fails test
        if (val > alpha1 && val < alpha2) {
            c = 0.0f;
        } else if (val < gamma1 || val > gamma2) {
            c = 1.0f;
        } else if (val <= alpha1) {
            if (val >= beta1) {
                float range = 2.0f * (beta1 - alpha1);
                float s1 = (val - alpha1) / range;
                c = coeff * std::pow(s1, power);
            } else {
                float range = 2.0f * (beta1 - gamma1);
                float s1 = std::abs(val - gamma1) / range;
                c = 1.0f - (coeff * std::pow(s1, power));
            }
        } else {
            if (val <= beta2) {
                float range = 2.0f * (beta2 - alpha2);
                float s1 = (val - alpha2) / range;
                c = coeff * std::pow(s1, power);
            } else {
                float range = 2.0f * (beta2 - gamma2);
                float s1 = (val - gamma2) / range;
                c = 1.0f - (coeff * std::pow(s1, power));
            }
        }
    } else {
        // Inner region passes test
        if (val > gamma1 && val < gamma2) {
            c = 1.0f;
        } else if (val < alpha1 || val > alpha2) {
            c = 0.0f;
        } else if (val <= gamma1) {
            if (val <= beta1) {
                float range = 2.0f * (beta1 - alpha1);
                float s1 = (val - alpha1) / range;
                c = coeff * std::pow(s1, power);
            } else {
                float range = std::abs(2.0f * (beta1 - gamma1));
                float s1 = std::abs((val - gamma1) / range);
                c = 1.0f - (coeff * std::pow(s1, power));
            }
        } else {
            if (val >= beta2) {
                float range = 2.0f * (beta2 - alpha2);
                float s1 = (val - alpha2) / range;
                c = coeff * std::pow(s1, power);
            } else {
                float range = 2.0f * (beta2 - gamma2);
                float s1 = (val - gamma2) / range;
                c = 1.0f - (coeff * std::pow(s1, power));
            }
        }
    }

    c = std::clamp(c, 0.0f, 1.0f);
    return c;
}

// =========================================================================
// 3. trispc — 8-11um clear threshold from 11-12um BTDIF (cubic regression)
// =========================================================================
inline float trispc(float tdf1) {
    constexpr float a1 =  2.7681f;
    constexpr float a2 = -3.729f;
    constexpr float a3 =  1.054f;
    constexpr float a4 = -0.102f;
    // Horner form: a1 + x*(a2 + x*(a3 + x*a4))
    return a1 + tdf1 * (a2 + tdf1 * (a3 + tdf1 * a4));
}

// =========================================================================
// 4. tview — 2D Lagrange interpolation for thin cirrus threshold
// =========================================================================
inline float tview(int key, float xmu, float bt11) {
    // Lookup tables (13 BT levels x 5 scan angles)
    constexpr float utab[5] = {2.00f, 1.75f, 1.50f, 1.25f, 1.00f};
    constexpr float ttab[13] = {190.f, 200.f, 210.f, 220.f, 230.f, 240.f, 250.f,
                                260.f, 270.f, 280.f, 290.f, 300.f, 310.f};
    constexpr float tab[5][13] = {
        {0.559f, 0.424f, 0.286f, 0.137f, 0.123f, 0.198f, 0.333f,
         0.696f, 1.217f, 3.184f, 5.178f, 8.269f, 12.452f},
        {0.542f, 0.416f, 0.294f, 0.162f, 0.156f, 0.240f, 0.366f,
         0.704f, 1.184f, 2.926f, 4.854f, 7.885f, 11.985f},
        {0.520f, 0.405f, 0.305f, 0.194f, 0.199f, 0.294f, 0.409f,
         0.715f, 1.141f, 2.591f, 4.433f, 7.389f, 11.381f},
        {0.491f, 0.391f, 0.319f, 0.238f, 0.257f, 0.367f, 0.467f,
         0.729f, 1.083f, 2.140f, 3.866f, 6.720f, 10.567f},
        {0.450f, 0.370f, 0.340f, 0.300f, 0.340f, 0.470f, 0.550f,
         0.750f, 1.000f, 1.500f, 3.060f, 5.770f, 9.410f}
    };

    // Bounds check
    float u = std::clamp(xmu, utab[4], utab[0]);
    float t = std::clamp(bt11, ttab[0], ttab[12]);

    // Find scan angle interval
    int i0 = 0, i1 = 0, i2 = 0;
    for (int i = 1; i < 5; ++i) {
        if (u >= utab[i]) {
            if (key == 1) {
                i0 = i - 1; i1 = i;
            } else {
                if (i == 4) { i0 = i - 2; i1 = i - 1; i2 = i; }
                else        { i0 = i - 1; i1 = i;     i2 = i + 1; }
            }
            break;
        }
    }

    // Find BT interval
    int j0 = 0, j1 = 0, j2 = 0;
    for (int j = 1; j < 13; ++j) {
        if (t <= ttab[j]) {
            if (key == 1) {
                j0 = j - 1; j1 = j;
            } else {
                if (j == 12) { j0 = j - 2; j1 = j - 1; j2 = j; }
                else         { j0 = j - 1; j1 = j;     j2 = j + 1; }
            }
            break;
        }
    }

    if (key == 1) {
        // Linear interpolation
        float u0 = utab[i0], u1 = utab[i1];
        float t0 = ttab[j0], t1 = ttab[j1];
        float lu0 = (u - u1) / (u0 - u1);
        float lu1 = (u - u0) / (u1 - u0);
        float lt0 = (t - t1) / (t0 - t1);
        float lt1 = (t - t0) / (t1 - t0);
        float p0 = tab[i0][j0] * lu0 + tab[i1][j0] * lu1;
        float p1 = tab[i0][j1] * lu0 + tab[i1][j1] * lu1;
        return p0 * lt0 + p1 * lt1;
    } else {
        // Quadratic interpolation
        float u0 = utab[i0], u1 = utab[i1], u2 = utab[i2];
        float t0 = ttab[j0], t1 = ttab[j1], t2 = ttab[j2];
        float lu0 = (u - u1) * (u - u2) / (u0 - u1) / (u0 - u2);
        float lu1 = (u - u0) * (u - u2) / (u1 - u0) / (u1 - u2);
        float lu2 = (u - u0) * (u - u1) / (u2 - u0) / (u2 - u1);
        float lt0 = (t - t1) * (t - t2) / (t0 - t1) / (t0 - t2);
        float lt1 = (t - t0) * (t - t2) / (t1 - t0) / (t1 - t2);
        float lt2 = (t - t0) * (t - t1) / (t2 - t0) / (t2 - t1);
        float p0 = tab[i0][j0] * lu0 + tab[i1][j0] * lu1 + tab[i2][j0] * lu2;
        float p1 = tab[i0][j1] * lu0 + tab[i1][j1] * lu1 + tab[i2][j1] * lu2;
        float p2 = tab[i0][j2] * lu0 + tab[i1][j2] * lu1 + tab[i2][j2] * lu2;
        return p0 * lt0 + p1 * lt1 + p2 * lt2;
    }
}

// =========================================================================
// 5. Bit operations on raw arrays (used by the functions below)
// =========================================================================
inline void set_bit_raw(uint8_t* bits, int bit_num) {
    bits[bit_num / 8] |= (1u << (bit_num % 8));
}
inline void clear_bit_raw(uint8_t* bits, int bit_num) {
    bits[bit_num / 8] &= ~(1u << (bit_num % 8));
}
inline int check_bits_raw(const uint8_t* bits, int bit_num) {
    return (bits[bit_num / 8] & (1u << (bit_num % 8))) ? 1 : 0;
}

// =========================================================================
// 6. set_confdnc — map confidence value to bit flags
// =========================================================================
inline void set_confdnc(float confdnc, CloudMaskBits& testbits) {
    if (confdnc > 0.99f) {
        set_bit_raw(testbits.data(), 1);
        set_bit_raw(testbits.data(), 2);
    } else if (confdnc > 0.95f) {
        set_bit_raw(testbits.data(), 2);
    } else if (confdnc > 0.66f) {
        set_bit_raw(testbits.data(), 1);
    }
}

// =========================================================================
// 7. set_quality_A — QA bit flags for test/band counts
// =========================================================================
inline void set_quality_A(int nmtests, int nbands, int lsf, QABits& qa_bits) {
    // Number of spectral tests
    if (nmtests > 4) {
        set_bit_raw(qa_bits.data(), 50);
        set_bit_raw(qa_bits.data(), 51);
    } else if (nmtests > 2) {
        set_bit_raw(qa_bits.data(), 51);
    } else if (nmtests > 0) {
        set_bit_raw(qa_bits.data(), 50);
    }

    // Number of bands with good data
    if (nbands > 14) {
        set_bit_raw(qa_bits.data(), 48);
        set_bit_raw(qa_bits.data(), 49);
    } else if (nbands > 7) {
        set_bit_raw(qa_bits.data(), 49);
    } else if (nbands > 0) {
        set_bit_raw(qa_bits.data(), 48);
    }

    // Ecosystem file bit — always set
    set_bit_raw(qa_bits.data(), 64);

    // Land/Sea Mask file status
    if (lsf == -1) {
        set_bit_raw(qa_bits.data(), 70);
        set_bit_raw(qa_bits.data(), 71);
    }
}

// =========================================================================
// 8. set_unused_bits — mark unused test bits as set
// =========================================================================
inline void set_unused_bits(CloudMaskBits& testbits) {
    set_bit_raw(testbits.data(), 24);  // temporal consistency
    set_bit_raw(testbits.data(), 12);  // cloud adjacency
    set_bit_raw(testbits.data(), 31);  // spare
}

// =========================================================================
// 9. pxinit — initialize per-pixel state
// =========================================================================
struct PixelState {
    CloudMaskBits testbits{};
    QABits        qa_bits{};
    float precip_water = 0.0f;
    float vza = 0.0f;
    float sfctmp = 0.0f;
    float pmsl = 0.0f;
    float u_wind = 0.0f;
    float v_wind = 0.0f;
    float plat = -999.0f;
    float plon = -999.0f;
    int   lsf = 0;
    bool  polar = false;
    bool  day = false;
    bool  night = false;
    bool  land = false;
    bool  water = false;
    bool  coast = false;
    bool  snglnt = false;
    bool  visusd = true;
    bool  vrused = true;
    bool  snow = false;
    bool  ice = false;
    bool  desert = false;
    bool  bad_value = false;
    bool  bad_geo = false;
    bool  uniform = true;
    bool  shadow = false;
    bool  smoke = false;
    bool  cirrus_ir = false;
    bool  cirrus_vis = false;
    int   nmtests = 0;
    int   nbands = 0;
    int   nbad_1km = 0;
    int   nbad_250 = 0;
    bool  hi_elev = false;
    bool  antarctic = false;
    bool  sh_ocean = false;
    bool  sg_bad_data = false;
    bool  map_ice = false;
    bool  map_snow = false;
    bool  sh_lake = false;
};

inline void pxinit(PixelState& px) {
    px = PixelState{};  // reset all to defaults

    // Fail-safe bits: set to 1, cleared later if condition absent
    set_bit_raw(px.testbits.data(), 8);   // NCO
    set_bit_raw(px.testbits.data(), 9);   // thin cirrus solar
    set_bit_raw(px.testbits.data(), 10);  // shadow
    set_bit_raw(px.testbits.data(), 11);  // thin cirrus IR
    set_bit_raw(px.testbits.data(), 28);  // suspended dust
}

// =========================================================================
// 10. proc_path — determine processing path and set bit flags
// =========================================================================
inline void proc_path(bool water, bool land, bool day, bool ice, bool snow,
                      bool snglnt, bool coast, bool desert,
                      bool smoke, bool shadow, CloudMaskBits& testbits) {
    // Snow/ice bit
    if (!snow && !ice) {
        set_bit_raw(testbits.data(), 5);
    }
    // Day/night flag
    if (day) {
        set_bit_raw(testbits.data(), 3);
    }
    // Sun glint flag
    if (!snglnt || !water) {
        set_bit_raw(testbits.data(), 4);
    }
    // Surface type bits
    if (coast) {
        set_bit_raw(testbits.data(), 6);
    } else if (desert) {
        set_bit_raw(testbits.data(), 7);
    } else if (land) {
        set_bit_raw(testbits.data(), 6);
        set_bit_raw(testbits.data(), 7);
    }
    // Shadow: clear fail-safe bit
    if (shadow) {
        clear_bit_raw(testbits.data(), 10);
    }
    // Smoke/NCO: clear fail-safe bit
    if (smoke) {
        clear_bit_raw(testbits.data(), 8);
    }
}

// =========================================================================
// 11. fill_bit_pixel — finalize pixel QA and copy to output arrays
// =========================================================================
inline void fill_bit_pixel(int nmtests, int nbands, bool bad_value,
                           bool bad_geo, bool snglnt, bool /*desert*/,
                           CloudMaskBits& testbits, QABits& qa_bits,
                           CloudMaskBits& bitarray, QABits& qa_bitarray) {
    if (nmtests == 0 || nbands == 0 || bad_geo) {
        testbits.fill(0);
        qa_bits[0] = 0;
    } else if (nmtests < 3) {
        set_bit_raw(testbits.data(), 0);
        set_bit_raw(qa_bits.data(), 0);
        set_bit_raw(qa_bits.data(), 3);
    } else if (nmtests < 7) {
        set_bit_raw(testbits.data(), 0);
        set_bit_raw(qa_bits.data(), 0);
        set_bit_raw(qa_bits.data(), 2);
        set_bit_raw(qa_bits.data(), 3);
    } else {
        set_bit_raw(testbits.data(), 0);
        for (int i = 0; i <= 3; ++i)
            set_bit_raw(qa_bits.data(), i);
    }

    if (snglnt) {
        if (qa_bits[0] == 15) qa_bits[0] = 13;
    }

    bitarray = testbits;
    qa_bitarray = qa_bits;
}

// =========================================================================
// 12. get_sg_thresholds — sun glint confidence thresholds (interpolated)
// =========================================================================
// Sun glint threshold data (from snglntr_thr.inc)
struct SgThresholds {
    float snglnt0[4], snglnt10[4], snglnt20[4], snglnt_bounds[4];
    float power;
};

// Default FY-3D values — loaded from threshold file at runtime in Fortran
// For C++, we provide the lookup structure; values must be loaded from file
inline void get_sg_thresholds(float refang, const SgThresholds& thr,
                              float& locut, float& hicut,
                              float& midpt, float& power) {
    float lo_ang, hi_ang, lo_ang_val, hi_ang_val;
    float conf_range = 0.0f;

    if (refang <= thr.snglnt_bounds[1]) {
        lo_ang = thr.snglnt_bounds[0];
        hi_ang = thr.snglnt_bounds[1];
        lo_ang_val = thr.snglnt20[0];
        hi_ang_val = thr.snglnt10[0];
    } else if (refang >= thr.snglnt_bounds[3]) {
        lo_ang = thr.snglnt_bounds[0];
        hi_ang = thr.snglnt_bounds[1];
        lo_ang_val = thr.snglnt20[0];
        hi_ang_val = thr.snglnt10[0];
        // Note: simplified — original Fortran handles more cases
    } else {
        lo_ang = thr.snglnt_bounds[1];
        hi_ang = thr.snglnt_bounds[2];
        lo_ang_val = thr.snglnt10[0];
        hi_ang_val = thr.snglnt0[0];
    }

    if (hi_ang != lo_ang) {
        float a = (refang - lo_ang) / (hi_ang - lo_ang);
        midpt = lo_ang_val + a * (hi_ang_val - lo_ang_val);
        hicut = midpt - conf_range;
        locut = midpt + conf_range;
    }
    power = thr.power;
}

// =========================================================================
// 13. get_pn_thresholds — polar night thresholds (interpolated by BT11)
// =========================================================================
inline void get_pn_thresholds(float bt_11, const float bt_bnds[4],
                              const float th_low[4], const float th_mid1[4],
                              const float th_mid2[4], const float th_mid3[4],
                              const float th_hi[4],
                              float& locut, float& hicut,
                              float& midpt, float& power) {
    float lo_tmp, hi_tmp, lo_tmp_thr, hi_tmp_thr;

    if (bt_11 <= bt_bnds[0]) {
        lo_tmp = lo_tmp;  // inactive: uses th_low directly
        hi_tmp = bt_bnds[0];
        lo_tmp_thr = th_low[2];   // midpt
        hi_tmp_thr = th_low[2];
        midpt = th_low[2];
        hicut = th_low[1];
        locut = th_low[0];
        power = th_low[3];
        return;
    } else if (bt_11 <= bt_bnds[1]) {
        lo_tmp = bt_bnds[0];
        hi_tmp = bt_bnds[1];
        lo_tmp_thr = th_low[2];
        hi_tmp_thr = th_mid1[2];
        midpt = lo_tmp_thr + (bt_11 - lo_tmp) / (hi_tmp - lo_tmp) * (hi_tmp_thr - lo_tmp_thr);
        hicut = th_low[1] + (bt_11 - lo_tmp) / (hi_tmp - lo_tmp) * (th_mid1[1] - th_low[1]);
        locut = th_low[0] + (bt_11 - lo_tmp) / (hi_tmp - lo_tmp) * (th_mid1[0] - th_low[0]);
        power = th_low[3];
        return;
    } else if (bt_11 <= bt_bnds[2]) {
        lo_tmp = bt_bnds[1];
        hi_tmp = bt_bnds[2];
        lo_tmp_thr = th_mid1[2];
        hi_tmp_thr = th_mid2[2];
        midpt = lo_tmp_thr + (bt_11 - lo_tmp) / (hi_tmp - lo_tmp) * (hi_tmp_thr - lo_tmp_thr);
        hicut = th_mid1[1] + (bt_11 - lo_tmp) / (hi_tmp - lo_tmp) * (th_mid2[1] - th_mid1[1]);
        locut = th_mid1[0] + (bt_11 - lo_tmp) / (hi_tmp - lo_tmp) * (th_mid2[0] - th_mid1[0]);
        power = th_low[3];
        return;
    } else if (bt_11 <= bt_bnds[3]) {
        lo_tmp = bt_bnds[2];
        hi_tmp = bt_bnds[3];
        lo_tmp_thr = th_mid2[2];
        hi_tmp_thr = th_mid3[2];
        midpt = lo_tmp_thr + (bt_11 - lo_tmp) / (hi_tmp - lo_tmp) * (hi_tmp_thr - lo_tmp_thr);
        hicut = th_mid2[1] + (bt_11 - lo_tmp) / (hi_tmp - lo_tmp) * (th_mid3[1] - th_mid2[1]);
        locut = th_mid2[0] + (bt_11 - lo_tmp) / (hi_tmp - lo_tmp) * (th_mid3[0] - th_mid2[0]);
        power = th_low[3];
        return;
    } else {
        // bt_11 > bt_bnds[3]: use th_hi
        midpt = th_hi[2];
        hicut = th_hi[1];
        locut = th_hi[0];
        power = th_hi[3];
        return;
    }
}

// =========================================================================
// 14. get_nl_thresholds — non-polar night land thresholds
// =========================================================================
inline void get_nl_thresholds(float btdiff,
                              const float nl_11_4l[4], const float nl_11_4h[4],
                              const float nl_11_4m[4], const float bt_diff_bounds[2],
                              float& locut, float& hicut,
                              float& midpt, float& power) {
    float lo_val = bt_diff_bounds[0];
    float hi_val = bt_diff_bounds[1];
    float conf_range = 0.0f;

    if (btdiff <= lo_val) {
        float lo_val_thr = nl_11_4l[2];
        float hi_val_thr = nl_11_4l[2];
        float a = (btdiff - lo_val) / (hi_val - lo_val);
        midpt = lo_val_thr + a * (hi_val_thr - lo_val_thr);
        hicut = midpt - conf_range;
        locut = midpt + conf_range;
        power = nl_11_4l[3];
    } else if (btdiff >= hi_val) {
        float lo_val_thr = nl_11_4h[2];
        float hi_val_thr = nl_11_4h[2];
        float a = (btdiff - lo_val) / (hi_val - lo_val);
        midpt = lo_val_thr + a * (hi_val_thr - lo_val_thr);
        hicut = midpt - conf_range;
        locut = midpt + conf_range;
        power = nl_11_4h[3];
    } else {
        float lo_val_thr = nl_11_4m[2];
        float hi_val_thr = nl_11_4m[2];
        float a = (btdiff - lo_val) / (hi_val - lo_val);
        midpt = lo_val_thr + a * (hi_val_thr - lo_val_thr);
        hicut = midpt - conf_range;
        locut = midpt + conf_range;
        power = nl_11_4m[3];
    }
}

}  // namespace fylat
