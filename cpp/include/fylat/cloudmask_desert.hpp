#pragma once
// Desert, coast, and Antarctic scene test modules
// Translated from:
//   LandDay_desert.f90, LandDay_desert_c.f90, LandDay_coast.f90,
//   PolarDay_desert.f90, PolarDay_desert_c.f90, PolarDay_coast.f90,
//   Antarctic_day.f90

#include <algorithm>
#include <cmath>
#include <cstdint>

#include "fylat/cloudmask_core.hpp"
#include "fylat/global_params.hpp"
#include "fylat/types.hpp"

namespace fylat {

// =========================================================================
// Common pattern for all desert/coast day scenes (Groups 1-4):
//   Group 1: pf/nf mft tests
//   Group 2: 11-12 thin cirrus (APOLLO) + 11-4 fog/low cloud
//   Group 3: VIS/NIR reflectance test
//   Group 4: 1.38um NIR high cloud + thin cirrus flag
// Confidence: (cmin1*cmin2*cmin3*cmin4)^(1/groups)
// =========================================================================

// Helper to reduce repetition in desert/coast group 1
inline void desert_group1(const float pxldat[INBAND],
                          const PfmftNfmftThr& pf, int is_cold_sfc,
                          const float btclr[7],
                          CloudMaskBits& testbits, QABits& qa_bits,
                          int& nmtests, int& ngtests, float& cmin1,
                          const float* pfmft_land_override = nullptr,
                          const float* nfmft_land_override = nullptr) {
    float m31 = pxldat[23], m32 = pxldat[24];
    const float* pf_land = pfmft_land_override ? pfmft_land_override : pf.pfmft_land;
    const float* nf_land = nfmft_land_override ? nfmft_land_override : pf.nfmft_land;

    // pfmft
    if (m31 > 0 && m32 > 0 && m31 < pf.pfmft_11maxthre[0] &&
        (btclr[4]-btclr[5]) > pf.pfmft_btd_min[0]) {
        float tv11_12;
        if (m31 > 270.0f && btclr[4] > 270.0f)
            tv11_12 = (m31-m32) - (btclr[4]-btclr[5])*(m31-260.0f)/(btclr[4]-260.0f);
        else
            tv11_12 = m31 - m32;
        set_bit_raw(testbits.data(), 14); set_bit_raw(qa_bits.data(), 14);
        const float* thr = (is_cold_sfc == 1) ? pf.pfmft_cold : pf_land;
        cmin1 = std::min(cmin1, conf_test(tv11_12, thr[0], thr[2], thr[3], thr[1], 1));
        ++ngtests;
    }
    // nfmft
    if (m31 > 0 && m32 > 0 && (m31-m32) <= pf.nfmft_maxthre[0]) {
        float tv11_12 = (m31-m32) - (btclr[4]-btclr[5]);
        set_bit_raw(testbits.data(), 15); set_bit_raw(qa_bits.data(), 15);
        cmin1 = std::min(cmin1, conf_test(tv11_12, nf_land[0], nf_land[2], nf_land[3], nf_land[1], 1));
        ++ngtests;
    }
    if (ngtests > 0) nmtests += ngtests;
}

// Helper for group 2: thin cirrus + 11-4
inline void desert_group2(const float pxldat[INBAND], float vza, bool visusd,
                          float cirrus_hi_thr, const float* lo4_thr,
                          bool use_2val, const float* lo4_hi = nullptr,
                          CloudMaskBits& testbits = CloudMaskBits{},
                          QABits& qa_bits = QABits{},
                          int& nmtests = *(int*)nullptr, int& ngtests = *(int*)nullptr,
                          float& cmin2 = *(float*)nullptr) {
    // This is called with reference semantics — caller provides the refs
    (void)testbits; (void)qa_bits; (void)nmtests; (void)ngtests; (void)cmin2;
    // For the actual implementations, we inline the logic per scene
}

// =========================================================================
// 1. LandDay_desert — daytime desert
// Key difference: uses conf_test_2val for 11-12um (3 thresholds)
// =========================================================================
inline void land_day_desert(const float pxldat[INBAND], float vza,
                            bool visusd, bool& cirrus_vis,
                            const DesertDayThr& thr, const PfmftNfmftThr& pf,
                            int is_cold_sfc,
                            CloudMaskBits& testbits, QABits& qa_bits,
                            int& nmtests, float& confdnc, const float btclr[7]) {
    float m31 = pxldat[23], m32 = pxldat[24], m20 = pxldat[19];
    float b2 = pxldat[1], b5 = pxldat[17];  // 0.87um, 1.38um

    // Group 1
    float cmin1 = 1.0f; int ng1 = 0;
    // pfmft/nfmft (desert uses pf.nfmft_desert)
    const float nf_desert[4] = {pf.nfmft_desert[0], pf.nfmft_desert[1], pf.nfmft_desert[2], pf.nfmft_desert[3]};
    desert_group1(pxldat, pf, is_cold_sfc, btclr, testbits, qa_bits, nmtests, ng1, cmin1, nullptr, nf_desert);

    // Group 2
    float cmin2 = 1.0f; int ng2 = 0;
    if (m31 > 0 && m32 > 0 && vza > 0) {
        float cv = std::cos(vza * DTOR);
        if (std::abs(cv) > 1e-6f) {
            float dt = tview(1, 1.0f/cv, m31);
            float df = (dt >= 0.1f) ? dt : thr.lds11_12hi[0];
            if (m31-m32 <= df) { set_bit_raw(testbits.data(), 18); set_bit_raw(qa_bits.data(), 18); }
            float lo[2] = {df + 0.3f*df, thr.lds11_12hi[1]};
            float hi[2] = {df - 0.3f*df, thr.lds11_12hi[2]};
            float md[2] = {df, df};
            cmin2 = std::min(cmin2, conf_test_2val(m31-m32, lo, hi, 1.0f, md, 2));
            ++ng2; ++nmtests;
        }
    }
    if (visusd && m31 > 0 && m20 > 0) {
        float mas11_4 = m31 - m20;
        set_bit_raw(qa_bits.data(), 19); ++nmtests;
        if (mas11_4 >= thr.lds11_4lo[1]) set_bit_raw(testbits.data(), 19);
        float lo2[2] = {thr.lds11_4lo[0], thr.lds11_4hi[0]};
        float hi2[2] = {thr.lds11_4lo[2], thr.lds11_4hi[2]};
        float md2[2] = {thr.lds11_4lo[1], thr.lds11_4hi[1]};
        cmin2 = std::min(cmin2, conf_test_2val(mas11_4, lo2, hi2, 4.0f, md2, 2));
        ++ng2;
    }
    if (ng2 > 0) nmtests += (ng2 > 1 ? 2 : 1);  // approximate

    // Group 3
    float cmin3 = 1.0f; int ng3 = 0;
    if (visusd && b2 > 0) {
        set_bit_raw(qa_bits.data(), 20); ++nmtests;
        if (b2 <= thr.ldsref2[1]) set_bit_raw(testbits.data(), 20);
        cmin3 = std::min(cmin3, conf_test(b2, thr.ldsref2[0], thr.ldsref2[2], thr.ldsref2[3], thr.ldsref2[1], 1));
        ++ng3;
    }

    // Group 4
    float cmin4 = 1.0f;
    if (visusd && b5 > 0) {
        set_bit_raw(qa_bits.data(), 16); ++nmtests;
        if (b5 <= thr.ldsref3[1]) set_bit_raw(testbits.data(), 16);
        cmin4 = std::min(cmin4, conf_test(b5, thr.ldsref3[0], thr.ldsref3[2], thr.ldsref3[3], thr.ldsref3[1], 1));
        if (b5 >= thr.ldstci[1] && b5 < thr.ldstci[0]) {
            clear_bit_raw(testbits.data(), 9); cirrus_vis = true;
        }
    }

    int groups = 0; float prod = 1.0f;
    if (ng1>0){++groups;prod*=cmin1;} if (ng2>0){++groups;prod*=cmin2;}
    if (ng3>0){++groups;prod*=cmin3;} if (cmin4<1.0f){++groups;prod*=cmin4;}
    if (groups>0) confdnc = std::pow(prod, 1.0f/groups);
}

// =========================================================================
// 2. LandDay_coast — daytime coastal land
// =========================================================================
inline void land_day_coast(const float pxldat[INBAND], float vza,
                           bool visusd, bool& cirrus_vis,
                           const LandDayThr& thr, const PfmftNfmftThr& pf,
                           int is_cold_sfc,
                           CloudMaskBits& testbits, QABits& qa_bits,
                           int& nmtests, float& confdnc, const float btclr[7]) {
    // Same structure as LandDay but uses _t2 suffixed thresholds
    // For FY-3D, the coastal thresholds (_t2) are identical to LandDay
    land_day(pxldat, vza, visusd, true, cirrus_vis, false, thr, pf, is_cold_sfc,
             testbits, qa_bits, nmtests, confdnc, btclr);
}

// =========================================================================
// 3. PolarDay_desert / PolarDay_desert_c / PolarDay_coast
//    These reuse polar_day_land with polar-specific desert thresholds
// =========================================================================

// PolarDay_desert: same as LandDay_desert but with polar thresholds (pds*)
inline void polar_day_desert(const float pxldat[INBAND], float vza,
                             bool visusd, bool& cirrus_vis,
                             const PolarDesertDayThr& thr, const PfmftNfmftThr& pf,
                             int is_cold_sfc,
                             CloudMaskBits& testbits, QABits& qa_bits,
                             int& nmtests, float& confdnc, const float btclr[7]) {
    float m31 = pxldat[23], m32 = pxldat[24], m20 = pxldat[19];
    float b2 = pxldat[1], b5 = pxldat[17];

    float cmin1 = 1.0f; int ng1 = 0;
    const float nf_desert[4] = {pf.nfmft_desert[0], pf.nfmft_desert[1], pf.nfmft_desert[2], pf.nfmft_desert[3]};
    desert_group1(pxldat, pf, is_cold_sfc, btclr, testbits, qa_bits, nmtests, ng1, cmin1, nullptr, nf_desert);

    float cmin2 = 1.0f; int ng2 = 0;
    if (m31 > 0 && m32 > 0 && vza > 0) {
        float cv = std::cos(vza * DTOR);
        if (std::abs(cv) > 1e-6f) {
            float dt = tview(1, 1.0f/cv, m31);
            float df = (dt >= 0.1f) ? dt : thr.pds11_12hi[0];
            if (m31-m32 <= df) { set_bit_raw(testbits.data(), 18); set_bit_raw(qa_bits.data(), 18); }
            float lo[2] = {df + 0.3f*df, thr.pds11_12hi[1]};
            float hi[2] = {df - 0.3f*df, thr.pds11_12hi[2]};
            float md[2] = {df, df};
            cmin2 = std::min(cmin2, conf_test_2val(m31-m32, lo, hi, 1.0f, md, 2));
            ++ng2; ++nmtests;
        }
    }
    if (visusd && m31 > 0 && m20 > 0) {
        float mas11_4 = m31 - m20;
        set_bit_raw(qa_bits.data(), 19); ++nmtests;
        if (mas11_4 >= thr.pds11_4lo[1]) set_bit_raw(testbits.data(), 19);
        float lo2[2] = {thr.pds11_4lo[0], thr.pds11_4hi[0]};
        float hi2[2] = {thr.pds11_4lo[2], thr.pds11_4hi[2]};
        float md2[2] = {thr.pds11_4lo[1], thr.pds11_4hi[1]};
        cmin2 = std::min(cmin2, conf_test_2val(mas11_4, lo2, hi2, 4.0f, md2, 2));
        ++ng2;
    }
    if (ng2 > 0) nmtests += 2;

    float cmin3 = 1.0f; int ng3 = 0;
    if (visusd && b2 > 0) {
        set_bit_raw(qa_bits.data(), 20); ++nmtests;
        if (b2 <= thr.pdsref2[1]) set_bit_raw(testbits.data(), 20);
        cmin3 = std::min(cmin3, conf_test(b2, thr.pdsref2[0], thr.pdsref2[2], thr.pdsref2[3], thr.pdsref2[1], 1));
        ++ng3;
    }

    float cmin4 = 1.0f;
    if (visusd && b5 > 0) {
        set_bit_raw(qa_bits.data(), 16); ++nmtests;
        if (b5 <= thr.pdsref3[1]) set_bit_raw(testbits.data(), 16);
        cmin4 = std::min(cmin4, conf_test(b5, thr.pdsref3[0], thr.pdsref3[2], thr.pdsref3[3], thr.pdsref3[1], 1));
        if (b5 >= thr.pdstci[1] && b5 < thr.pdstci[0]) { clear_bit_raw(testbits.data(), 9); cirrus_vis = true; }
    }

    int groups = 0; float prod = 1.0f;
    if (ng1>0){++groups;prod*=cmin1;} if (ng2>0){++groups;prod*=cmin2;}
    if (ng3>0){++groups;prod*=cmin3;} if (cmin4<1.0f){++groups;prod*=cmin4;}
    if (groups>0) confdnc = std::pow(prod, 1.0f/groups);
}

// Desert_c and Coast variants reuse polar_day_land with appropriate thresholds
inline void polar_day_desert_c(const float pxldat[INBAND], float vza,
                               bool visusd, bool& cirrus_vis,
                               const PolarDesertDayThr& thr, const PfmftNfmftThr& pf,
                               int is_cold_sfc,
                               CloudMaskBits& testbits, QABits& qa_bits,
                               int& nmtests, float& confdnc, const float btclr[7]) {
    polar_day_desert(pxldat, vza, visusd, cirrus_vis, thr, pf, is_cold_sfc,
                     testbits, qa_bits, nmtests, confdnc, btclr);
}

inline void polar_day_coast(const float pxldat[INBAND], float vza,
                            bool visusd, bool vrused, bool& cirrus_vis,
                            bool hi_elev,
                            const PolarDayLandThr& thr, const PfmftNfmftThr& pf,
                            int is_cold_sfc,
                            CloudMaskBits& testbits, QABits& qa_bits,
                            int& nmtests, float& confdnc, const float btclr[7]) {
    polar_day_land(pxldat, vza, visusd, vrused, cirrus_vis, hi_elev, thr, pf,
                   is_cold_sfc, testbits, qa_bits, nmtests, confdnc, btclr);
}

// =========================================================================
// 4. Antarctic_day — daytime Antarctic snow
// =========================================================================
inline void antarctic_day(const float pxldat[INBAND], bool visusd,
                          const AntarcticDayThr& thr, const PfmftNfmftThr& pf,
                          int is_cold_sfc,
                          CloudMaskBits& testbits, QABits& qa_bits,
                          int& nmtests, float& confdnc, const float btclr[7]) {
    float m31 = pxldat[23], m32 = pxldat[24], m20 = pxldat[19];
    float cmin1 = 1.0f, cmin2 = 1.0f;
    int ng1 = 0, ng2 = 0;

    // pfmft
    if (m31 > 0 && m32 > 0 && m31 < pf.pfmft_11maxthre[0] &&
        (btclr[4]-btclr[5]) > pf.pfmft_btd_min[0]) {
        float tv11_12;
        if (m31 > 270.0f && btclr[4] > 270.0f)
            tv11_12 = (m31-m32) - (btclr[4]-btclr[5])*(m31-260.0f)/(btclr[4]-260.0f);
        else
            tv11_12 = m31 - m32;
        set_bit_raw(testbits.data(), 14); set_bit_raw(qa_bits.data(), 14);
        const float* thr_pf = (is_cold_sfc == 1) ? pf.pfmft_cold : pf.pfmft_land;
        conf_test(tv11_12, thr_pf[0], thr_pf[2], thr_pf[3], thr_pf[1], 1);
    }
    // nfmft
    if (m31 > 0 && m32 > 0 && (m31-m32) <= pf.nfmft_maxthre[0]) {
        float tv11_12 = (m31-m32) - (btclr[4]-btclr[5]);
        set_bit_raw(testbits.data(), 15); set_bit_raw(qa_bits.data(), 15);
        conf_test(tv11_12, pf.nfmft_land[0], pf.nfmft_land[2], pf.nfmft_land[3], pf.nfmft_land[1], 1);
    }

    // Group 2: 4-11 um BT with dynamic thresholds
    if (visusd && m31 > 0 && m20 > 0 && m31 > 230.0f) {
        ++nmtests;
        set_bit_raw(qa_bits.data(), 19);
        float mas4_11 = m20 - m31;
        float lo, hi, mid, pw;
        get_pn_thresholds(m31, thr.bt_11_bnds4, thr.ant4_11l, thr.ant4_11m1,
                          thr.ant4_11m1, thr.ant4_11m1, thr.ant4_11h,
                          lo, hi, mid, pw);
        if (mas4_11 <= mid) set_bit_raw(testbits.data(), 19);
        cmin2 = std::min(cmin2, conf_test(mas4_11, lo, hi, pw, mid, 1));
        ++ng2;
    }
    if (ng2 > 0) nmtests += ng2;

    int groups = 0; float prod = 1.0f;
    if (ng1>0){++groups;prod*=cmin1;} if (ng2>0){++groups;prod*=cmin2;}
    if (groups>0) confdnc = std::pow(prod, 1.0f/groups);
}

}  // namespace fylat
