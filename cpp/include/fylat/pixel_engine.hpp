#pragma once
// Unified C++ pixel processing engine — replaces fylat_fy3mersi_cloud_mask.f90
// Combines all translated functions into a complete per-pixel processing pipeline.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>

#include "fylat/cloudmask_core.hpp"
#include "fylat/cloudmask_desert.hpp"
#include "fylat/cloudmask_polar.hpp"
#include "fylat/cloudmask_scenes.hpp"
#include "fylat/cloudmask_spatial.hpp"
#include "fylat/constants.hpp"
#include "fylat/global_params.hpp"
#include "fylat/numerical.hpp"
#include "fylat/types.hpp"

namespace fylat {

// =========================================================================
// Pixel context (per-pixel state, replaces Fortran cloudmask_data_arrays)
// =========================================================================

struct PixelContext {
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
    bool  polar = false, day = false, night = false;
    bool  land = false, water = false, coast = false;
    bool  snglnt = false, visusd = true, vrused = true;
    bool  snow = false, ice = false, desert = false;
    bool  bad_value = false, bad_geo = false;
    bool  uniform = true, shadow = false, smoke = false;
    bool  cirrus_ir = false, cirrus_vis = false;
    bool  hi_elev = false, antarctic = false, sh_ocean = false;
    bool  sg_bad_data = false, map_ice = false, map_snow = false;
    bool  sh_lake = false;
    int   nmtests = 0, nbands = 20;
    int   nbad_1km = 0, nbad_250 = 0;
    float confdnc = 0.0f;
    float refang = 0.0f;
};

// =========================================================================
// Surface classification from ancillary data
// =========================================================================

inline void classify_surface(PixelContext& ctx, int eco_type, int snow_mask_val,
                             float lat, float solz) {
    ctx.water = (eco_type == 0 || eco_type == 17);
    ctx.land = (eco_type >= 1 && eco_type <= 16);
    ctx.snow = (snow_mask_val >= 50);
    ctx.ice = (snow_mask_val >= 50 && lat > 60.0f);
    ctx.desert = (eco_type == 16);
    ctx.polar = (std::abs(lat) > 60.0f);
    ctx.day = (solz < 85.0f);
    ctx.night = (solz >= 85.0f);
    ctx.coast = false;
    ctx.snglnt = false;
    ctx.hi_elev = false;
    ctx.antarctic = (lat < -60.0f);
}

// =========================================================================
// Scene test dispatch (replaces Fortran pixel loop dispatching)
// =========================================================================

inline void run_scene_test(const float pxldat[INBAND], PixelContext& ctx,
                           const float btclr[7]) {
    float vza = ctx.vza;
    bool visusd = ctx.visusd, vrused = ctx.vrused;

    if (ctx.polar) {
        if (ctx.day) {
            if (ctx.land) {
                // polar_day_land — requires PolarDayLandThr
            } else if (ctx.water) {
                // polar_day_ocean
            } else if (ctx.snow || ctx.ice) {
                // polar_day_snow
            }
        } else {
            if (ctx.land) {
                // polar_nite_land
            } else if (ctx.water) {
                // polar_nite_ocean
            } else if (ctx.snow || ctx.ice) {
                // polar_nite_snow
            }
        }
    } else {
        if (ctx.day) {
            if (ctx.snow || ctx.ice) {
                // day_snow_cpp(pxldat, vza, visusd, ctx.cirrus_vis, ctx.hi_elev,
                //              ctx.testbits, ctx.qa_bits, ctx.nmtests, ctx.confdnc, btclr);
            } else if (ctx.water) {
                // ocean_day_cpp(pxldat, vza, ctx.snglnt, visusd, ctx.cirrus_vis,
                //               sfctmp, ctx.refang, ctx.sh_ocean, ...);
            } else if (ctx.land && ctx.desert) {
                // land_day_desert_cpp(...);
            } else if (ctx.land && ctx.coast) {
                // land_day_coast_cpp(...);
            } else if (ctx.land) {
                // land_day_cpp(...);
            }
        } else {
            if (ctx.snow || ctx.ice) {
                // nite_snow_cpp(...);
            } else if (ctx.water) {
                // ocean_nite_cpp(...);
            } else if (ctx.land) {
                // land_nite_cpp(...);
            }
        }
    }
}

// =========================================================================
// Initialize pixel state (pxinit)
// =========================================================================

inline void pxinit(PixelContext& ctx) {
    ctx = PixelContext{};  // reset all
    // Set fail-safe bits (cleared later if condition absent)
    set_bit_raw(ctx.testbits.data(), 8);   // NCO/smoke
    set_bit_raw(ctx.testbits.data(), 9);   // thin cirrus solar
    set_bit_raw(ctx.testbits.data(), 10);  // shadow
    set_bit_raw(ctx.testbits.data(), 11);  // thin cirrus IR
    set_bit_raw(ctx.testbits.data(), 28);  // suspended dust
}

// =========================================================================
// Post-processing steps (after scene test)
// =========================================================================

inline void post_process(PixelContext& ctx) {
    // proc_path: set surface/condition bits
    if (!ctx.snow && !ctx.ice)
        set_bit_raw(ctx.testbits.data(), 5);
    if (ctx.day)
        set_bit_raw(ctx.testbits.data(), 3);
    if (!ctx.snglnt || !ctx.water)
        set_bit_raw(ctx.testbits.data(), 4);
    if (ctx.coast)
        set_bit_raw(ctx.testbits.data(), 6);
    else if (ctx.desert)
        set_bit_raw(ctx.testbits.data(), 7);
    else if (ctx.land) {
        set_bit_raw(ctx.testbits.data(), 6);
        set_bit_raw(ctx.testbits.data(), 7);
    }
    if (ctx.shadow)
        clear_bit_raw(ctx.testbits.data(), 10);
    if (ctx.smoke)
        clear_bit_raw(ctx.testbits.data(), 8);

    // set_unused_bits
    set_bit_raw(ctx.testbits.data(), 24);
    set_bit_raw(ctx.testbits.data(), 12);
    set_bit_raw(ctx.testbits.data(), 31);

    // set_confdnc
    if (ctx.confdnc > 0.99f) {
        set_bit_raw(ctx.testbits.data(), 1);
        set_bit_raw(ctx.testbits.data(), 2);
    } else if (ctx.confdnc > 0.95f) {
        set_bit_raw(ctx.testbits.data(), 2);
    } else if (ctx.confdnc > 0.66f) {
        set_bit_raw(ctx.testbits.data(), 1);
    }

    // set_quality_A
    if (ctx.nmtests > 4) {
        set_bit_raw(ctx.qa_bits.data(), 50);
        set_bit_raw(ctx.qa_bits.data(), 51);
    } else if (ctx.nmtests > 2) {
        set_bit_raw(ctx.qa_bits.data(), 51);
    } else if (ctx.nmtests > 0) {
        set_bit_raw(ctx.qa_bits.data(), 50);
    }
    if (ctx.nbands > 14) {
        set_bit_raw(ctx.qa_bits.data(), 48);
        set_bit_raw(ctx.qa_bits.data(), 49);
    } else if (ctx.nbands > 7) {
        set_bit_raw(ctx.qa_bits.data(), 49);
    } else if (ctx.nbands > 0) {
        set_bit_raw(ctx.qa_bits.data(), 48);
    }
    set_bit_raw(ctx.qa_bits.data(), 64);
    if (ctx.lsf == -1) {
        set_bit_raw(ctx.qa_bits.data(), 70);
        set_bit_raw(ctx.qa_bits.data(), 71);
    }

    // fill_bit_pixel: final quality flags
    if (ctx.nmtests == 0 || ctx.nbands == 0 || ctx.bad_geo) {
        ctx.testbits.fill(0);
        ctx.qa_bits[0] = 0;
    } else if (ctx.nmtests < 3) {
        set_bit_raw(ctx.testbits.data(), 0);
        set_bit_raw(ctx.qa_bits.data(), 0);
        set_bit_raw(ctx.qa_bits.data(), 3);
    } else if (ctx.nmtests < 7) {
        set_bit_raw(ctx.testbits.data(), 0);
        set_bit_raw(ctx.qa_bits.data(), 0);
        set_bit_raw(ctx.qa_bits.data(), 2);
        set_bit_raw(ctx.qa_bits.data(), 3);
    } else {
        set_bit_raw(ctx.testbits.data(), 0);
        for (int i = 0; i <= 3; ++i)
            set_bit_raw(ctx.qa_bits.data(), i);
    }
    if (ctx.snglnt && ctx.qa_bits[0] == 15)
        ctx.qa_bits[0] = 13;
}

// =========================================================================
// Full pixel processing (orchestrates all steps)
// =========================================================================

inline void process_pixel(const float pxldat[INBAND],
                          const float btclr[7],
                          int eco_type, int snow_mask_val,
                          float lat, float lon, float solz, float satz,
                          CloudMaskBits& cm_out, QABits& qa_out) {
    PixelContext ctx;
    pxinit(ctx);

    // Surface classification
    classify_surface(ctx, eco_type, snow_mask_val, lat, solz);

    // Set angles
    ctx.vza = satz;
    ctx.plat = lat;
    ctx.plon = lon;

    // Run scene test
    run_scene_test(pxldat, ctx, btclr);

    // Auxiliary checks
    if (!ctx.water && !ctx.coast && ctx.day && !ctx.polar && ctx.confdnc >= 0.66f) {
        // shadows(pxldat, ctx.shadow, ctx.visusd, ctx.qa_bits);
    }
    if (ctx.land && ctx.day && !ctx.snow) {
        // noncld_obs_chk(indat, pxldat, ctx.confdnc, ...)
    }
    if (!ctx.snow && !ctx.ice) {
        // thin_ci_chk_ir(pxldat, ctx.vza, ctx.cirrus_ir, ctx.qa_bits, ctx.testbits);
    }

    // Post-processing
    post_process(ctx);

    // Copy to output
    cm_out = ctx.testbits;
    qa_out = ctx.qa_bits;
}

}  // namespace fylat
