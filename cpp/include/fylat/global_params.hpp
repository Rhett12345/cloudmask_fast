#pragma once
// Cloud mask algorithm parameters — translated from global.inc

namespace fylat {

// Processing context box dimensions (3x3 window)
constexpr int NLCNTX = 3;
constexpr int NECNTX = 3;

// Channel / band counts
constexpr int INBAND       = 25;   // total bands used in cloud mask
constexpr int VIS_BAND     = 4;    // 250 m visible band reflectances
constexpr int BANDS_USED   = 19;   // actual number of bands used
constexpr int SG_BANDS_USED = 5;   // bands used in sun-glint algorithm
constexpr int IR_BAND      = 6;    // IR channels
constexpr int VAR_BAND     = 2;    // bands for spatial variability tests
constexpr int NUM250_PER_1KM = 16; // 250 m FOVs within 1 km footprint

// Cloud mask byte / bit dimensions
constexpr int CM_BYTE_DIM = 6;   // cloud mask bit array (48 bits)
constexpr int CM_QA_DIM   = 10;  // QA bit array (80 bits)

// Sentinel values
constexpr float BAD_DATA   = -999.0f;
constexpr float MISG       = -99.0f;
constexpr float META_MISG  = -99999.0f;

}  // namespace fylat
