#pragma once
// Strongly-typed enum replacements for symbol_struct in constant.f90
// Fortran integer parameters → C++ enum class

#include <cstdint>

namespace fylat {

// ---------------------------------------------------------------------------
// Cloud mask confidence flags
// ---------------------------------------------------------------------------
enum class CloudConfidence : int8_t {
    CONFIDENT_CLOUDY   = 0,
    PROBABLY_CLOUDY    = 1,
    PROBABLY_CLEAR     = 2,
    CONFIDENT_CLEAR    = 3
};

// ---------------------------------------------------------------------------
// Cloud phase / type
// ---------------------------------------------------------------------------
enum class CloudPhase : int8_t {
    CLEAR_SKY    = 0,
    WATER_CLOUD  = 1,
    ICE_CLOUD    = 2,
    MIXED_PHASE  = 3,
    UNKNOWN      = 4
};

// ---------------------------------------------------------------------------
// Surface types (IGBP)
// ---------------------------------------------------------------------------
enum class SurfaceType : int8_t {
    WATER            = 0,
    EVERGREEN_NEEDLE = 1,
    EVERGREEN_BROAD  = 2,
    DECIDUOUS_NEEDLE = 3,
    DECIDUOUS_BROAD  = 4,
    MIXED_FOREST     = 5,
    CLOSED_SHRUB     = 6,
    OPEN_SHRUB       = 7,
    WOODY_SAVANNA    = 8,
    SAVANNA          = 9,
    GRASSLAND        = 10,
    WETLAND          = 11,
    CROPLAND         = 12,
    URBAN            = 13,
    CROP_MOSAIC      = 14,
    SNOW_ICE         = 15,
    BARREN           = 16,
    WATER_BODY       = 17,
    TUNDRA           = 18  // unused but kept for IGBP completeness
};

// ---------------------------------------------------------------------------
// Binary flags
// ---------------------------------------------------------------------------
enum class YesNo : int8_t { NO = 0, YES = 1 };
enum class SpaceMask : int8_t { NO_SPACE = 0, SPACE = 1 };
enum class LandWater : int8_t { WATER = 0, LAND = 1, COAST = 2 };
enum class DesertMask : int8_t { NO_DESERT = 0, DESERT = 1 };
enum class SnowMask : int8_t { NO_SNOW = 0, SNOW = 1 };
enum class VolcanoMask : int8_t { NO_VOLCANO = 0, VOLCANO = 1 };
enum class CoastMask : int8_t { NO_COAST = 0, COAST = 1 };

// ---------------------------------------------------------------------------
// Scale method
// ---------------------------------------------------------------------------
enum class ScaleMethod : int8_t {
    NONE        = 0,
    LINEAR      = 1,
    LOGARITHMIC = 2
};

// ---------------------------------------------------------------------------
// System / endian
// ---------------------------------------------------------------------------
enum class Endian : int8_t { LITTLE = 0, BIG = 1 };

// ---------------------------------------------------------------------------
// Aerosol types
// ---------------------------------------------------------------------------
enum class AerosolType : int8_t {
    NONE         = 0,
    SMOKE        = 1,
    DUST         = 2,
    VOLCANIC_ASH = 3,
    URBAN        = 4
};

// ---------------------------------------------------------------------------
// Emissivity / albedo type
// ---------------------------------------------------------------------------
enum class EmissivityType : int8_t { NONE = 0, STATIC = 1, DYNAMIC = 2 };
enum class AlbedoType    : int8_t { NONE = 0, STATIC = 1, DYNAMIC = 2 };

// ---------------------------------------------------------------------------
// Nighttime quality flags
// ---------------------------------------------------------------------------
enum class NightQuality : int8_t { LOW = 0, MEDIUM = 1, HIGH = 2 };

// ---------------------------------------------------------------------------
// Generic surface type (simplified)
// ---------------------------------------------------------------------------
enum class GenericSurface : int8_t {
    LAND     = 0,
    WATER    = 1,
    SNOW     = 2,
    ICE      = 3,
    DESERT   = 4,
    COAST    = 5
};

}  // namespace fylat
