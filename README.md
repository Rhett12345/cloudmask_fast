# FYLAT — FY-3D MERSI-II Cloud Mask Retrieval System

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Version](https://img.shields.io/badge/version-2.3.2-green.svg)](https://github.com/Rhett12345/cloudmask_fast)

**FYLAT (FengYun Land-Atmosphere Toolkit)** is a cloud detection retrieval system for the FY-3D/E MERSI-II medium-resolution spectral imager. It ingests L1 HDF5 satellite observations, combines numerical weather prediction (NWP) and ancillary data, and applies a multi-spectral decision-tree algorithm (derived from the MODIS MOD35 framework) to produce L2 cloud mask (CLM) products.

Developed by the **Min Min team at Sun Yat-sen University (中山大学)**.

## Quick Start

### Requirements

- Linux x86_64
- Intel Fortran (`ifort`) + Intel C (`icc`)
- HDF5 / HDF4 libraries
- ecCodes for default GRIB2 decoding; wgrib2 is retained as an explicit fallback
- Python 3, conda environment `cloudmask`

### One-Command Run

```bash
conda activate cloudmask
python run_fylat.py --date 20220803
```

This single command:
1. Auto-discovers all time slots (GEO + L1B pairs)
2. Auto-discovers NWP GRIB2 forecast files
3. Builds the Fortran executable (skipped if already built)
4. Runs cloud mask retrieval for **BUSINESS** (onboard) AND **RECALI** (external) calibrations
5. Verifies all output files

```bash
python run_fylat.py --date 20220803 --cores 4    # use 4 CPU cores
python run_fylat.py --date 20220803 --dry-run    # preview config only
python run_fylat.py --date 20220803 --skip-build # skip compilation
python run_fylat.py --date 20220803 --verify-only  # verify existing outputs only
```

### Calibration Modes

| Mode | Description |
|------|-------------|
| `business` | HDF5 built-in calibration coefficients |
| `recali` | External recalibration from `../fy3d_recali/YYYYMM/RAD_*.h5` (auto-loaded per date) |

### Manual Build & Run

```bash
conda activate cloudmask
./build.sh                                    # compile Fortran
./fylat_FY3_MERSI_II_PGS config.nml          # single run with namelist
python scripts/run_single.py --date 20220803 --time 0740 --calibration business
```

## Processing Pipeline

```
FY-3D MERSI-II L1 HDF5 ──┬── NWP (GRIB2) interpolation
                          ├── Ancillary data (emissivity, albedo, snow/ice, ecosystem)
                          ├── Radiative transfer model (PFAAST)
                          └── Decision-tree cloud mask ──→ L2 CLM (HDF5)
```

1. Read namelist config (`.nml`) — sensor ID, NWP source, algorithm toggles, I/O paths
2. Load L1 data — 25 MERSI-II channels (19 reflective + 6 IR bands), geometry
3. Load NWP (GRIB2 → Fortran-compatible binary via Python ecCodes by default, wgrib2 fallback optional) — interpolate to 101 pressure levels, temporally between forecast steps
4. Load ancillary data — snow/ice mask (NISE), ecosystem (IGBP), OISST, emissivity/albedo → PFAAST RTM for clear-sky BTs
5. Run cloud mask — per-pixel decision tree, 3×3 spatial tests, thin cirrus/shadow/sun-glint checks, QA bit packing
6. Output L2 HDF5 — cloud mask (1 km, 6-byte bit array + QA)

## Project Structure

```
├── src/                         # Fortran/C source
│   ├── fylat_FY3_MERSI_II_PGS_Driver.f90   # Main driver
│   ├── fylat_makefile_cldmask              # Makefile
│   ├── cloudmask/               # Cloud mask algorithm (~20 scene modules)
│   ├── cloudamount/             # Cloud amount (not built by default)
│   ├── cloudphase/              # Cloud phase (not built by default)
│   ├── cloudheight/             # Cloud top height (not built by default)
│   ├── cloudod_day/             # Cloud optical depth (not built by default)
│   ├── sea_surface_temperature/ # SST (not built by default)
│   └── *.f90, *.f, *.c         # Support: IO, RTM, numerics, platform
├── coeff/                       # Coefficient files (thresholds, RTM LUTs, ancillary)
├── wgrib/                       # Optional legacy NWP GRIB2 fallback scripts
├── python/fylat/                # Python config & calibration management
│   ├── config.py                # YAML → .nml namelist generation
│   └── calibration.py           # Recalibration coefficient discovery & loading
├── config/                      # YAML configuration
│   ├── default.yaml             # Default settings
│   └── scenes/                  # Per-scene overrides
├── scripts/
│   └── run_single.py            # Single-scene runner (YAML-based)
├── visualize/                   # Validation & visualization tools
├── build.sh                     # Build script
├── clean.sh                     # Clean build artifacts
└── run_fylat.py                 # One-command end-to-end driver
```

## Cloud Mask Algorithm

Derived from MODIS MOD35 framework, adapted for FY-3D MERSI-II 25-channel configuration.

### Decision Tree

| Surface Type | Region | Illumination | Modules |
|-------------|--------|--------------|---------|
| Land/Desert/Coast/Snow | Non-polar | Day | LandDay, LandDay_desert, LandDay_coast |
| Land/Desert/Coast/Snow | Non-polar | Night | LandNite |
| Ocean | Non-polar | Day/Night | ocean_day, ocean_nite |
| Land/Ocean/Snow | Polar (>60°) | Day/Night | PolarDay_*, PolarNite_* |
| — | Antarctica (<-60°) | Day | Antarctic_day |

### Key Spectral Tests

- Visible reflectance (0.66, 0.86, 1.38 μm)
- NDSI snow detection (0.55 vs 1.64/2.13 μm)
- IR BTD (8.5–11, 3.8–11, 11–12 μm)
- Visible ratio test (R0.86/R0.65)
- 3×3 spatial uniformity
- Thin cirrus detection (1.38 μm + IR BTD)
- Sun glint (angle ≤ 36°)

### Cloud Confidence

| Level | Meaning |
|-------|---------|
| 0 | Cloudy |
| 1 | Probably Cloudy |
| 2 | Probably Clear |
| 3 | Clear |

## Team

- **PI**: Min Min (闵敏), Sun Yat-sen University
- **Developer**: Yu Qiang (余强), Sun Yat-sen University — yuqiang6@mail2.sysu.edu.cn

## License

MIT License — see [LICENSE](LICENSE).

## Citation

```
Min, M., & Yu, Q. (2025). FYLAT: FY-3D MERSI-II Cloud Mask Retrieval System (Version 3.1).
Sun Yat-sen University. https://github.com/Rhett12345/cloudmask_fast
```
