# FYLAT — FY-3D MERSI-II Cloud Mask Retrieval System

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Version](https://img.shields.io/badge/version-0.1.0-green.svg)](https://github.com/Rhett12345/cloudmask_fast)

**FYLAT (FengYun Land-Atmosphere Toolkit)** is a cloud detection retrieval system for the FY-3D MERSI-II medium-resolution spectral imager. It ingests L1 HDF5 satellite observations, combines numerical weather prediction (NWP) and ancillary data, and applies a multi-spectral decision-tree algorithm (derived from the MODIS MOD35 framework) to produce L2 cloud mask (CLM) products.

**Developed by the Min Min team at Sun Yat-sen University (中山大学).**

## Processing Pipeline

```
FY-3D MERSI-II L1 HDF5 ──┬── NWP (GRIB2) interpolation
                          ├── Ancillary data (emissivity, albedo, snow/ice, ecosystem)
                          ├── Radiative transfer model (PFAAST)
                          └── Decision-tree cloud mask ──→ L2 CLM/CLA/SST (HDF5)
```

1. **Read namelist config** (`.nml`) — sensor ID, NWP source, algorithm toggles, I/O paths
2. **Load L1 data** — 25 MERSI-II channels (19 reflective + 6 IR bands), geometry
3. **Load NWP** (GRIB2 → binary via wgrib2) — interpolate spatially to 101 pressure levels, temporally between two forecast steps
4. **Load ancillary data** — snow/ice mask (NISE), ecosystem type (IGBP), OISST, surface emissivity/albedo → run PFAAST RTM for clear-sky brightness temperatures
5. **Run cloud mask** — per-pixel decision tree (polar/land/water modules, ~20 scene-type subroutines), 3×3 spatial uniformity tests, thin cirrus/shadow/sun-glint checks, QA bit packing
6. **Output L2 HDF5** — cloud mask (1 km, 6-byte bit array + QA), optional cloud amount (5 km) and SST

### Supported Sensors

| Sensor ID | Description |
|-----------|-------------|
| 21 | FY-3D MERSI-II |
| 22 | FY-3E MERSI-II |
| 1–3 | MODIS / VIIRS → MERSI-II format conversion |

### Supported NWP Sources (10)

GFS (various resolutions), NCEP Reanalysis, T639, GRAPES GFS, GDAS — GRIB1/GRIB2, 26–41 layers.

## Quick Start

### Requirements

- Linux x86_64
- Intel Fortran (ifort) + Intel C (icc)
- HDF5 / HDF4 libraries
- wgrib / wgrib2
- Python 3 (for batch processing)

### Build

```bash
conda activate cloudmask
./build.sh              # runs: cd src && make -f fylat_makefile_cldmask
```

### Run (Single)

```bash
./fylat_FY3_MERSI_II_PGS config.nml
```

### Run (Batch)

```bash
python3 paral_bat_driver_mersi_ii_fylat.py
```

## Project Structure

```
├── src/                    # Fortran/C source (~37,000 lines)
│   ├── fylat_FY3_MERSI_II_PGS_Driver.f90  # Main driver
│   ├── cloudmask/          # Cloud mask algorithm core
│   ├── cloudamount/        # Cloud amount algorithm
│   └── *.f90, *.f, *.c     # Support modules (IO, RTM, numerics)
├── coeff/                  # Coefficient files (thresholds, RTM tables)
├── include/                # Fortran include files
├── wgrib/                  # NWP GRIB2 → binary conversion scripts
├── drivefile/              # Namelist config files
├── convert_oisst/          # OISST NC → HDF5 conversion
└── paral_bat_driver_mersi_ii_fylat.py  # Batch processing driver
```

## Cloud Mask Algorithm

Based on the MODIS MOD35 framework, adapted for the FY-3D MERSI-II 25-channel configuration.

### Decision Tree Logic

| Surface Type | Region | Illumination | Modules |
|-------------|--------|--------------|---------|
| Land/Desert/Coast/Snow | Non-polar | Day | LandDay, LandDay_desert, LandDay_coast |
| Land/Desert/Coast/Snow | Non-polar | Night | LandNite |
| Ocean | Non-polar | Day/Night | ocean_day, ocean_nite |
| Land/Ocean/Snow | Polar (>60°) | Day/Night | PolarDay_*, PolarNite_* |
| — | Antarctica (<-60°) | Day | Antarctic_day |

### Key Spectral Tests

- Visible reflectance thresholds (0.66, 0.86, 1.38 μm)
- NDSI snow detection (0.55 μm vs 1.64/2.13 μm)
- IR brightness temperature differences (BTD 8.5–11, 3.8–11, 11–12)
- Visible ratio test (VRAT = R0.86/R0.65)
- 3×3 spatial uniformity test
- Thin cirrus detection (1.38 μm reflectance + IR BTD)
- Sun glint detection (glint angle ≤ 36°)

### Cloud Confidence Levels

| Level | Meaning |
|-------|---------|
| 0 | Cloudy |
| 1 | Probably Cloudy |
| 2 | Probably Clear |
| 3 | Clear |

## Team

- **Principal Investigator**: Min Min (闵敏), Sun Yat-sen University
- **Developer**: Yu Qiang (余强), Sun Yat-sen University — yuqiang6@mail2.sysu.edu.cn, 2942204121@qq.com

## License

MIT License — see [LICENSE](LICENSE) for details.

## Citation

If you use this software in your research, please cite:

```
Min, M., & Yu, Q. (2025). FYLAT: FY-3D MERSI-II Cloud Mask Retrieval System (Version 3.1).
Sun Yat-sen University. https://github.com/Rhett12345/cloudmask_fast
```
