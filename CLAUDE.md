# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

FYLAT — FY-3D MERSI-II Cloud Mask Retrieval System V3.1. Reads FY-3D/E MERSI-II satellite L1 HDF5 data, runs a MODIS MOD35-derived decision-tree cloud detection algorithm, and outputs L2 cloud mask (CLM) and related products in HDF5. Developed by the Min Min team at Sun Yat-sen University (中山大学). Team: Min Min (闵敏), Yu Qiang (余强, yuqiang6@mail2.sysu.edu.cn).

## Build & Run

```bash
# Conda environment required
conda activate cloudmask

# Build (Intel ifort + icc + HDF5 required)
./build.sh            # runs: cd src && make -f fylat_makefile_cldmask

# Clean
cd src && make clean -f fylat_makefile_cldmask && cd ..

# One-command end-to-end run (auto-discover, build, both calibrations, verify)
python run_fylat.py --date 20220803
python run_fylat.py --date 20220803 --cores 4

# Single run (manual)
./fylat_FY3_MERSI_II_PGS config.nml

# Single run (YAML-based)
python scripts/run_single.py --date 20220803 --time 0740 --calibration business
```

Compiler toolchain: `ifort` (Fortran), `icc` (C), `h5fc` (HDF5 Fortran wrapper). HDF4/HDF5 libs at `/opt/intel_lib/lib/` — hardcoded in the makefile.

## Architecture

**Processing flow** (6 steps orchestrated by `src/fylat_FY3_MERSI_II_PGS_Driver.f90`):

1. Read Fortran namelist config (`.nml`) — sensor ID, NWP source, I/O paths, algorithm toggles
2. Load FY-3D MERSI-II GEO + L1b HDF5 (25 channels: 19 reflective + 6 IR) — `io_module.f90`
3. Load NWP data (GRIB2 → Fortran-compatible binary via Python ecCodes by default; wgrib2 fallback optional), interpolate to observation time and 101 pressure levels — `python/fylat/nwp_reader.py`, `read_nwp_data_module.f90`, `nwp_utils_module.f90`
4. Load ancillary data (snow/ice, ecosystem, emissivity, albedo, OISST), run RTM (PFAAST) for clear-sky brightness temperatures — `get_ancil_data_module.f90`, `rtm_tran_module.f90`, `rtm_utils_module.f90`
5. Run cloud mask decision tree per-pixel — `fylat_fy3mersi_cloud_mask.f90`
6. Write L2 HDF5 output, optionally cloud amount (5 km) and SST

**Cloud mask algorithm** (`src/cloudmask/`):

- Top-level: `fylat_fy3mersi_cloud_mask.f90` → dispatches to `polar_module.f90`, `land_module.f90`, `water_module.f90`
- Each module routes to scene-specific subroutines: e.g. `LandDay.f90`, `LandNite.f90`, `PolarDay_land.f90`, `ocean_day.f90`, etc. (~20 scene modules)
- Threshold parameters defined in paired `*_thr.inc` include files
- Bit-level QA operations in `check_bits.f`, `set_bit.f`, `set_qa_bit.f`, `clear_bit.f`, `set_confdnc.f`
- Ancillary checks: `shadows.f90` (cloud shadows), `thin_ci_chk_ir.f90` (thin cirrus), `chk_sunglint.f90` (sun glint), `noncld_obs_chk.f90`

**Key support modules:**

| Module | Purpose |
|---|---|
| `platform_module.f90` | Sensor channel parameters (wavelengths, wavenumbers) — hardcoded per sensor_id |
| `data_arrays_module.f90` | Shared global derived types: `sat`, `geo`, `nwp_*` (7 NWP variants), `ancil` |
| `planck_module.f90` | Planck radiance ↔ brightness temperature conversion |
| `frontend_module.f90` | Solar/sensor geometry, scattering angles, sun-earth distance |
| `numerical.f90` | Interpolation routines |
| `constant.f90` | Physical constants |
| `names_module.f90` | Magic number parameters (channel counts, array dimensions) |

**Algorithm modules not currently built** (source exists but rules commented out in makefile):
`cloudamount/`, `cloudphase/`, `cloudheight/`, `cloudod_day/`, `sea_surface_temperature/`

## Configuration

Namelist format (`.nml`). Key parameters:
- `fylat_sensor_id` — 21=FY-3D, 22=FY-3E (primary); 1/2/3 for MODIS/VIIRS conversion
- `fylat_nwp_opt` — 1–10 selecting NWP source (GFS, NCEP, T639, GRAPES, GDAS)
- `cloudmask_id` / `cloudamount_id` / `surface_sst_id` — algorithm toggles (0/1)
- See `temp_fy3d_config_20200330_0000.nml` for a complete example
- `drivefile/` has many scenario-specific configs

## Important Codebase Conventions

- **Mixed Fortran 77/90**: `.f` files are fixed-format F77 (many lack `implicit none`); `.f90` files are free-format. The makefile uses `ifort` for both but `h5fc` for HDF5-dependent modules.
- **No version control**: No `.git` directory. Consider this before destructive edits.
- **Hardcoded paths**: Makefile uses `FYLAT_ROOT_PATH` and `HDF5_ROOT` env vars (with defaults). Namelist files contain absolute data paths. Changing machines requires updating these.
- **Debug artifacts**: Variables prefixed `lyj` (e.g., `out_*` arrays in cloudmask) and comments marked `jincheng test` are developer debug markers, not production code.
- **Comment-out convention**: Disabled code is commented out (not deleted or `#ifdef`-gated), across both Fortran and the makefile.

## 特别说明
- 每次回答前都需要说一句“打报告”，然后再进行回答
- 测试的日期为20220803，该日期所对应的MYD35数据集为标准验证集
- 测试的输出目录为/data/Data_yuq/fy3_cloud,在该目录下创建yyyymmdd目录进行存放
- /data/Data_yuq/mersi和/data/Data_yuq/aqua_modis/MYD35_L2/,/data/nwp/为现阶段所需数据的存放目录
- 每次都需要测试两个定标系数，针对于太阳反射波段的七个通道，分别分为数据文件内置的业务定标以及我存放在../fy3d_recali按日期分类的再定标系数
- 虚拟环境为cloudmask
