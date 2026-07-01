# FYLAT 云检测系统使用说明

## 环境准备

```bash
conda activate cloudmask
```

依赖工具链：

| 工具 | 用途 | 路径 |
|------|------|------|
| `ifort` | Fortran 编译器 | PATH 中 |
| `icc` | C 编译器 | PATH 中 |
| `h5fc` | HDF5 Fortran 包装器 | PATH 中 |
| `ecCodes` | 默认 GRIB2 解码 | conda `cloudmask` 环境 |
| `wgrib2` | GRIB2 兼容 fallback | `/opt/software/grib2/wgrib2/wgrib2` |
| HDF5/HDF4 库 | 数据 I/O | `/opt/intel_lib/` |

## 一键运行（推荐）

```bash
python run_fylat.py --date 20220803
```

一条命令自动完成：

1. **发现时次** — 扫描 L1B 数据目录，找到所有 GEO + L1B 配对的观测时次
2. **发现 NWP** — 自动选择对应的 GFS 预报场 GRIB2 文件
3. **编译检查** — 若可执行文件不存在则自动编译
4. **双定标反演** — 同时运行 BUSINESS（业务）和 RECALI（再定标）模式
5. **输出验证** — 检查所有输出文件完整性

### 常用参数

```bash
# 指定并行核心数（默认全部可用）
python run_fylat.py --date 20220803 --cores 4

# 跳过编译（已编译过可节省时间）
python run_fylat.py --date 20220803 --skip-build

# 仅预览配置，不实际运行
python run_fylat.py --date 20220803 --dry-run

# 仅验证已有输出
python run_fylat.py --date 20220803 --verify-only

# 只跑一种定标
python run_fylat.py --date 20220803 --calibrations business
python run_fylat.py --date 20220803 --calibrations recali
```

### 完整参数列表

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--date` | 必填 | 观测日期，格式 YYYYMMDD |
| `--cores` | 全部可用 | 并行核心数 |
| `--calibrations` | `business,recali` | 定标模式，逗号分隔 |
| `--l1b-path` | `/data/Data_yuq/mersi/` | L1B 数据根目录 |
| `--nwp-path` | `/data/nwp/` | NWP 数据根目录 |
| `--oisst-path` | `/data/Data_minmin/oisst/` | OISST 海温数据目录 |
| `--output-path` | `/data/Data_yuq/fy3_cloud/` | 输出根目录 |
| `--skip-build` | false | 跳过编译 |
| `--verify-only` | false | 仅验证输出，不运行反演 |
| `--dry-run` | false | 仅生成配置，不运行 |

## 手动运行（调试用）

### 编译

```bash
./build.sh                          # 编译
cd src && make clean -f fylat_makefile_cldmask && cd ..  # 清理
```

编译产物：`fylat_FY3_MERSI_II_PGS`（项目根目录）

### YAML 方式（推荐的手动方式）

```bash
# 列出可用的定标模式
python scripts/run_single.py --list-calibrations

# 运行单个场景
python scripts/run_single.py --date 20220803 --time 0740 --calibration business

# 预览配置不运行
python scripts/run_single.py --date 20220803 --time 0740 --dry-run
```

配置体系：
- `config/default.yaml` — 默认参数
- `config/scenes/YYYYMMDD_HHMM.yaml` — 场景专属覆盖

### C++ IO 后端（M1/P1 迁移）

Python 侧的 `fylat.mersi_io` 已接入 C++ HDF5/pybind11 后端。构建后会优先用 C++ 读取 GEO/L1B 业务级 payload，属性读取、缩放和缺省兜底仍由 h5py/Python 负责。

```bash
mkdir -p build_migration
cd build_migration
pybind11_DIR=$(python -c "import pybind11; print(pybind11.get_cmake_dir())")
cmake .. -Dpybind11_DIR="$pybind11_DIR"
cmake --build . --target fylat_py -j2
ctest --output-on-failure
```

远端 `liusy2020@172.16.102.82` 推荐直接运行：

```bash
conda activate cloudmask
bash scripts/check_migration_build.sh /tmp/fylat_build_check
```

常用开关：

```bash
FYLAT_IO_BACKEND=auto  python scripts/run_single.py ...  # 默认：有 C++ 后端则使用，否则回退 h5py
FYLAT_IO_BACKEND=cpp   python -c "from fylat.mersi_io import io_backend_name; print(io_backend_name())"
FYLAT_IO_BACKEND=h5py  python scripts/run_single.py ...  # 强制回退旧 Python h5py 读取
```

### NWP GRIB2 后端（M2 迁移）

批处理入口会在运行 Fortran 前预生成 Fortran 兼容的 NWP 二进制文件。默认后端是 Python ecCodes，不再依赖旧的 `wgrib/` shell 脚本链；`wgrib2` 仅作为显式兼容 fallback。

```bash
FYLAT_NWP_BACKEND=eccodes python run_fylat.py --date 20220803  # 默认推荐
FYLAT_NWP_BACKEND=wgrib2  python run_fylat.py --date 20220803  # 显式兼容旧路径
```

如果 Fortran 当前读取的单时次 `gfs0p25_41L_{YYYYMMDD}_{HH}_00` 二进制文件已经存在，预处理会直接复用，不重复解析 GRIB2。

### 云检测核心后端（M3 预留）

当前生产默认仍是 Fortran 决策树。M3 的 C++ `ocean_day` 迁移会通过显式环境变量接入，默认行为不变。

```bash
FYLAT_CLOUDMASK_BACKEND=fortran        python run_fylat.py --date 20220803  # 默认
FYLAT_CLOUDMASK_BACKEND=cpp_ocean_day  python run_fylat.py --date 20220803  # M3 实验开关，占位
FYLAT_CLOUDMASK_BACKEND=auto           python run_fylat.py --date 20220803  # 后续自动路由，占位
```

### 传统 namelist 方式

```bash
./fylat_FY3_MERSI_II_PGS config.nml
```

`.nml` 文件可手动编写或由 `run_fylat.py --dry-run` 自动生成。

## 定标模式

| 模式 | 配置 | 定标系数来源 |
|------|------|-------------|
| BUSINESS | 默认 | HDF5 L1B 文件内置的业务定标系数 |
| RECALI | 运行时自动 | `../fy3d_recali/YYYYMM/RAD_YYYYMMDD.h5` |

RECALI 模式下，系统自动将 `recali` 参数映射为观测日期对应的 `YYYYMM` 目录，搜索再定标系数文件。若找不到对应月份的数据，自动回退到 BUSINESS 模式。

## 数据目录结构

### 输入数据

```
/data/Data_yuq/mersi/{YYYYMMDD}/
  FY3D_MERSI_GBAL_L1_{date}_{HHMM}_GEO1K_MS.HDF    ← 几何定位文件
  FY3D_MERSI_GBAL_L1_{date}_{HHMM}_1000M_MS.HDF    ← L1B 辐亮度文件

/data/nwp/{YYYYMMDD}/ORG/
  gfs.t{CC}z.pgrb2.0p25.f{LLL}                     ← GFS 0.25° 预报场 GRIB2

/data/Data_minmin/oisst/
  sst.day.mean.{YYYYMMDD}.hdf5                      ← OISST 逐日海温

../fy3d_recali/{YYYYMM}/
  RAD_{YYYYMMDD}.h5                                 ← 再定标系数文件
```

### 输出数据

```
/data/Data_yuq/fy3_cloud/{YYYYMMDD}/
  FY3D_MERSI_ORBT_L2_CLM_MLT_NUL_{date}_{HHMM}_1000M_MS_{CAL}.HDF
```

每个 HDF5 文件包含：

| 数据集 | 说明 |
|--------|------|
| `Cloud_Mask` | 16 级云掩膜比特阵列（2000×2048×6 bytes） |
| `Quality_Assurance` | 逐像素 QA 质量标记 |

### 常用测试日期

测试基准日期为 `20220803`，该日期对应的 MYD35（MODIS 官方云产品）为标准验证集。

该日期有 3 个观测时次：`0740`、`0830`、`0920`（UTC）。

## 运行时间参考

| 场景 | 单时次 | 3 时次串行 | 3 时次并行（3核） |
|------|--------|-----------|-----------------|
| BUSINESS | ~70s | ~3.5min | ~1.2min |
| RECALI | ~70s | ~3.5min | ~1.2min |
| 双定标合计 | — | ~7min | ~1.2min |

以上时间基于 Intel Xeon 处理器测试，实际运行时间因硬件而异。

## 验证

### 自动验证

```bash
python run_fylat.py --date 20220803 --verify-only
```

脚本会检查所有输出文件是否存在且大小合理（>1 MB）。

### 手动验证

```bash
# 对比 BUSINESS 和 RECALI 的差异
h5diff \
  /data/Data_yuq/fy3_cloud/20220803/FY3D_MERSI_ORBT_L2_CLM_MLT_NUL_20220803_0740_1000M_MS_BUSINESS.HDF \
  /data/Data_yuq/fy3_cloud/20220803/FY3D_MERSI_ORBT_L2_CLM_MLT_NUL_20220803_0740_1000M_MS_RECALI.HDF

# 检查可执行文件 MD5
md5sum fylat_FY3_MERSI_II_PGS
# 基线: 937af04506a1ffbb9fcf7a78579e11d7
```

BUSINESS 和 RECALI 之间应有约 189 万像素的云掩膜差异（正常，因反射波段定标系数不同）。

### 迁移回归

固定三轨、两套定标的可重复比较：

```bash
python scripts/compare_clm_regression.py \
  --date 20220803 \
  --times 0715,0720,0740 \
  --calibrations business,recali \
  --output-root /data/Data_yuq/fy3_cloud
```

如果正在验证 M1/M2 plumbing 是否完全不改变输出，应准备一份旧版本输出目录，再做同定标 old/new 对比：

```bash
python scripts/compare_clm_regression.py \
  --date 20220803 \
  --times 0715,0720,0740 \
  --calibrations business,recali \
  --baseline-root /data/Data_yuq/fy3_cloud_baseline \
  --output-root /data/Data_yuq/fy3_cloud \
  --require-identical
```

`--require-identical` 会在任一 SHA256 或 dataset 差异出现时返回非零退出码。

### MYD35 精度与分层诊断

```bash
python analyze_accuracy.py --date 20220803 --output accuracy_20220803.json
```

默认输出 Overall、纬度带、day/night、SZA 分箱、BT11 温度分箱、海岸/内陆水体、NWP forecast age 和重点光谱特征诊断。若只需要快速二分类指标，可跳过 L1B 光谱统计：

```bash
python analyze_accuracy.py --date 20220803 --no-spectral-diagnostics
```

## 常见问题

### 编译失败

1. 确认 `conda activate cloudmask` 已执行
2. 确认 `ifort`、`icc`、`h5fc` 在 PATH 中
3. 确认 `/opt/intel_lib/` 下 HDF5/HDF4 库完整

### NWP 数据找不到

确认 `/data/nwp/{date}/ORG/` 下有 GFS GRIB2 文件。批处理入口会先用 Python ecCodes 预生成 Fortran 兼容二进制缓存文件（约 1.1 GB/文件），存储在同级目录下；只有显式设置 `FYLAT_NWP_BACKEND=wgrib2` 时才走旧 wgrib2 兼容路径。

### NISE 雪冰掩膜报错

确认 `coeff/sfc_snow_ice/` 目录下有对应月份的 `NISE_SSMIF13_EASEGRID_M{MM}.HDF` 文件。程序读取时依赖 `code_root_path` 末尾的 `/`，新版本已自动处理。

### 并行运行时内存不足

每个 Fortran 进程约占用 3-4 GB 内存（主要来自 NWP 1.1 GB × 2 + 图像数据）。建议核心数不超过可用内存 / 4 GB。
