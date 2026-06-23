# FYLAT — FY-3D MERSI-II Cloud Mask Retrieval System V3.1

## 工程概述

本项目是 **FYLAT (FengYun Land-Atmosphere Toolkit)** 的风云三号D星（FY-3D）MERSI-II 中分辨率光谱成像仪的云检测反演系统，由中国气象局国家卫星气象中心（NSMC/CMA）Min Min 团队开发。

系统读取 FY-3D/MERSI-II 卫星 L1 级 HDF5 观测数据，结合数值天气预报（NWP）和辅助数据，通过辐射传输模型（RTM）和多光谱决策树算法，生成 L2 级云掩膜产品（Cloud Mask, CLM），并可扩展输出云量（CLA）、云相态（CLP）、云顶高度（CTP）、云光学厚度（COT）和海表温度（SST）产品。

## 技术架构

```
retrieval_system_V3.1_cldmask/
├── src/                        # Fortran/C 核心源码（~37,000 行）
│   ├── fylat_FY3_MERSI_II_PGS_Driver.f90  # 主驱动程序
│   ├── cloudmask/              # 云掩膜算法核心（~60+ 个文件）
│   ├── cloudamount/            # 云量算法
│   ├── cloudphase/             # 云相态算法（已注释/未启用）
│   ├── cloudheight/            # 云顶高度算法（已注释/未启用）
│   ├── cloudod_day/            # 云光学厚度算法（已注释/未启用）
│   ├── sea_surface_temperature/ # SST 算法（已注释/未启用）
│   └── *.f90, *.f, *.c         # 支撑模块（IO、RTM、数值计算等）
├── coeff/                      # 辅助数据系数文件
│   ├── fylat_thresholds.*      # 云检测阈值表
│   ├── fylat_ecosystem.hdf     # 生态系统地图
│   └── sfc_emiss/, sfc_snow_ice/, sfc_albedo/  # 地表参数
├── include/                    # 头文件
│   └── cloudmask.inc           # Fortran include 定义
├── wgrib/                      # NWP GRIB2→二进制转换工具
│   ├── wgrib, wgrib2           # 二进制工具
│   └── *.sh                    # 多套 NWP 数据转换脚本
├── drivefile/                  # 驱动配置文件
│   ├── *.nml                   # Namelist 配置文件（多场景）
│   └── *.py                    # Python 批处理驱动
├── convert_oisst/              # OISST 海温数据 NC→HDF5 转换
├── fylat_FY3_MERSI_II_PGS      # 编译后的可执行文件
├── paral_bat_driver_mersi_ii_fylat.py  # 并行批处理主驱动
├── build.sh / clean.sh         # 编译/清理脚本
├── delete_file.py              # 文件清理工具
└── VIS_Cal_Coeff.xcfg          # 可见光定标系数
```

## 处理流程

```text
STEP 0: 读取 namelist 配置文件 (.nml)
  ├── 指定传感器 ID（21=FY-3D, 22=FY-3E）
  ├── 指定 NWP 数据源（GFS/NCEP/T639/GRAPES/GDAS 等 10 种）
  ├── 指定启用算法（cloudmask/cloudamount/sst 等）
  └── 指定输入输出文件路径

STEP 1: 读取 FY-3D MERSI-II GEO + L1b HDF5 数据
  ├── 经纬度、太阳/传感器天顶角/方位角
  ├── 19 个可见光/近红外通道反射率（0.41–2.13 μm）
  ├── 6 个红外通道亮温（3.8, 4.05, 7.3, 8.5, 11.0, 12.0 μm）
  └── 计算散射角、太阳-地球距离

STEP 2: 读取 NWP 数据并插值到卫星观测时间
  ├── GRIB2→二进制转换（wgrib2）
  ├── 时间插值（相邻两个 NWP 时次）
  └── 空间插值到 101 层气压层

STEP 3: 读取辅助数据集
  ├── 雪/冰掩膜 (NISE)
  ├── 生态系统类型 (IGBP)
  ├── OISST 海表温度
  ├── 地表发射率（6 个 IR 通道）
  ├── 地表反照率（5 个 NIR/SWIR 通道）
  └── 运行辐射传输模型（PFAAST），计算各 IR 通道的晴空亮温

STEP 4: 运行云掩膜算法
  ├── 逐像素决策树分类
  │   ├── 按地表类型分：陆地/水体/海岸/沙漠/雪冰
  │   ├── 按区域分：极区（>60°）/非极区/南极（<-60°）
  │   ├── 按太阳天顶角分：白天/夜间（SZA>85°）
  │   └── 执行对应的光谱测试组合
  ├── 辅助检测：阴影/非云障碍物/薄卷云
  ├── 设置质量控制位（QA bits）
  └── 生成云掩膜比特阵列（6 字节/像素，含置信度）

STEP 5: 写入 L2 输出产品（HDF5 格式）
  └── 可选：云量（5km）、SST 等扩展产品

STEP 6: 释放内存
```

## 算法详解

### 云掩膜（Cloud Mask）核心算法

算法源自 MODIS MOD35 云检测框架，针对 FY-3D MERSI-II 25 通道特性做了适配改造。

**决策树分叉逻辑（polar_module.f90 / land_module.f90 / water_module.f90）：**

| 地表类型 | 区域 | 光照 | 处理模块 |
|---------|------|------|---------|
| 极区 | 陆地 | 白天 | PolarDay_land, PolarDay_desert_c, PolarDay_desert, PolarDay_coast, PolarDay_snow |
| 极区 | 陆地 | 夜间 | PolarNite_land, PolarNite_snow |
| 极区 | 海洋 | 白天 | PolarDay_ocean |
| 极区 | 海洋 | 夜间 | PolarNite_ocean |
| 极区 | 雪盖 | 夜间 | PolarNite_snow |
| 南极 | — | 白天 | Antarctic_day |
| 非极区 | 陆地 | 白天 | LandDay (含沙漠/海岸/雪盖子类型) |
| 非极区 | 陆地 | 夜间 | LandNite (含沙漠/海岸/雪盖子类型) |
| 非极区 | 海洋 | 白天 | ocean_day |
| 非极区 | 海洋 | 夜间 | ocean_nite |
| — | 雪盖 | 白天 | Day_snow |
| — | 雪盖 | 夜间 | Nite_snow |

**主要光谱测试（从阈值 include 文件中提取）：**
- 可见光反射率测试（0.66/0.86/1.38 μm 等）
- NDSI 雪检测（0.55 μm vs 1.64/2.13 μm）
- 红外亮温差测试（BTD 8.5-11, 3.8-11, 11-12）
- 可见光比率测试（VRAT = R0.86/R0.65）
- 空间均匀性测试（3x3 滑动窗口）
- 薄卷云检测（1.38 μm 反射率 + IR BTD）
- 太阳耀斑检测（glint angle ≤ 36°）

**云置信度输出（4 级）：**
- 0: 有云（Cloudy）
- 1: 可能云（Probably Cloudy）
- 2: 可能晴空（Probably Clear）
- 3: 晴空（Clear）

### 支持的 NWP 数据源（10 种）

| 编号 | 数据源 | 分辨率 | 层数 | GRIB 版本 |
|------|--------|--------|------|-----------|
| 1 | NCEP Reanalysis | 1.0° | 26 | GRIB1 |
| 2 | GFS | 1.0° | 26 | GRIB2 |
| 3 | T639 | 0.125° | 36 | GRIB2 |
| 4 | NCEP Reanalysis | 1.0° | 26 | GRIB2 |
| 5 | GFS | 0.5° | 26 | GRIB2 |
| 6 | GRAPES GFS | 0.25° | 40 | GRIB2 |
| 7 | GDAS | 0.25° | 31 | GRIB2 |
| 8 | GFS | 0.25° | 31 | GRIB2 |
| 9 | GFS | 0.5° | 41 | GRIB2 |
| 10 | GFS | 0.25° | 41 | GRIB2 |

## 技术栈

| 组件 | 技术 |
|------|------|
| 核心算法 | Fortran 90/77 |
| 编译 | Intel Fortran (ifort) + Intel C (icc) |
| 数据格式 | HDF5, HDF4, GRIB2, NetCDF |
| 辐射传输 | PFAAST (Radiative Transfer Model) |
| Python 驱动 | Python 3 (multiprocessing) |
| NWP 转换 | wgrib/wgrib2 shell 脚本 |
| 运行环境 | Linux x86_64, Intel + HDF5 库链 |

## 构建与运行

### 编译
```bash
# 修改 src/fylat_makefile_cldmask 中的路径
#   root_path = 项目 src 路径
#   HDF4MYLIB / HDF5MYLIB = HDF4/HDF5 库路径
./build.sh
```

### 运行（单次处理）
```bash
./fylat_FY3_MERSI_II_PGS config.nml
```

### 运行（批量处理）
```bash
# 修改 paral_bat_driver_mersi_ii_fylat.py 中的路径和时间范围
python3 paral_bat_driver_mersi_ii_fylat.py
```
支持的传感器：
- `sensor_id = 1`:  MODIS → MERSI-II 格式转换
- `sensor_id = 2`:  MODIS → MERSI-II（第二种转换）
- `sensor_id = 3`:  VIIRS → MERSI-II 格式转换
- `sensor_id = 21`: FY-3D 真实 MERSI-II
- `sensor_id = 22`: FY-3E 真实 MERSI-II

## 当前状态

### 已启用功能
- 云掩膜（Cloud Mask）算法 — **主要运行的算法**
- 云量（Cloud Amount）算法 — 5km 网格
- SST 反演 — 可选
- 中间结果输出 — 调试用

### 已注释/未启用功能
- 云相态（Cloud Phase）算法 — 被注释
- 云顶高度（Cloud Top Height）算法 — 被注释
- 云光学厚度白天/夜间（Cloud Optical Depth）算法 — 被注释
- 云类型 II 算法 — 被注释

## 可优化方向

### 1. 代码架构与组织

- **重复的 NWP 数据结构**：`data_arrays_module.f90` 中定义了 7 组几乎相同的 NWP 插值类型（nwp26/nwp31/nwp36/nwp40/nwp41/T639/grapes），每组 ~25 个成员变量，仅分辨率和层数不同。可使用参数化派生类型统一为一个模板，减少 ~500 行重复定义。
- **重复的 NWP deallocate 逻辑**：`fylat_FY3_MERSI_II_PGS_Driver.f90` 中每类 NWP 有一个 deallocate 子程序，结构和逻辑高度雷同，可统一处理。
- **Python 驱动中 6 个几乎相同的 `convert_nwp` 函数**：`paral_bat_driver_mersi_ii_fylat.py`（约 470 行）中 `convert_nwp` 到 `convert_nwp6` 逻辑雷同，差异仅在文件名模板和输出格式，可合并为一个带参数的函数。
- **Python 中 5 个 `find_nwp_name` 变体**：逻辑几乎相同，仅查找的文件名模式不同，可统一。

### 2. 性能优化

- **主循环可做 OpenMP 并行化**：云掩膜算法的主循环 `line_loop_1`/`element_loop_1` 逐像素处理 2048×2000 图像，当前为串行。各像素的决策树判断在计算上独立（除 3×3 窗口均匀性检测），可用 OpenMP 指令实现线程级并行加速 4–8 倍。
- **Python 驱动的多进程利用不充分**：`paral_bat_driver_mersi_ii_fylat.py` 中 `multiprocessing.Pool(nthread)` 默认 `nthread=1`，未利用多核。
- **3×3 邻域数据读取效率**：`get_pxl3X3` 逐像素从大数组中提取 3×3 邻域，边界判断逻辑简单但每个像素都重复检查。可预定义边界掩膜或使用滑动窗口缓存。
- **阈值文件每次运行时重新读取**：`thresholds_read()` 在主循环外只调用一次，这点是正确的。但如果一天处理多轨数据，每次启动可执行文件都重新读取相同的阈值文件，可考虑将阈值数据预加载并持久化。

### 3. 内存管理

- **存在注释掉的分配/释放不对称**：多个算法的 allocate/deallocate 被注释掉（cloudphase/cloudheight/cloudod_day/sst），如果未来启用需要确认内存管理无泄漏。
- **cm_tmp 数组的分配时机**：`cm_tmp` 仅在 `cloudmask_id==1` 且其他算法启用时才分配，但释放逻辑条件匹配，需确保所有条件下路径一致。
- **io_module.f90 中的临时缓冲区**：每次读取 HDF5 数据都分配/释放 `var_int1/var_int4/var_real4/var_char` 临时数组，可考虑复用。

### 4. 可维护性

- **混合使用 Fortran 77 和 Fortran 90**：文件中 `.f` 和 `.f90` 混用，Fortran 77 子程序使用固定格式，缺少 `implicit none`，变量命名较短且无类型声明，增加了维护难度。可逐步迁移到 Fortran 90 自由格式，并统一引入 `implicit none`。
- **硬编码路径和魔术数字**：
  - Python 驱动中 `L1_path`、`GEO_path`、`code_path` 等硬编码，可通过命令行参数或配置文件管理。
  - `temp_fy3d_config_20200330_0000.nml` 和 `.xcfg` 文件中的路径硬编码为特定用户目录。
  - `platform_module.f90` 中的传感器通道参数（中心波长、波数等）以数组形式硬编码在 case 分支中，新增传感器需要重新编译。
- **大量注释掉的代码**：在多个文件中存在成块注释掉的代码（如 `io_module.f90` 中的 MODIS 转换分支、`platform_module.f90` 中的定标系数、`fylat_FY3_MERSI_II_PGS_Driver.f90` 中的多个算法调用），建议清理或通过编译选项控制。
- **调试标记分散**：代码中存在 `lyj`、`jincheng test` 等个人调试标记（如 cloudmask 模块中 13 个 `out_*` 数组的分配和赋值），不应出现在生产代码中。
- **缺少版本控制**：项目无 `.git` 目录，无法追踪变更历史，也无 CI/CD 配置。

### 5. 安全性

- **Python 中的 `os.system()` 调用**：`paral_bat_driver_mersi_ii_fylat.py` 大量使用 `os.system(ww + ' ' + filename + ' ' + bin_name)` 拼接 shell 命令，虽 filename 来自 glob 结果而非用户输入，但更好的做法是使用 `subprocess.run()` 带参数列表调用。
- **Fortran 中的数组越界风险**：部分 `.f` 文件没有 `implicit none`，变量类型可能被隐式推导，且数组边界检查依赖编译选项（当前 Makefile 中已注释的 `-check all` 调试编译选项）。

### 6. 可扩展性

- **NWP 类型扩展困难**：每新增一种 NWP 需要修改至少 5 处代码（data_arrays_module、read_nwp_data_module、nwp_utils_module、Driver、Makefile），建议抽象 NWP 数据源为统一接口。
- **传感器扩展困难**：`platform_module.f90` 每个新传感器需新增一个 case 分支并硬编码参数，传感器参数可外部化为配置文件。
- **算法模块耦合度高**：各算法模块虽通过 `use` 语句做了模块化，但共享全局状态（`sat`, `geo`, `nwp` 等类型变量），增加了独立单元测试的难度。

### 7. 编译配置

- **Makefile 中的路径硬编码**：`fylat_makefile_cldmask` 中 `root_path`、HDF 库路径均为绝对路径，移植到其他机器需修改 Makefile。
- **依赖 Intel 编译器**：完全依赖 ifort/icc 编译器链，未测试 GNU gfortran/gcc 编译兼容性。
- **缺少自动化构建系统**：没有 CMake/autotools 支持，依赖手动修改 Makefile。

### 8. 文档

- **原始 readme 仅 3 行**：只说明了命令行调用方式，无项目介绍、算法描述、依赖说明。
- **Fortran 源码中有详细注释头**：每个模块/子程序有作者、输入输出参数说明，但缺少整体算法流程文档和阈值物理意义说明。

## 输出产品文件说明

| 文件名缩写 | 产品 | 说明 |
|-----------|------|------|
| CLM | Cloud Mask | 1km 云掩膜（6 字节比特阵列 + 10 字节 QA） |
| CLA | Cloud Amount | 5km 云量 |
| CLP | Cloud Phase | 云相态和类型 |
| CTP | Cloud Top Properties | 云顶温度/气压/高度 |
| COT | Cloud Optical Thickness | 云光学厚度（白天） |
| CON | Cloud Optical Properties (Night) | 云光学特性（夜间） |
| SST | Sea Surface Temperature | 海表温度 |
| INTERMED | Intermediate Results | 中间调试输出 |

## 作者

- 负责人：Min Min (minmin@cma.gov.cn)
- 合作者：Wu Xiao, Zheng Zhaojun, Liu Ruixia, Zhang Miao, Yang Changjun, Qiu Hong
- 单位：中国气象局国家卫星气象中心（NSMC/CMA）
