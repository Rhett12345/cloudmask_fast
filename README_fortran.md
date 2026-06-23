# FYLAT — FY-3D MERSI-II Cloud Mask Retrieval System V3.1

## 项目概述

**FYLAT (FengYun Land-Atmosphere Toolkit)** 是风云三号D星（FY-3D）MERSI-II 中分辨率光谱成像仪的云检测反演系统，由中国气象局国家卫星气象中心（NSMC/CMA）开发。

系统读取 FY-3D/MERSI-II 卫星 L1 级 HDF5 观测数据，结合数值天气预报（NWP）和辅助数据，通过辐射传输模型（RTM）和多光谱决策树算法，生成 L2 级云掩膜产品（Cloud Mask, CLM），并可扩展输出云量（CLA）、云相态（CLP）、云顶高度（CTP）、云光学厚度（COT）和海表温度（SST）产品。

## 目录结构

```
retrieval_system_V3.1_cldmask/
├── src/                            # Fortran/C 核心源码（~37,000 行）
│   ├── fylat_FY3_MERSI_II_PGS_Driver.f90    # 主驱动程序
│   ├── cloudmask/                  # 云掩膜算法核心（~60+ 个文件）
│   ├── cloudamount/                # 云量算法
│   ├── cloudphase/                 # 云相态算法（已注释/未启用）
│   ├── cloudheight/                # 云顶高度算法（已注释/未启用）
│   ├── cloudod_day/                # 云光学厚度算法（已注释/未启用）
│   ├── sea_surface_temperature/    # SST 算法（已注释/未启用）
│   └── *.f90, *.f, *.c             # 支撑模块（IO、RTM、数值计算等）
├── coeff/                          # 辅助数据系数文件
│   ├── fylat_thresholds.*          # 云检测阈值表
│   ├── fylat_ecosystem.hdf         # 生态系统地图
│   └── sfc_emiss/, sfc_snow_ice/, sfc_albedo/  # 地表参数
├── include/                        # 头文件
│   └── cloudmask.inc               # Fortran include 定义
├── wgrib/                          # NWP GRIB2→二进制转换工具
│   ├── wgrib, wgrib2               # 二进制工具
│   └── *.sh                        # 多套 NWP 数据转换脚本
├── drivefile/                      # 驱动配置文件
│   ├── *.nml                       # Namelist 配置文件（多场景）
│   └── *.py                        # Python 批处理驱动
├── convert_oisst/                  # OISST 海温数据 NC→HDF5 转换
├── fylat_FY3_MERSI_II_PGS          # 编译后的可执行文件
├── paral_bat_driver_mersi_ii_fylat.py  # 并行批处理主驱动
├── build.sh / clean.sh             # 编译/清理脚本
├── delete_file.py                  # 文件清理工具
└── VIS_Cal_Coeff.xcfg              # 可见光定标系数
```

## 运行环境

### 依赖环境

- **操作系统**: Linux x86_64
- **编译器**: Intel Fortran (ifort) + Intel C (icc)
- **Conda 环境**: `cloudmask`
- **数据格式库**: HDF5, HDF4, GRIB2, NetCDF
- **Python**: Python 3 (multiprocessing)

### 激活 Conda 环境

```bash
conda activate cloudmask
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

**决策树分叉逻辑：**

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

**主要光谱测试：**
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

## 构建与运行

### 编译

```bash
# 激活 conda 环境
conda activate cloudmask

# 进入项目目录
cd /home/liusy2020/yuq/cloudmask/retrieval_system_V3.1_cldmask

# 修改 src/fylat_makefile_cldmask 中的路径（如需要）
#   root_path = 项目 src 路径
#   HDF4MYLIB / HDF5MYLIB = HDF4/HDF5 库路径

# 编译
./build.sh
```

### 运行（单次处理）

```bash
# 激活 conda 环境
conda activate cloudmask

# 运行
./fylat_FY3_MERSI_II_PGS config.nml
```

### 运行（批量处理）

```bash
# 激活 conda 环境
conda activate cloudmask

# 修改 paral_bat_driver_mersi_ii_fylat.py 中的路径和时间范围
python3 paral_bat_driver_mersi_ii_fylat.py
```

### 支持的传感器

| sensor_id | 说明 |
|-----------|------|
| 1 | MODIS → MERSI-II 格式转换 |
| 2 | MODIS → MERSI-II（第二种转换） |
| 3 | VIIRS → MERSI-II 格式转换 |
| 21 | FY-3D 真实 MERSI-II |
| 22 | FY-3E 真实 MERSI-II |

## 配置文件说明

Namelist 配置文件示例 (`*.nml`):

```fortran
&config
  fylat_sensor_id     = 21,                    ! 传感器ID: 21=FY-3D, 22=FY-3E
  code_root_path      = "/path/to/code/",      ! 代码根目录
  L1b_data_path       = "/path/to/L1/",        ! L1数据路径
  nwp_data_path       = "/path/to/nwp/",       ! NWP数据路径
  oisst_data_path     = "/path/to/oisst/",     ! OISST数据路径
  fy3_mersi_GEO_data  = "/path/to/GEO.HDF",   ! GEO数据文件
  fy3_mersi_L1b_data  = "/path/to/L1b.HDF",   ! L1b数据文件
  fy3_mersi_CLM_data  = "/path/to/CLM.HDF",   ! 云掩膜输出
  fy3_mersi_CLA_data  = "/path/to/CLA.HDF",   ! 云量输出
  fy3_mersi_SST_data  = "/path/to/SST.HDF",   ! SST输出
  fylat_nwp_opt       = 5,                     ! NWP数据源选项
  cloudmask_id        = 1,                     ! 启用云掩膜算法
  surface_sst_id      = 1,                     ! 启用SST反演
/
```

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

## 当前状态

### 已启用功能
- 云掩膜（Cloud Mask）算法 — **主要运行的算法**
- 云量（Cloud Amount）算法 — 5km 网格
- SST 反演 — 可选
- 中间结果输出 — 调试用

### 已注释/未启用功能
- 云相态（Cloud Phase）算法
- 云顶高度（Cloud Top Height）算法
- 云光学厚度白天/夜间（Cloud Optical Depth）算法
- 云类型 II 算法

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

## 作者

- 负责人：Min Min (minmin@cma.gov.cn)
- 合作者：Wu Xiao, Zheng Zhaojun, Liu Ruixia, Zhang Miao, Yang Changjun, Qiu Hong
- 单位：中国气象局国家卫星气象中心（NSMC/CMA）

## 相关目录

源码位置：`/home/liusy2020/yuq/cloudmask/retrieval_system_V3.1_cldmask/`
