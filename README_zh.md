# FYLAT — 风云三号D星 MERSI-II 云检测反演系统

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Version](https://img.shields.io/badge/version-2.3.2-green.svg)](https://github.com/Rhett12345/cloudmask_fast)

**FYLAT (FengYun Land-Atmosphere Toolkit)** 是风云三号D/E星 MERSI-II 中分辨率光谱成像仪的云检测反演系统。系统读取 L1 级 HDF5 卫星观测数据，结合数值天气预报（NWP）和辅助数据，通过多光谱决策树算法（源自 MODIS MOD35 框架）生成 L2 级云掩膜（CLM）产品。

由**中山大学闵敏团队**开发维护。

## 快速开始

### 环境依赖

- Linux x86_64
- Intel Fortran (`ifort`) + Intel C (`icc`)
- HDF5 / HDF4 库
- wgrib2
- Python 3, conda 虚拟环境 `cloudmask`

### 一键运行

```bash
conda activate cloudmask
python run_fylat.py --date 20220803
```

一条命令自动完成：
1. 扫描 L1B 目录，发现所有观测时次
2. 自动发现 NWP GRIB2 预报文件
3. 编译 Fortran 可执行文件（已编译则跳过）
4. 并行运行**业务定标**（BUSINESS）和**再定标**（RECALI）反演
5. 验证全部输出文件

```bash
python run_fylat.py --date 20220803 --cores 4    # 指定 4 核并行
python run_fylat.py --date 20220803 --dry-run    # 仅预览配置，不运行
python run_fylat.py --date 20220803 --skip-build # 跳过编译
python run_fylat.py --date 20220803 --verify-only  # 仅验证已有输出
```

### 定标模式

| 模式 | 说明 |
|------|------|
| `business` | 使用 HDF5 文件内置的业务定标系数 |
| `recali` | 使用外部再定标系数 `../fy3d_recali/YYYYMM/RAD_*.h5`（按日期自动加载） |

### 手动编译运行

```bash
conda activate cloudmask
./build.sh                                    # 编译
./fylat_FY3_MERSI_II_PGS config.nml          # 单次运行
python scripts/run_single.py --date 20220803 --time 0740 --calibration business
```

## 处理流程

```
FY-3D MERSI-II L1 HDF5 ──┬── NWP (GRIB2) 插值
                          ├── 辅助数据（发射率、反照率、雪冰、生态系统）
                          ├── 辐射传输模型 (PFAAST)
                          └── 决策树云检测 ──→ L2 CLM (HDF5)
```

1. 读取 namelist 配置 (`.nml`) — 传感器 ID、NWP 数据源、算法开关、I/O 路径
2. 读取 L1 数据 — 25 个 MERSI-II 通道（19 可见光/近红外 + 6 红外）、几何信息
3. 读取 NWP 数据（GRIB2 → 二进制，通过 wgrib2）— 空间插值到 101 层气压，时间插值于两预报时次间
4. 读取辅助数据 — 雪/冰掩膜 (NISE)、生态系统 (IGBP)、OISST 海温、发射率/反照率 → PFAAST 辐射传输模式计算晴空亮温
5. 运行云掩膜算法 — 逐像素决策树，3×3 空间均匀性检测，薄卷云/阴影/太阳耀斑检测，QA 位打包
6. 输出 L2 HDF5 — 云掩膜（1 km，6 字节比特阵列 + QA）

## 项目结构

```
├── src/                         # Fortran/C 核心源码
│   ├── fylat_FY3_MERSI_II_PGS_Driver.f90   # 主驱动程序
│   ├── fylat_makefile_cldmask              # Makefile
│   ├── cloudmask/               # 云掩膜算法（~20 个场景模块）
│   ├── cloudamount/             # 云量算法（默认未启用）
│   ├── cloudphase/              # 云相态（默认未启用）
│   ├── cloudheight/             # 云顶高度（默认未启用）
│   ├── cloudod_day/             # 云光学厚度（默认未启用）
│   ├── sea_surface_temperature/ # 海表温度（默认未启用）
│   └── *.f90, *.f, *.c         # 支撑模块：IO、RTM、数值计算、平台配置
├── coeff/                       # 辅助系数文件（阈值表、RTM 查找表、辅助数据）
├── wgrib/                       # NWP GRIB2 → 二进制转换脚本
├── python/fylat/                # Python 配置与定标管理
│   ├── config.py                # YAML → .nml namelist 生成
│   └── calibration.py           # 再定标系数发现与加载
├── config/                      # YAML 配置
│   ├── default.yaml             # 默认配置
│   └── scenes/                  # 场景专属覆盖配置
├── scripts/
│   └── run_single.py            # 单场景运行器（基于 YAML）
├── visualize/                   # 验证与可视化工具
├── build.sh                     # 编译脚本
├── clean.sh                     # 清理编译产物
└── run_fylat.py                 # 一键全链条自动化反演
```

## 云掩膜算法

算法源自 MODIS MOD35 云检测框架，针对 FY-3D MERSI-II 25 通道特性适配。

### 决策树逻辑

| 地表类型 | 区域 | 光照 | 处理模块 |
|---------|------|------|---------|
| 陆地/沙漠/海岸/雪盖 | 非极区 | 白天 | LandDay, LandDay_desert, LandDay_coast |
| 陆地/沙漠/海岸/雪盖 | 非极区 | 夜间 | LandNite |
| 海洋 | 非极区 | 白天/夜间 | ocean_day, ocean_nite |
| 陆地/海洋/雪盖 | 极区 (>60°) | 白天/夜间 | PolarDay_*, PolarNite_* |
| — | 南极 (<-60°) | 白天 | Antarctic_day |

### 主要光谱测试

- 可见光反射率阈值（0.66/0.86/1.38 μm）
- NDSI 雪检测（0.55 μm vs 1.64/2.13 μm）
- 红外亮温差（BTD 8.5–11, 3.8–11, 11–12 μm）
- 可见光比率测试（R0.86/R0.65）
- 3×3 空间均匀性测试
- 薄卷云检测（1.38 μm + IR BTD）
- 太阳耀斑检测（glint angle ≤ 36°）

### 云置信度

| 等级 | 含义 |
|------|------|
| 0 | 有云 (Cloudy) |
| 1 | 可能云 (Probably Cloudy) |
| 2 | 可能晴空 (Probably Clear) |
| 3 | 晴空 (Clear) |

## 团队

- **负责人**: 闵敏 (Min Min)，中山大学
- **开发者**: 余强 (Yu Qiang)，中山大学 — yuqiang6@mail2.sysu.edu.cn

## 开源许可

MIT 许可证 — 详见 [LICENSE](LICENSE)。

## 引用

```
Min, M., & Yu, Q. (2025). FYLAT: FY-3D MERSI-II Cloud Mask Retrieval System (Version 3.1).
Sun Yat-sen University. https://github.com/Rhett12345/cloudmask_fast
```
