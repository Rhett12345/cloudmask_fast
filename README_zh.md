# FYLAT — 风云三号D星 MERSI-II 云检测反演系统

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Version](https://img.shields.io/badge/version-0.1.0-green.svg)](https://github.com/Rhett12345/cloudmask_fast)

**FYLAT (FengYun Land-Atmosphere Toolkit)** 是风云三号D星（FY-3D）MERSI-II 中分辨率光谱成像仪的云检测反演系统。系统读取 L1 级 HDF5 卫星观测数据，结合数值天气预报（NWP）和辅助数据，通过多光谱决策树算法（源自 MODIS MOD35 框架）生成 L2 级云掩膜产品。

**由中山大学闵敏团队开发维护。**

## 处理流程

```
FY-3D MERSI-II L1 HDF5 ──┬── NWP (GRIB2) 插值
                          ├── 辅助数据（发射率、反照率、雪冰、生态系统）
                          ├── 辐射传输模型 (PFAAST)
                          └── 决策树云检测 ──→ L2 CLM/CLA/SST (HDF5)
```

1. **读取 namelist 配置** (`.nml`) — 传感器 ID、NWP 数据源、算法开关、I/O 路径
2. **读取 L1 数据** — 25 个 MERSI-II 通道（19 可见光/近红外 + 6 红外）、几何信息
3. **读取 NWP 数据**（GRIB2 → 二进制，通过 wgrib2）— 空间插值到 101 层气压层，时间插值于两个预报时次之间
4. **读取辅助数据** — 雪/冰掩膜 (NISE)、生态系统类型 (IGBP)、OISST 海温、地表发射率/反照率 → 运行 PFAAST 辐射传输模型计算晴空亮温
5. **运行云掩膜算法** — 逐像素决策树分类（极区/陆地/水体三大模块，约 20 种场景子程序），3×3 空间均匀性检测，薄卷云/阴影/太阳耀斑检测，QA 位打包
6. **输出 L2 HDF5** — 云掩膜（1 km，6 字节比特阵列 + QA），可选云量（5 km）和 SST

### 支持的传感器

| 传感器 ID | 说明 |
|-----------|------|
| 21 | FY-3D MERSI-II |
| 22 | FY-3E MERSI-II |
| 1–3 | MODIS / VIIRS → MERSI-II 格式转换 |

### 支持的 NWP 数据源（10 种）

GFS（多种分辨率）、NCEP 再分析、T639、GRAPES GFS、GDAS — GRIB1/GRIB2，26–41 层。

## 快速开始

### 环境依赖

- Linux x86_64
- Intel Fortran (ifort) + Intel C (icc)
- HDF5 / HDF4 库
- wgrib / wgrib2
- Python 3（批量处理用）

### 编译

```bash
conda activate cloudmask
./build.sh              # 执行: cd src && make -f fylat_makefile_cldmask
```

### 运行（单次处理）

```bash
./fylat_FY3_MERSI_II_PGS config.nml
```

### 运行（批量处理）

```bash
python3 paral_bat_driver_mersi_ii_fylat.py
```

## 项目结构

```
├── src/                    # Fortran/C 核心源码（~37,000 行）
│   ├── fylat_FY3_MERSI_II_PGS_Driver.f90  # 主驱动程序
│   ├── cloudmask/          # 云掩膜算法核心
│   ├── cloudamount/        # 云量算法
│   └── *.f90, *.f, *.c     # 支撑模块（IO、RTM、数值计算等）
├── coeff/                  # 辅助系数文件（阈值表、RTM 查找表）
├── include/                # Fortran include 头文件
├── wgrib/                  # NWP GRIB2 → 二进制转换脚本
├── drivefile/              # Namelist 配置文件
├── convert_oisst/          # OISST 海温 NC → HDF5 转换
└── paral_bat_driver_mersi_ii_fylat.py  # 并行批处理驱动
```

## 云掩膜算法

算法源自 MODIS MOD35 云检测框架，针对 FY-3D MERSI-II 25 通道特性做了适配改造。

### 决策树分叉逻辑

| 地表类型 | 区域 | 光照 | 处理模块 |
|---------|------|------|---------|
| 陆地/沙漠/海岸/雪盖 | 非极区 | 白天 | LandDay, LandDay_desert, LandDay_coast |
| 陆地/沙漠/海岸/雪盖 | 非极区 | 夜间 | LandNite |
| 海洋 | 非极区 | 白天/夜间 | ocean_day, ocean_nite |
| 陆地/海洋/雪盖 | 极区 (>60°) | 白天/夜间 | PolarDay_*, PolarNite_* |
| — | 南极 (<-60°) | 白天 | Antarctic_day |

### 主要光谱测试

- 可见光反射率阈值测试（0.66/0.86/1.38 μm）
- NDSI 雪检测（0.55 μm vs 1.64/2.13 μm）
- 红外亮温差测试（BTD 8.5–11, 3.8–11, 11–12）
- 可见光比率测试（VRAT = R0.86/R0.65）
- 3×3 空间均匀性测试
- 薄卷云检测（1.38 μm 反射率 + IR BTD）
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
- **开发者**: 余强 (Yu Qiang)，中山大学 — yuqiang6@mail2.sysu.edu.cn, 2942204121@qq.com

## 开源许可

MIT 许可证 — 详见 [LICENSE](LICENSE)。

## 引用

如果您在研究中使用了本软件，请引用：

```
Min, M., & Yu, Q. (2025). FYLAT: FY-3D MERSI-II Cloud Mask Retrieval System (Version 3.1).
Sun Yat-sen University. https://github.com/Rhett12345/cloudmask_fast
```
