# FYLAT 工程化改造计划

> **项目**：cloudmask_fast — FY-3D MERSI-II Cloud Mask Retrieval System  
> **目标**：Fortran 保留计算核心，C++ 重写 I/O 与业务层，Python 负责调度与分析  
> **文档版本**：v1.0 | 2025-07

---

## 目录

1. [现状评估](#1-现状评估)
2. [目标架构](#2-目标架构)
3. [阶段一：基础工程化（低风险）](#3-阶段一基础工程化低风险)
4. [阶段二：C++ 胶水层建立](#4-阶段二c-胶水层建立)
5. [阶段三：Python 调度层升级](#5-阶段三python-调度层升级)
6. [阶段四：验证与收尾](#6-阶段四验证与收尾)
7. [技术决策说明](#7-技术决策说明)
8. [里程碑与时间线](#8-里程碑与时间线)
9. [风险与缓解措施](#9-风险与缓解措施)

---

## 1. 现状评估

### 1.1 语言占比

| 语言 | 占比 | 行数（估算） | 主要职责 |
|------|------|-------------|---------|
| Fortran | 84.8% | ~37,000 行 | 全部逻辑（I/O、计算、业务） |
| Python | 10.0% | ~4,000 行 | 批量调度脚本 |
| Shell | 2.3% | ~900 行 | 构建与辅助脚本 |
| C | 1.0% | ~400 行 | 少量 C 辅助模块 |

### 1.2 已知工程问题清单

| 问题 | 位置 | 严重程度 | 备注 |
|------|------|---------|------|
| `clean.sh` 引用错误 makefile 名 | `clean.sh` | 高 | 应为 `fylat_makefile_cldmask` |
| HDF 库路径硬编码 | `fylat_makefile_cldmask` | 高 | `/opt/intel_lib/lib/` |
| 数据路径硬编码 | 各 `.nml` 文件 | 高 | 换机器或用户即失效 |
| 无构建隔离 | makefile | 中 | 没有 out-of-source build |
| debug 变量未清理 | `cloudmask/` 各模块 | 中 | `lyj` 前缀变量、`jincheng test` 注释 |
| 注释掉的未完成模块 | makefile | 中 | `cloudamount/`、`cloudphase/` 等 |
| 无单元测试 | 全局 | 中 | 数值结果无回归基线 |
| 再定标系数无统一管理 | Python 批量脚本 | 高 | 手动对应日期目录，易出错 |
| F77 固定格式混用 | `src/*.f` | 低 | 缺 `implicit none`，调试困难 |
| 无版本标记基线 | 仓库 | 低 | 改造前应打 tag |

### 1.3 再定标系数现状

- 系数存放于 `../fy3d_recali/<yyyymmdd>/` 按日期分目录
- 太阳反射波段 7 个通道，每次需同时测试**业务定标**与**再定标**两套
- 当前靠 Python 脚本手写日期字符串对应，无自动索引，无缺失日期处理
- 验证基准：MYD35 数据集，测试日期基准 `20220803`

---

## 2. 目标架构

```
┌──────────────────────────────────────────────────────────────┐
│                     Python 调度层                             │
│  BatchScheduler | RecaliManager | Visualizer | CLI / Config   │
│  (pybind11 调用 C++ 扩展模块，不再 subprocess 调可执行文件)    │
└────────────────────────┬─────────────────────────────────────┘
                         │  pybind11 / ctypes
┌────────────────────────▼─────────────────────────────────────┐
│                     C++ 胶水/业务层                            │
│  IOManager | NWPInterp | AncilLoader | RecaliApplier          │
│  (HDF5/GRIB2 I/O、时空插值、辅助数据、系数管理)               │
└────────────────────────┬─────────────────────────────────────┘
                         │  ISO_C_BINDING / extern "C"
┌────────────────────────▼─────────────────────────────────────┐
│                  Fortran 计算核心（保留）                       │
│  cloudmask/ | rtm_tran_module | planck_module | numerical      │
│  (决策树云检测、PFAAST RTM、Planck 辐射、数值插值)             │
└──────────────────────────────────────────────────────────────┘
```

### 2.1 各层职责边界

**Fortran 层（只做计算，不碰 I/O）**

- 保留：`cloudmask/` 所有决策树模块（约 20 个场景子程序）
- 保留：`rtm_tran_module.f90`、`rtm_utils_module.f90`
- 保留：`planck_module.f90`、`numerical.f90`、`constant.f90`
- 暴露接口：通过 `ISO_C_BINDING` 导出 C 兼容函数签名
- 移除职责：HDF5 读写、NWP 文件解析、namelist 配置读取

**C++ 层（I/O + 业务逻辑 + 系数管理）**

- `IOManager`：读写 L1 HDF5、输出 L2 HDF5，使用 C++ HDF5 API
- `NWPInterp`：GRIB2 解析（调用 eccodes/wgrib2）、时空插值到 101 压力层
- `AncilLoader`：NISE、IGBP、OISST、emissivity/albedo 加载
- `RecaliManager`：再定标系数统一索引与加载（核心新增模块）
- `RTMWrapper`：封装对 Fortran RTM 的调用，处理内存布局转换
- 通过 `pybind11` 向 Python 暴露 C++ 接口

**Python 层（调度 + 分析 + 可视化）**

- `BatchScheduler`：多进程批量任务，自动对齐再定标日期
- `RecaliComparer`：同日期两套系数结果对比，输出统计报告
- `Visualizer`：云掩膜可视化，对比 MYD35 差异图
- `ConfigBuilder`：生成结构化配置，替代手写 `.nml`
- CLI 入口：`fylat run`、`fylat batch`、`fylat compare`

---

## 3. 阶段一：基础工程化（低风险）

**目标**：不改任何算法逻辑，只修工程层面的问题。可以在 1~2 周内完成。

### 3.1 打基线 Tag

```bash
git tag v0.1-fortran-baseline
git push origin v0.1-fortran-baseline
```

所有后续改造都在新分支上进行，基线随时可回退。

### 3.2 修复 clean.sh

```bash
# clean.sh 原内容（错误）
make clean -f fylat_makefile

# 修复为
cd src && make clean -f fylat_makefile_cldmask && cd ..
```

### 3.3 迁移到 CMake

替换手写 makefile，目标结构：

```
cmake/
  FindHDF5Fortran.cmake
  FindIntelFortran.cmake
CMakeLists.txt          ← 顶层
src/
  CMakeLists.txt        ← Fortran 核心静态库
  cpp/
    CMakeLists.txt      ← C++ 胶水层（阶段二填充）
  python/
    CMakeLists.txt      ← pybind11 扩展（阶段三填充）
```

顶层 `CMakeLists.txt` 骨架：

```cmake
cmake_minimum_required(VERSION 3.20)
project(fylat LANGUAGES Fortran C CXX)

set(CMAKE_Fortran_STANDARD 2008)
set(CMAKE_CXX_STANDARD 17)

find_package(HDF5 REQUIRED COMPONENTS Fortran C HL)
find_package(Python3 3.8 COMPONENTS Development Interpreter REQUIRED)

# 可选模块开关（替代 makefile 注释）
option(BUILD_CLOUD_AMOUNT  "Build cloud amount module"  OFF)
option(BUILD_CLOUD_PHASE   "Build cloud phase module"   OFF)
option(BUILD_CLOUD_HEIGHT  "Build cloud height module"  OFF)
option(BUILD_SST           "Build SST module"           OFF)

add_subdirectory(src)
add_subdirectory(src/cpp)
add_subdirectory(src/python)
```

### 3.4 提取硬编码路径

新建 `env.sh`（开发环境）和 CMake cache 变量（CI 环境）：

```bash
# env.sh
export FYLAT_HDF_ROOT=/opt/intel_lib
export FYLAT_DATA_ROOT=/data/Data_yuq
export FYLAT_NWP_ROOT=/data/nwp
export FYLAT_RECALI_ROOT=/data/Data_yuq/fy3d_recali
export FYLAT_OUTPUT_ROOT=/data/Data_yuq/fy3_cloud
```

makefile 过渡期可用 `$(FYLAT_HDF_ROOT)` 引用，CMake 迁移后统一用 `find_package`。

### 3.5 清理 debug 残留

```bash
# 扫描所有 debug 标记
grep -rn "lyj\|jincheng test\|out_lyj\|test_lyj" src/ > debug_markers.txt

# 人工确认后批量删除（不要机械替换，部分可能是有效变量名前缀）
```

处理原则：注释掉的 debug 输出语句直接删除；`out_*` 数组如果只用于 debug 输出则一并删除；保留带有物理意义的变量（即使前缀是 `lyj`）。

### 3.6 统一 implicit none

对 `src/*.f`（F77 固定格式）统一补充 `IMPLICIT NONE`：

```bash
# 找出缺少 implicit none 的 .f 文件
grep -rL "IMPLICIT NONE\|implicit none" src/*.f
```

逐文件添加后跑一次完整编译验证，对比 `20220803` 日期的输出 HDF5 文件 MD5，确保数值结果不变。

### 3.7 阶段一验收标准

- `./build.sh` 在全新 conda 环境下一次成功，无绝对路径依赖
- `clean.sh` 正常清理编译产物
- 对 `20220803` 的 CLM 输出与 v0.1 基线 byte-for-byte 一致
- `debug_markers.txt` 中的条目全部处理完毕

---

## 4. 阶段二：C++ 胶水层建立

**目标**：用 C++ 重写 I/O 和业务逻辑，Fortran 计算核心通过 `ISO_C_BINDING` 对接。核心新增 `RecaliManager`。

### 4.1 Fortran 端暴露 C 接口

在每个需要被 C++ 调用的 Fortran 模块中添加 `ISO_C_BINDING` wrapper：

```fortran
! src/cloudmask/cloudmask_c_api.f90
module cloudmask_c_api
  use iso_c_binding
  use cloudmask_module
  implicit none
contains

  ! 暴露给 C++ 的云掩膜入口
  subroutine run_cloudmask_c( &
      n_lines, n_pixels,      &
      sat_data, geo_data,     &
      nwp_data, ancil_data,   &
      clm_result, qa_result   &
  ) bind(C, name="run_cloudmask")
    integer(c_int), intent(in), value :: n_lines, n_pixels
    real(c_float),  intent(in)  :: sat_data(*)
    real(c_float),  intent(in)  :: geo_data(*)
    real(c_float),  intent(in)  :: nwp_data(*)
    real(c_float),  intent(in)  :: ancil_data(*)
    integer(c_int), intent(out) :: clm_result(*)
    integer(c_int), intent(out) :: qa_result(*)
    ! ... 调用原有 Fortran 子程序
  end subroutine

end module cloudmask_c_api
```

**注意**：Fortran 数组默认列优先（column-major），C++/Python 默认行优先（row-major）。在 C++ wrapper 中统一处理转置，不要在 Fortran 端改数组布局。

### 4.2 RecaliManager（再定标核心模块）

```
src/cpp/
  recali/
    RecaliManager.hpp
    RecaliManager.cpp
    RecaliCoeff.hpp     ← 系数数据结构
    RecaliLoader.cpp    ← 从文件加载系数
```

`RecaliManager.hpp` 接口设计：

```cpp
#pragma once
#include <chrono>
#include <filesystem>
#include <map>
#include <optional>
#include <vector>

// 7 个太阳反射波段的再定标系数
struct RecaliCoeff {
    std::array<double, 7> gain;    // 增益系数
    std::array<double, 7> offset;  // 偏置系数
    std::chrono::year_month_day valid_date;
    std::filesystem::path source_file;
};

class RecaliManager {
public:
    // 扫描 recali_root 下所有 yyyymmdd 目录，建立索引
    explicit RecaliManager(const std::filesystem::path& recali_root);

    // 精确匹配指定日期的系数；无则返回 nullopt
    std::optional<RecaliCoeff> get_exact(
        const std::chrono::year_month_day& date) const;

    // 返回最近可用日期的系数（向前查找，最多回退 n_days 天）
    std::optional<RecaliCoeff> get_nearest(
        const std::chrono::year_month_day& date,
        int n_days = 30) const;

    // 列出所有已索引的日期
    std::vector<std::chrono::year_month_day> available_dates() const;

    // 诊断：列出指定日期范围内缺失再定标系数的日期
    std::vector<std::chrono::year_month_day> missing_dates(
        const std::chrono::year_month_day& start,
        const std::chrono::year_month_day& end) const;

private:
    std::map<std::chrono::year_month_day, RecaliCoeff> index_;
    void scan_directory(const std::filesystem::path& root);
    RecaliCoeff load_coeff_file(const std::filesystem::path& path) const;
};
```

### 4.3 IOManager（HDF5 读写）

替换 `io_module.f90`，使用 HDF5 C++ API：

```
src/cpp/
  io/
    IOManager.hpp
    L1Reader.cpp        ← 读 FY-3D L1 HDF5（25 通道）
    L2Writer.cpp        ← 写 CLM/CLA/SST L2 HDF5
    GribReader.cpp      ← GRIB2 → 内存数组（调用 eccodes）
```

关键设计决策：

- L1 数据读入后统一存为 `std::vector<float>`，行优先，与 Python numpy 兼容
- L2 写出时按 FY-3D 标准格式打包 6-byte bit array + QA 字段
- GRIB2 读取优先用 eccodes C API，备选方案是继续调用 wgrib2 二进制并解析输出

### 4.4 NWPInterp（NWP 时空插值）

替换 `read_nwp_data_module.f90` 和 `nwp_utils_module.f90` 的非计算部分：

```cpp
// 时间插值：在两个预报时次之间线性插值到观测时刻
// 空间插值：从 NWP 格点双线性插值到卫星像元位置
// 垂直插值：从 NWP 层次插值到 101 标准压力层
class NWPInterp {
public:
    void load_forecast(const std::filesystem::path& t0_file,
                       const std::filesystem::path& t1_file,
                       double obs_time_fraction);  // 0.0~1.0

    // 批量空间插值，返回 [n_pixels × n_levels] 的数组
    std::vector<float> interp_to_pixels(
        const std::vector<float>& lat,
        const std::vector<float>& lon,
        const std::string& variable_name) const;
};
```

**注意**：原 Fortran 插值子程序（`numerical.f90`）计算精度已经验证，C++ 层只做文件解析和坐标路由，不重新实现插值数学。

### 4.5 双套系数并行运行

在 C++ 层封装一个 `DualRunContext`，一次调用同时执行业务定标和再定标两套：

```cpp
struct DualRunResult {
    std::vector<int> clm_biz;     // 业务定标结果
    std::vector<int> clm_recali;  // 再定标结果
    std::vector<int> qa_biz;
    std::vector<int> qa_recali;
    std::optional<RecaliCoeff> used_coeff;  // 记录实际使用的系数
};

DualRunResult run_dual(const L1Data& l1,
                       const NWPData& nwp,
                       const AncilData& ancil,
                       const RecaliManager& recali_mgr,
                       const std::chrono::year_month_day& obs_date);
```

两套结果写到同一 `yyyymmdd/` 目录下的不同子路径：

```
/data/Data_yuq/fy3_cloud/
  20220803/
    biz/          ← 业务定标结果
      CLM_*.HDF
    recali/       ← 再定标结果
      CLM_*.HDF
    metadata.json ← 使用的系数文件路径、版本、运行时间
```

### 4.6 pybind11 绑定

```cpp
// src/python_bindings.cpp
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include "recali/RecaliManager.hpp"
#include "io/IOManager.hpp"

namespace py = pybind11;

PYBIND11_MODULE(fylat_py, m) {
    m.doc() = "FYLAT C++ extension for Python";

    py::class_<RecaliManager>(m, "RecaliManager")
        .def(py::init<std::string>())
        .def("available_dates", &RecaliManager::available_dates_str)
        .def("missing_dates",   &RecaliManager::missing_dates_str)
        .def("get_nearest",     &RecaliManager::get_nearest_py);

    py::class_<DualRunResult>(m, "DualRunResult")
        .def_readonly("clm_biz",    &DualRunResult::clm_biz)
        .def_readonly("clm_recali", &DualRunResult::clm_recali)
        .def_readonly("used_coeff", &DualRunResult::used_coeff_info);

    m.def("run_dual", &run_dual_py,
          py::arg("l1_path"), py::arg("nwp_paths"),
          py::arg("ancil_cfg"), py::arg("recali_root"), py::arg("date"));
}
```

### 4.7 阶段二验收标准

- C++ 单元测试覆盖：`RecaliManager` 的索引、查找、缺失日期报告
- 对 `20220803` 的输出，C++/Fortran 混合链路与纯 Fortran 基线数值一致（CLM bit 级别）
- 双套系数结果正确写入两个子目录
- `pybind11` 模块可在 conda 环境中 `import fylat_py`

---

## 5. 阶段三：Python 调度层升级

**目标**：将现有 `paral_bat_driver_mersi_ii_fylat.py` 重构为结构化的调度框架，消除手写日期字符串和 subprocess 调用。

### 5.1 目录结构

```
python/fylat/
  __init__.py
  cli.py              ← 命令行入口（argparse / click）
  batch/
    scheduler.py      ← BatchScheduler
    task.py           ← FylatTask 数据类
  recali/
    manager.py        ← Python 侧的 RecaliManager 封装
    comparer.py       ← 双套结果对比统计
  config/
    builder.py        ← 替代手写 .nml
    schema.py         ← dataclass 配置结构
  viz/
    cloudmask_plot.py ← 云掩膜可视化
    compare_plot.py   ← 与 MYD35 对比图
  tests/
    test_scheduler.py
    test_recali.py
    test_comparer.py
```

### 5.2 配置结构（替代 .nml）

```python
# python/fylat/config/schema.py
from dataclasses import dataclass, field
from datetime import date
from pathlib import Path
from typing import Optional

@dataclass
class FylatConfig:
    # 传感器
    sensor_id: int                        # 21=FY-3D, 22=FY-3E
    nwp_opt: int                          # 1-10，选择 NWP 数据源

    # 算法开关
    enable_cloudmask: bool = True
    enable_cloudamount: bool = False
    enable_sst: bool = False

    # 路径
    l1_root: Path = Path("/data/Data_yuq/mersi")
    nwp_root: Path = Path("/data/nwp")
    myd35_root: Path = Path("/data/Data_yuq/aqua_modis/MYD35_L2")
    recali_root: Path = Path("/data/Data_yuq/fy3d_recali")
    output_root: Path = Path("/data/Data_yuq/fy3_cloud")

    # 再定标
    use_recali: bool = True
    recali_fallback_days: int = 30        # 找不到精确日期时最多回退天数

    @property
    def output_dir(self) -> Path:
        return self.output_root           # 具体日期子目录由 Scheduler 管理

    def to_dict(self) -> dict:
        """序列化为 JSON/TOML，随结果文件一起存档"""
        ...
```

### 5.3 BatchScheduler

```python
# python/fylat/batch/scheduler.py
from concurrent.futures import ProcessPoolExecutor
from datetime import date, timedelta
from pathlib import Path
from typing import Iterator
import logging

import fylat_py  # pybind11 扩展

from .task import FylatTask, TaskResult
from ..config.schema import FylatConfig
from ..recali.manager import PythonRecaliManager

class BatchScheduler:
    def __init__(self, config: FylatConfig, max_workers: int = 4):
        self.config = config
        self.max_workers = max_workers
        self.recali_mgr = PythonRecaliManager(config.recali_root)
        self._log = logging.getLogger(__name__)

    def date_range(self, start: date, end: date) -> Iterator[date]:
        """生成日期序列，自动跳过无 L1 数据的日期"""
        current = start
        while current <= end:
            l1_files = list(self._find_l1_files(current))
            if not l1_files:
                self._log.warning("No L1 data for %s, skipping", current)
            else:
                yield current
            current += timedelta(days=1)

    def missing_recali_dates(self, start: date, end: date) -> list[date]:
        """提前报告哪些日期没有再定标系数（而非运行时静默降级）"""
        return self.recali_mgr.missing_dates(start, end)

    def run_date(self, obs_date: date) -> TaskResult:
        """单日期处理，双套系数并行，结果写入 yyyymmdd/ 子目录"""
        out_dir = self.config.output_root / obs_date.strftime("%Y%m%d")
        out_dir.mkdir(parents=True, exist_ok=True)

        result = fylat_py.run_dual(
            l1_path=str(self._find_l1_files(obs_date)[0]),
            nwp_paths=self._resolve_nwp_paths(obs_date),
            ancil_cfg=self.config.to_dict(),
            recali_root=str(self.config.recali_root),
            date=obs_date.strftime("%Y%m%d"),
        )

        # 写出结果 + 元数据
        self._write_result(result, out_dir, obs_date)
        return TaskResult(date=obs_date, success=True, output_dir=out_dir)

    def run_batch(self, start: date, end: date) -> list[TaskResult]:
        """多进程批量运行"""
        # 提前报告缺失再定标系数
        missing = self.missing_recali_dates(start, end)
        if missing:
            self._log.warning(
                "%d dates missing recali coeffs (will fall back to biz): %s",
                len(missing), [d.strftime("%Y%m%d") for d in missing]
            )

        dates = list(self.date_range(start, end))
        with ProcessPoolExecutor(max_workers=self.max_workers) as pool:
            futures = {pool.submit(self.run_date, d): d for d in dates}
            results = []
            for fut in futures:
                try:
                    results.append(fut.result())
                except Exception as e:
                    self._log.error("Failed on %s: %s", futures[fut], e)
                    results.append(TaskResult(date=futures[fut],
                                              success=False, error=str(e)))
        return results
```

### 5.4 RecaliComparer（双套结果对比）

```python
# python/fylat/recali/comparer.py
import numpy as np
import h5py
from pathlib import Path
from datetime import date

class RecaliComparer:
    """对比同一日期业务定标与再定标的 CLM 结果"""

    def compare_date(self, obs_date: date, output_root: Path) -> dict:
        date_str = obs_date.strftime("%Y%m%d")
        biz_clm   = self._load_clm(output_root / date_str / "biz")
        recali_clm = self._load_clm(output_root / date_str / "recali")

        diff = recali_clm.astype(int) - biz_clm.astype(int)
        n_total = diff.size
        n_changed = np.count_nonzero(diff)

        return {
            "date": date_str,
            "n_pixels": n_total,
            "n_changed": n_changed,
            "change_pct": 100.0 * n_changed / n_total,
            # 细分：哪些像元从 clear 变成了 cloud，反之亦然
            "clear_to_cloud": int(np.sum((biz_clm >= 2) & (recali_clm < 2))),
            "cloud_to_clear": int(np.sum((biz_clm < 2) & (recali_clm >= 2))),
        }

    def compare_vs_myd35(self, obs_date: date,
                          output_root: Path,
                          myd35_root: Path) -> dict:
        """与 MYD35 验证集对比，分别计算两套系数的 POD/FAR/HSS"""
        ...
```

### 5.5 CLI 入口

```python
# python/fylat/cli.py
import click
from datetime import date

@click.group()
def cli():
    """FYLAT — FY-3D Cloud Mask Retrieval System"""

@cli.command()
@click.argument("date_str")                          # yyyymmdd
@click.option("--config", default="config.toml")
def run(date_str, config):
    """处理单个日期"""
    ...

@cli.command()
@click.option("--start", required=True)
@click.option("--end",   required=True)
@click.option("--workers", default=4)
def batch(start, end, workers):
    """批量处理日期范围，自动管理再定标系数"""
    ...

@cli.command()
@click.option("--start", required=True)
@click.option("--end",   required=True)
def compare(start, end):
    """对比业务定标与再定标结果，输出统计报告"""
    ...

@cli.command()
@click.option("--start", required=True)
@click.option("--end",   required=True)
def check_recali(start, end):
    """检查日期范围内缺失再定标系数的日期"""
    ...

if __name__ == "__main__":
    cli()
```

### 5.6 阶段三验收标准

- `fylat batch --start 20220801 --end 20220810` 正常运行，输出结构化目录
- 缺失再定标系数时打印警告而非静默降级，`check_recali` 命令可提前核查
- `fylat compare` 输出 CSV/JSON 格式的对比统计
- 对 `20220803` 的结果与 MYD35 对比指标与原 Fortran 版本一致

---

## 6. 阶段四：验证与收尾

### 6.1 数值回归测试框架

建立基于 `20220803` 的黄金输出文件，作为持续集成的回归基准：

```python
# tests/regression/test_baseline.py
def test_20220803_clm_matches_baseline():
    """核心回归测试：确保重构后 CLM 输出与 Fortran 基线一致"""
    result_path = Path("tests/fixtures/20220803/biz/CLM_expected.HDF")
    actual_path = run_single_date("20220803", use_recali=False)

    expected = load_clm(result_path)
    actual   = load_clm(actual_path)

    # CLM 4 级置信度的分布应完全一致
    assert np.array_equal(expected, actual), \
        f"CLM mismatch: {np.sum(expected != actual)} pixels differ"
```

### 6.2 性能基准

记录各阶段改造后的处理时间，确保重构没有带来性能退化：

| 指标 | 基线（纯 Fortran） | 目标（C++/Fortran 混合） |
|------|-------------------|------------------------|
| 单景处理时间 | TBD（运行基线后填写） | ≤ 基线 × 1.1 |
| 批量 10 天处理时间 | TBD | ≤ 基线 × 0.8（多进程） |
| 内存峰值 | TBD | ≤ 基线 × 1.2 |

### 6.3 文档更新

- `README.md`：更新构建步骤（CMake 替代 make）
- `CLAUDE.md`：更新架构说明和测试日期约定
- `docs/recali_guide.md`：新增再定标系数目录规范文档
- `docs/api.md`：C++ API 和 Python API 文档

### 6.4 阶段四验收标准

- CI（GitHub Actions）自动运行回归测试
- 所有阶段的 Git tag 清晰：`v0.1-fortran-baseline`、`v0.2-cmake`、`v0.3-cpp-glue`、`v1.0-full-refactor`
- `PLAN.md` 中所有验收标准全部打勾

---

## 7. 技术决策说明

### 7.1 为什么不重写 Fortran 决策树

Fortran 的云检测决策树（约 20 个场景子程序 + ~100 个 `*_thr.inc` 阈值文件）是经过 MODIS MOD35 框架长期验证的算法实现。重写的风险远大于收益：

- 数值精度：Fortran 的浮点行为与 C++ 存在细微差异，可能引入系统性偏差
- 维护成本：算法本身还在持续调整（阈值、场景分支），双份代码同步维护成本高
- 性能收益：云掩膜计算不是瓶颈，I/O 和 NWP 插值才是

**结论**：Fortran 计算核心原地保留，通过 `ISO_C_BINDING` 暴露接口，不做逻辑重写。

### 7.2 为什么选 pybind11 而不是 f2py

| 方案 | 优点 | 缺点 |
|------|------|------|
| f2py | 直接包装 Fortran | 对 Fortran 派生类型支持差；数组布局隐式转换易出错 |
| pybind11 + ISO_C_BINDING | 接口清晰可控；C++ 中间层可以做内存布局转换 | 多一层封装 |
| ctypes | 零依赖 | 类型安全差；复杂数据结构难以传递 |

**结论**：pybind11 + ISO_C_BINDING 的两层结构，C++ 负责内存布局转换，Python 看到的是干净的 numpy 数组接口。

### 7.3 再定标系数格式约定

当前 `../fy3d_recali/<yyyymmdd>/` 目录下的文件格式需要文档化。建议统一为：

```
fy3d_recali/
  20220803/
    recali_coeff_ch01-07.txt   ← 7 通道 gain/offset，文本格式
    README.txt                 ← 系数来源、版本说明
```

`RecaliManager` 按此格式解析，如有历史遗留文件格式，在 loader 中做兼容处理。

### 7.4 GRIB2 解析方案

| 方案 | 说明 |
|------|------|
| 继续调用 wgrib2 二进制 | 零迁移成本，subprocess 开销可接受 |
| eccodes C API | 原生 C/C++ 集成，无进程创建开销 |

**阶段二推荐**：先用 subprocess 包装 wgrib2（等价于现有 Fortran 做法），后续如有性能需求再迁移到 eccodes。

---

## 8. 里程碑与时间线

```
Week 1-2   ████████░░░░░░░░░░░░░░░░░░  阶段一：基础工程化
           打 tag / 修 clean.sh / CMake / 路径提取 / debug 清理

Week 3-5   ░░░░░░░░████████████░░░░░░  阶段二：C++ 胶水层
           RecaliManager / IOManager / NWPInterp / pybind11 绑定

Week 6-8   ░░░░░░░░░░░░░░░░░░████████  阶段三：Python 调度层
           BatchScheduler / RecaliComparer / CLI / 配置重构

Week 9     ░░░░░░░░░░░░░░░░░░░░░░████  阶段四：验证与收尾
           回归测试 / 性能基准 / 文档更新 / v1.0 tag
```

### 关键里程碑

| 时间 | 里程碑 | 验证方式 |
|------|--------|---------|
| Week 2 末 | v0.2-cmake：CMake 构建成功，基线输出不变 | MD5 对比 |
| Week 5 末 | v0.3-cpp-glue：C++ 层对接 Fortran，双套系数运行 | 数值回归测试 |
| Week 8 末 | v0.4-python：`fylat batch` 命令可用 | 集成测试 |
| Week 9 末 | v1.0：全部验收标准通过 | CI 绿灯 |

---

## 9. 风险与缓解措施

| 风险 | 可能性 | 影响 | 缓解措施 |
|------|--------|------|---------|
| Fortran/C++ 内存布局转换引入数值误差 | 中 | 高 | 在 C++ wrapper 中明确指定列/行优先；回归测试 bit 级别对比 |
| HDF5 C++ API 与现有 Fortran HDF5 输出格式不兼容 | 低 | 高 | 阶段二先读后写，用 h5diff 对比输出 |
| 再定标系数目录格式不统一（历史文件） | 高 | 中 | RecaliLoader 做多格式兼容；扫描时记录无法解析的文件 |
| CMake 在当前 HPC 环境（Intel ifort + 老版 CMake）下兼容性问题 | 中 | 中 | 锁定 CMake 最低版本 3.20；保留原 makefile 作为备份 |
| pybind11 与 conda Python 版本不匹配 | 低 | 中 | 在 conda 环境中同步安装 pybind11；CI 固定 Python 版本 |
| 重构期间算法本身需要修改（阈值调整） | 中 | 低 | Fortran 核心不动，阈值修改只在 `*_thr.inc` 文件，与重构正交 |

---

## 附录 A：文件改造对应关系

| 现有 Fortran 文件 | 改造后归属 | 操作 |
|-----------------|-----------|------|
| `io_module.f90` | C++ `IOManager` | 重写 |
| `read_nwp_data_module.f90` | C++ `NWPInterp` | 重写 |
| `nwp_utils_module.f90` | C++ `NWPInterp` | 部分重写（计算保留 Fortran） |
| `get_ancil_data_module.f90` | C++ `AncilLoader` | 重写 |
| `fylat_FY3_MERSI_II_PGS_Driver.f90` | C++ `FylatPipeline` + Python CLI | 重写 |
| `cloudmask/*.f90` | Fortran 保留 | 添加 `ISO_C_BINDING` wrapper |
| `rtm_tran_module.f90` | Fortran 保留 | 添加 `ISO_C_BINDING` wrapper |
| `planck_module.f90` | Fortran 保留 | 不动 |
| `numerical.f90` | Fortran 保留 | 不动 |
| `constant.f90` | Fortran 保留 | 不动 |
| `platform_module.f90` | C++ `SensorConfig` | 重写（消除 magic number） |
| `data_arrays_module.f90` | C++ struct | 重写 |
| `paral_bat_driver_mersi_ii_fylat.py` | Python `BatchScheduler` | 重构 |

---

## 附录 B：依赖清单

```toml
# pyproject.toml（Python 依赖）
[project]
dependencies = [
  "h5py>=3.8",
  "numpy>=1.24",
  "click>=8.0",
  "tomllib",           # Python 3.11+ 内置，3.10 以下用 tomli
  "matplotlib>=3.7",   # 可视化
  "dacite>=1.8",       # dataclass 反序列化
]

[project.optional-dependencies]
dev = ["pytest", "pytest-cov", "mypy", "ruff"]
```

```cmake
# CMake 依赖（系统级）
# - Intel oneAPI Fortran (ifort) >= 2021
# - HDF5 >= 1.12（含 Fortran 和 C++ binding）
# - pybind11 >= 2.11
# - eccodes >= 2.25（可选，用于 GRIB2）
# - CMake >= 3.20
```

---

*本文档随改造进展持续更新。各阶段完成后在对应验收标准处打 ✅。*