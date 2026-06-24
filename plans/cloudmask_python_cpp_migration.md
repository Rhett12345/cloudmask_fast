# FY-3D MERSI-II Cloud Mask: Fortran → Python + C++ 混编迁移方案

**编写日期**: 2026-06-24
**代码规模**: ~37,200 行 Fortran 77/90 + C 混合代码
**目标**: 将 FYLAT V3.1 逐步迁移为 **Python 编排 + C++ 高性能计算内核** 的混编架构
**约束**: 性能最大化，精度零损失

---

## 一、架构总览

```
┌─────────────────────────────────────────────────────────┐
│                    Python 层 (编排 + 可视化)               │
│  ┌──────────┐  ┌──────────┐  ┌──────────┐  ┌─────────┐ │
│  │ YAML配置  │  │ 工作流调度 │  │ 可视化    │  │ 批处理   │ │
│  └────┬─────┘  └────┬─────┘  └──────────┘  └─────────┘ │
│       │              │                                   │
│  ┌────▼──────────────▼────────────────────────────────┐ │
│  │           pybind11 零拷贝桥接层                       │
│  │     numpy.ndarray ↔ C++ 指针 (共享内存)              │ │
│  └────────────────────┬───────────────────────────────┘ │
└───────────────────────┼─────────────────────────────────┘
                        │
┌───────────────────────┼─────────────────────────────────┐
│                   C++ 层 (计算内核)                       │
│  ┌──────────┐  ┌──────▼──────┐  ┌───────────────────┐  │
│  │ HDF5 I/O │  │ 像素循环引擎  │  │ NWP/RTM 计算      │  │
│  │ (hdf5-cpp)│  │ (OpenMP并行) │  │ (PFAAST/Planck)  │  │
│  └──────────┘  └─────────────┘  └───────────────────┘  │
│  ┌──────────┐  ┌─────────────┐  ┌───────────────────┐  │
│  │ 阈值管理  │  │ 空间一致性检验│  │ 中值滤波 (SIMD)   │  │
│  │ (编译时常量)│  │ (3×3窗口)   │  │ (optmed优化)     │  │
│  └──────────┘  └─────────────┘  └───────────────────┘  │
└─────────────────────────────────────────────────────────┘
```

### 设计原则

| 原则 | 说明 |
|------|------|
| **精度零损失** | C++ 使用与 Fortran 完全相同的 float/double 精度，阈值常量原样移植，不做任何近似 |
| **性能最大化** | 计算密集的像素循环全部在 C++ 中用 OpenMP 并行，Python 只做调度 |
| **零拷贝** | pybind11 让 numpy 数组直接映射到 C++ 内存指针，无数据复制开销 |
| **逐步迁移** | 每迁移一个模块，通过 pybind11 暴露给 Python，系统始终可运行 |
| **Fortran 共存** | 迁移过程中未迁移的 Fortran 模块通过 bind(C) + extern "C" 继续工作 |

---

## 二、性能瓶颈分析

### 关键热路径 (占总运行时间 >90%)

```
fy3mersi_cloud_mask()           ← 主循环：2048×10000+ 像素
  ├── pxinit()                  ← 像素初始化
  ├── classify_surface()        ← 表面类型分类
  ├── [场景子程序]              ← LandDay/ocean_day/PolarDay_* 等 ~20 个
  │   ├── conf_test()           ← 置信度计算 (每像素调用 4-8 次)
  │   ├── tview()               ← 角度插值查表
  │   ├── trispc()              ← 三光谱回归
  │   └── set_bit/clear_bit()   ← 位操作
  ├── spatial_var()             ← 3×3 空间方差
  ├── chk_spatial2()            ← 空间一致性检验
  ├── shadows()                 ← 阴影检测
  ├── thin_ci_chk_ir()          ← 卷云检测
  └── fill_bit_pixel()          ← 最终置信度填充
```

### 迁移优先级矩阵

| 优先级 | 模块 | 原因 |
|--------|------|------|
| P0 | 像素主循环 + 场景子程序 | 占 >80% 运行时间 |
| P0 | conf_test / tview / trispc | 每像素调用多次，热路径 |
| P0 | spatial_var / chk_spatial2 | 3×3 窗口操作，计算密集 |
| P0 | optmed 中值滤波 | 已有 C 优化版本，直接复用 |
| P1 | Planck 函数 | 辐射传输计算 |
| P1 | RTM (PFAAST) | 大气透过率计算 |
| P2 | NWP 数据读取与插值 | I/O 密集 + 空间插值 |
| P2 | HDF5 产品写入 | I/O 密集 |
| P3 | 配置读取 / 辅助数据加载 | 一次性操作，非瓶颈 |

---

## 三、目录结构

```
fy3_cloudmask/
├── CMakeLists.txt                    # 顶层构建
├── pyproject.toml                    # Python 包管理
├── setup.py                          # setuptools + pybind11
│
├── cpp/                              # C++ 计算内核
│   ├── CMakeLists.txt
│   ├── include/
│   │   └── fylat/
│   │       ├── types.hpp             # 数据类型
│   │       ├── constants.hpp         # 物理常量
│   │       ├── thresholds.hpp        # 云检测阈值
│   │       ├── config.hpp            # 配置结构体
│   │       └── arrays.hpp            # 全局数组容器
│   ├── src/
│   │   ├── core/
│   │   │   ├── numerical.cpp         # 数值工具
│   │   │   ├── planck.cpp            # Planck 函数
│   │   │   └── bitops.cpp            # 位操作
│   │   ├── io/
│   │   │   ├── hdf5_reader.cpp       # HDF5 读取
│   │   │   ├── hdf5_writer.cpp       # HDF5 写入
│   │   │   ├── nwp_reader.cpp        # NWP GRIB 读取
│   │   │   └── ancillary_reader.cpp  # 辅助数据
│   │   ├── rtm/
│   │   │   ├── rtm_tran.cpp          # 大气透过率
│   │   │   └── rtm_utils.cpp         # RTM 工具
│   │   ├── cloudmask/
│   │   │   ├── cloud_mask_engine.cpp # 主像素循环引擎
│   │   │   ├── scene_tests.cpp       # 场景测试
│   │   │   ├── confidence.cpp        # 置信度计算
│   │   │   ├── spatial.cpp           # 空间检验
│   │   │   ├── spectral.cpp          # 光谱检验
│   │   │   ├── surface_classify.cpp  # 表面分类
│   │   │   └── auxiliary_checks.cpp  # 辅助检验
│   │   └── nwp/
│   │       ├── nwp_interp.cpp        # NWP 插值
│   │       └── nwp_utils.cpp         # NWP 工具
│   └── bindings/
│       └── pybind_module.cpp         # pybind11 绑定
│
├── python/                           # Python 层
│   ├── fylat/
│   │   ├── __init__.py
│   │   ├── config.py                 # YAML 配置解析
│   │   ├── pipeline.py               # 主处理流水线
│   │   ├── cloudmask.py              # 云检测封装
│   │   ├── io/
│   │   │   ├── l1b_loader.py
│   │   │   └── product_writer.py
│   │   └── visualize/
│   │       ├── plot_cloudmask.py
│   │       ├── plot_diagnostics.py
│   │       └── compare_results.py
│   └── scripts/
│       ├── run_single.py
│       └── run_batch.py
│
├── config/                           # 配置文件
│   ├── default.yaml
│   ├── sensor_mersi2.yaml
│   └── nwp_gfs0p25.yaml
│
├── coeff/                            # 系数文件 (保持不变)
│
├── fortran/                          # 原始 Fortran 代码
│
├── tests/
│   ├── cpp/                          # C++ 单元测试
│   ├── python/                       # Python 测试
│   │   └── test_precision.py         # 精度对比测试
│   └── data/                         # 测试数据
│       └── expected_clm.bin          # Fortran 基准输出
│
└── scripts/
    ├── build.sh
    ├── run_precision_check.sh
    └── benchmark.sh
```

---

## 四、分阶段实施计划

### 阶段 0：基础设施 (1-2 周)

**目标**: 搭建 C++ + pybind11 构建系统，迁移常量和配置。

#### 0.1 构建系统
- CMake 顶层项目，同时编译 C++ 和 Fortran
- pybind11 绑定模块，可 import fylat 在 Python 中
- Google Test 测试框架

#### 0.2 常量与类型移植

| Fortran 源文件 | C++ 目标 | 内容 |
|---------------|---------|------|
| constant.f90 | constants.hpp | 物理常量 (PI, C_1, C_2 等) |
| names_module.f90 | config.hpp | 配置结构体 |
| data_arrays_module.f90 | arrays.hpp | 全局数据数组 |
| cloudmask_data_arrays.f90 | scene_context.hpp | 像素级上下文 |
| global.inc | constants.hpp | inband, bad_data 等 |
| symbol_struct | enum class | C++ 强类型枚举 |

#### 0.3 配置系统 (YAML 替代 namelist)

```yaml
sensor:
  id: 21
  nwp_opt: 8
  rtm_opt: 1
input:
  l1b_data: ""
  geo_data: ""
output:
  cloud_mask: ""
algorithms:
  cloudmask: true
```

#### 0.4 验证基线
用 Fortran 处理 3 个测试场景，保存输出作为精度基线。

---

### 阶段 1：计算内核移植 (3-4 周)

**目标**: 将像素级计算热路径全部移植到 C++。

#### 1.1 纯计算函数

| Fortran | C++ | 说明 |
|---------|-----|------|
| conf_test.f | confidence::conf_test() | 置信度计算 |
| conf_test_2val.f | confidence::conf_test_2val() | 双阈值置信度 |
| tview.f | spectral::tview() | 角度插值查表 |
| trispc.f | spectral::trispc() | 三光谱回归 |
| check_bits/set_bit/clear_bit | bitops:: | 位操作 |
| optmed.c / optmed_int1.c | median:: | 中值滤波 (直接复用C) |
| get_sg_thresholds.f90 | thresholds::get_sg() | 太阳耀斑阈值 |
| get_pn_thresholds.f | thresholds::get_pn() | 极地夜间阈值 |
| get_nl_thresholds.f90 | thresholds::get_nl() | 非极地夜间阈值 |

#### 1.2 阈值常量 (include -> constexpr)

```cpp
namespace fylat::thresholds::land_day {
    constexpr float dl11_12hi[4] = {-0.0075f, 2.5f, -0.25f, 1.0f};
    constexpr float dlref1[4] = {0.04f, 0.20f, 0.08f, 1.0f};
}
```

#### 1.3 场景测试子程序 (20个)

| Fortran | C++ | 场景 |
|---------|-----|------|
| LandDay.f90 | scene::land_day() | 白天陆地 |
| LandNite.f90 | scene::land_nite() | 夜间陆地 |
| LandDay_desert.f90 | scene::land_day_desert() | 白天沙漠 |
| LandDay_coast.f90 | scene::land_day_coast() | 白天海岸 |
| ocean_day.f90 | scene::ocean_day() | 白天海洋 |
| ocean_nite.f90 | scene::ocean_nite() | 夜间海洋 |
| PolarDay_land.f90 | scene::polar_day_land() | 极地白天陆地 |
| PolarDay_ocean.f90 | scene::polar_day_ocean() | 极地白天海洋 |
| PolarDay_snow.f90 | scene::polar_day_snow() | 极地白天雪 |
| PolarDay_desert.f90 | scene::polar_day_desert() | 极地白天沙漠 |
| PolarDay_coast.f90 | scene::polar_day_coast() | 极地白天海岸 |
| PolarNite_land.f90 | scene::polar_nite_land() | 极地夜间陆地 |
| PolarNite_ocean.f90 | scene::polar_nite_ocean() | 极地夜间海洋 |
| PolarNite_snow.f90 | scene::polar_nite_snow() | 极地夜间雪 |
| Antarctic_day.f90 | scene::antarctic_day() | 南极白天 |
| Day_snow.f90 | scene::day_snow() | 白天雪 |
| Nite_snow.f90 | scene::nite_snow() | 夜间雪 |

#### 1.4 pybind11 绑定

```cpp
void run_cloud_mask(
    py::array_t<float> radiance,
    py::array_t<float> geo_data,
    py::array_t<float> nwp_data,
    py::array_t<int8_t> ancillary,
    py::array_t<float> btclr,
    py::array_t<int8_t> cm_bitarray,
    py::array_t<int8_t> cm_qa_bitarray,
    int sensor_id, int nwp_opt
);

PYBIND11_MODULE(fylat_core, m) {
    m.def("run_cloud_mask", &run_cloud_mask);
}
```

#### 1.5 C++ 云检测引擎 (OpenMP 并行)

```cpp
void CloudMaskEngine::run(...) {
    #pragma omp parallel for schedule(dynamic, 16) collapse(2)
    for (int iline = 0; iline < nLine; ++iline) {
        for (int ielem = 0; ielem < nElem; ++ielem) {
            PixelContext ctx;
            pxinit(ctx);
            extract_pixel(radiance, iline, ielem, ctx.pxldat);
            classify_surface(geo, ancillary, iline, ielem, ctx);
            // 根据 surface type 选择场景测试
            // 空间一致性检验
            // 辅助检验
            // 写入输出
        }
    }
}
```

---

### 阶段 2：I/O 层移植 (2-3 周)

| Fortran | C++ | 数据 |
|---------|-----|------|
| io_module.f90 (L1b) | io::read_l1b() | 25通道辐射 |
| io_module.f90 (GEO) | io::read_geo() | 经纬度角度 |
| read_nwp_data_module.f90 | nwp::read_nwp_data() | NWP数据 |
| nwp_utils_module.f90 | nwp::interp() | 空间/时间插值 |
| get_ancil_data_module.f90 | io::read_ancillary() | 辅助数据 |
| io_module.f90 (write) | io::write_cloud_mask() | 产品写入 |

---

### 阶段 3：辐射传输移植 (1-2 周)

| Fortran | C++ |
|---------|-----|
| planck_module.f90 | rtm::planck_rad2tbb() / tbb2rad() |
| rtm_tran_module.f90 | rtm::tran_vmodis_d101() |
| rtm_utils_module.f90 | rtm::ir_rtm_driver() |
| atmtran.f | rtm::atmtran() |
| frontend_module.f90 | rtm::compute_angles() |

---

### 阶段 4：Python 编排层 (2-3 周)

```python
class FylatPipeline:
    def run(self):
        l1b, geo = read_l1b_geo(config)        # C++ HDF5
        nwp = read_and_interp_nwp(config, geo)  # C++
        ancil = read_ancillary(config)           # C++
        btclr = compute_clear_sky_bt(l1b, nwp)  # C++ RTM
        run_cloud_mask(l1b, geo, nwp, ancil, btclr, cm_out, qa_out)  # C++ 核心
        write_product(cm_out, config.output)     # C++ HDF5
        plot_cloudmask(cm_out, geo)              # Python matplotlib
```

---

### 阶段 5：驱动与批处理 (1 周)

```bash
python -m fylat.run_single config/business.yaml
python -m fylat.run_batch config/batch/ --workers 8
python -m fylat.visualize output/cloudmask.hdf5
```

---

### 阶段 6：精度验证与性能优化 (2-3 周)

#### 精度验证
- Fortran 输出 vs Python+C++ 输出逐比特对比
- 差异率必须 = 0% (比特级一致)
- 允许 <0.01% 浮点舍入差异 (IEEE 754, 1 ULP)

#### 性能目标

| 指标 | Fortran 基线 | Python+C++ 目标 |
|------|-------------|----------------|
| 单景处理时间 | ~X 秒 | <= X 秒 |
| OpenMP 线程 | 1 (串行) | 8+ (并行) |
| 并行加速比 | N/A | >= 6x |

#### 优化手段
1. OpenMP 并行: #pragma omp parallel for schedule(dynamic, 16)
2. SIMD 向量化: conf_test, tview 等用 #pragma omp simd
3. Cache 优化: 按行处理，邻域访问局部性
4. 内存预分配: 循环外分配，避免 malloc
5. 编译优化: -O3 -march=native

---

## 五、关键技术决策

### 5.1 数组维度
C++ 使用 0-based 行优先 (与 numpy 一致)，阈值索引做 -1 偏移。

```fortran
masv66 = pxldat(3)   ! Fortran: 1-based
```
```cpp
float masv66 = pxldat[2];  // C++: 0-based (3-1=2)
```

### 5.2 全局状态
C++ 用 struct PixelContext 封装每像素状态，通过函数参数传递，天然线程安全。

### 5.3 阈值管理
C++ 用 constexpr 数组 + 命名空间，编译时内联。

### 5.4 错误处理
C++ 用异常 (std::runtime_error)，Python 层捕获并记录日志。

---

## 六、风险与缓解

| 风险 | 缓解措施 |
|------|---------|
| 浮点精度差异 | 相同精度，避免 -ffast-math，比特级对比测试 |
| HDF5 API 差异 | 官方 C++ API，逐字段对比 |
| OpenMP 竞态 | PixelContext 纯栈分配，无共享写入 |
| 阈值移植错误 | 自动化脚本从 .inc 生成 .hpp，人工校验 |

---

## 七、时间线

```
周  1-2:  阶段 0 - 基础设施
周  3-6:  阶段 1 - 计算内核 (核心)
周  7-8:  阶段 2 - I/O 层
周  9-10: 阶段 3 - 辐射传输
周 11-13: 阶段 4 - Python 编排
周 14:    阶段 5 - 驱动与批处理
周 15-17: 阶段 6 - 精度验证与优化

总计: ~17 周 (约 4 个月)
```

---

## 八、最终形态

```bash
pip install -e .
python -m fylat.run_single config/business.yaml
python -m fylat.run_batch config/batch/ --workers 8
python -m fylat.visualize output/cloudmask.hdf5
python -m fylat.verify --fortran-baseline baseline/ --cpp-output output/
```

Fortran 代码迁移完成后移除，最终项目仅包含 Python + C++。

---

## 附录：迁移时长修订评估 (Planner Agent 评审)

### 总览

原始 17 周计划对**单人开发**偏乐观约 30-40%，对 **2 人团队**基本准确（需重构以利用更多并行）。

| 场景 | 周数 | 置信度 |
|------|------|--------|
| 原计划 (单人) | 17 | 低 — 可能超期 30-40% |
| **修订 (单人)** | **19** | 中 — 含 2 周缓冲 |
| 最坏 (单人) | 22 | HDF5 + 数值验证出问题时 |
| **修订 (2人)** | **11** | 中 — 并行度好 |
| 最坏 (2人) | 14 | 同风险因素 |

---

### 模块难度重评

#### 比看起来简单的模块 (可更快完成)

| 模块 | 原估 | 修订 | 原因 |
|------|------|------|------|
| F77 工具函数 (15个, ~4125行) | 1.5周 | 3-4天 | 纯计算，无I/O无状态，机械翻译，每天可翻3-5个 |
| 阈值 include 文件 (25个, ~2000行) | 1周 | 1-2天 | 常量定义，find-and-replace |
| 字符串工具 (~200行) | - | 半天 | 直接用 std::string 替代 |
| constant.f90 (~400行) | - | 1天 | 物理常量，机械翻译 |
| names_module.f90 (~100行) | - | 半天 | 变量声明，trivial |
| platform_module.f90 (~200行) | - | 半天 | 平台信息 |
| message_module.f90 (~100行) | - | 1天 | 用 spdlog 替代 |
| cloudmask_data_arrays.f90 (~80行) | - | 1天 | OMP threadprivate → C++ struct |
| optmed.c (314行) | - | 1天 | 已经是C，加 pybind11 包装即可 |

#### 比看起来难的模块 (需要更多时间)

| 模块 | 原估 | 修订 | 原因 |
|------|------|------|------|
| **io_module.f90 (~3500行)** | 1.5周 | **3周** | HDF5 Fortran API ≠ C API，函数名/类型/内存布局全不同，含 compound type |
| **read_nwp + nwp_utils (~2700行)** | 1周 | **2周** | wgrib2 C interop + 空间/时间插值的网格变换 |
| **get_ancil_data (~2000行)** | 1周 | **2周** | 5个独立数据源，5种格式，5种分辨率 |
| **rtm_tran + atmtran (~900行)** | 1周 | **1.5周** | PFAAST 系数表 + EQUIVALENCE/COMMON 块处理 |
| **像素主循环 (~1200行)** | 1周 | **2周** | 20个场景子程序调度 + OpenMP 并行化 |

#### 评估正确的模块

| 模块 | 周数 | 说明 |
|------|------|------|
| 20个场景测试 (~11000行) | 2-3周 | 核心科学，每个需独立验证 |
| planck_module.f90 (~800行) | 1周 | 定义明确的物理 |
| frontend_module.f90 (~400行) | 3-5天 | 角度计算 |
| numerical.f90 (~1200行) | 1周 | 中值滤波已有C实现 |
| data_arrays_module.f90 (~900行) | 3-5天 | 数据结构定义 |

---

### 关键路径

```
周1: 基础设施
  → 周3-6: io_module.f90 (瓶颈!)
    → 周7-8: RTM
      → 周9-13: 场景测试 + 像素循环
        → 周14-15: OpenMP + 集成
          → 周16-19: 验证
```

**io_module.f90 (3500行)** 是单点最大风险：最大文件、最复杂API翻译、下游全部依赖。

---

### 修订周计划 (单人开发)

#### 阶段 0: 基础设施 (第1周)
- CMake + pybind11 构建系统
- 翻译 constant.f90, names_module.f90, data_arrays_module.f90
- 翻译 platform_module.f90, message_module.f90
- 翻译 25个 .inc 阈值文件 → constexpr
- 翻译字符串工具 → std::string
- Google Test + pybind11 骨架

**里程碑**: C++ 库编译通过，pybind11 模块可加载，常量单测通过。

#### 阶段 1: F77 工具 + 数值 (第2-3周)
- 第2周: 翻译 conf_test, tview, trispc, set_bit 等 F77 函数 (batch 1)
- 第3周: 翻译剩余 F77 函数 + numerical.f90 + cloudmask_data_arrays.f90
- 包装 optmed.c

**里程碑**: 所有纯计算函数迁移完成并测试通过。

#### 阶段 2: I/O 层 (第3-6周) ← 关键瓶颈
- 第3周 (与第3周并行): io_module.f90 part 1 — HDF5 L1b 读取
- 第4周: io_module.f90 part 2 — 产品写入 + GEO 读取
- 第5周: NWP 模块 (read_nwp + nwp_utils)
- 第6周: 辅助数据 (get_ancil_data) + hdf4.f90

**里程碑**: 完整 I/O 流水线可用，能读写所有数据格式。

#### 阶段 3: 辐射传输 (第7-8周)
- 第7周: planck_module + frontend_module
- 第8周: rtm_tran_module + rtm_utils_module + atmtran.f

**里程碑**: RTM 输出与 Fortran 参考一致 (<0.01%)。

#### 阶段 4: 核心算法 (第9-13周)
- 第9-10周: 场景测试 batch 1 (10个: LandDay, ocean_day 等)
- 第11-12周: 场景测试 batch 2 (10个: PolarDay_*, spatial tests)
- 第13周: 像素主循环 + 所有场景测试集成

**里程碑**: 端到端流水线可产出单景云检测结果。

#### 阶段 5: OpenMP + 云量 (第14-15周)
- 第14周: OpenMP 并行化 + 性能基准
- 第15周: cloud_amount 算法 + 批处理

**里程碑**: 并行执行可用，批处理功能完整。

#### 阶段 6: 验证 + 打磨 (第16-19周)
- 第16-17周: 逐比特精度验证 + 边缘情况测试
- 第18周: Python CLI + 可视化集成 + YAML 配置
- 第19周: 文档 + 部署

**里程碑**: 生产就绪系统。

---

### 修订周计划 (2人团队)

| 周 | 开发者 1 | 开发者 2 |
|----|----------|----------|
| 1 | 基础设施 + 类型 + 常量 | F77 工具 batch 1 |
| 2 | io_module.f90 part 1 | F77 工具 batch 2 + numerical |
| 3 | io_module.f90 part 2 | cloudmask_data_arrays + 场景测试 batch 1 (5个) |
| 4 | NWP 模块 | 场景测试 batch 2 (5个) |
| 5 | 辅助数据 | 场景测试 batch 3 (5个) |
| 6 | RTM (planck + frontend) | 场景测试 batch 4 (5个) + 空间测试 |
| 7 | RTM (PFAAST + atmtran) | 像素主循环 |
| 8 | OpenMP + cloud amount | 像素主循环集成 |
| 9-10 | 逐比特精度验证 | Python 编排 + 可视化 |
| 11 | 文档 + 部署 | 边缘情况测试 |

---

### 高风险项

1. **HDF5 API 翻译** (io_module.f90): 可能额外 +1-2 周
2. **逐比特数值验证**: IEEE 754 舍入差异，可能 +1 周
3. **OpenMP 线程安全**: 可能发现隐藏的共享状态，可能 +1 周

### 关键建议

1. **io_module.f90 从第2周就开始**，它在关键路径上
2. **pybind11 桥接从第1天增量构建**，不要留到最后
3. **迁移前先生成 Fortran 参考输出**：5-10 个代表性场景
4. **阈值 .inc 文件最先翻译**：trivial 且解锁所有场景测试
5. **第一遍不追求逐比特一致**：先功能正确 (99% 像素一致)，再精修
6. **考虑临时保留 Fortran I/O 层**：用 f2py 调用 Fortran I/O，只迁移计算内核到 C++

---

## 附录二：Agent 辅助开发修订 (2026-06-24)

### 前提

开发环境配备 AI Agent (Claude Code) 全程辅助，负责：
- 批量 Fortran → C++ 代码翻译
- 生成 pybind11 绑定、CMake 构建系统
- 生成单元测试
- 写 Python 编排层和可视化模块

开发者的角色转变为：**审查、验证、提供测试数据、科学正确性判定**。

### Agent 加速倍数分析

| 工作类型 | Agent 擅长程度 | 加速倍数 |
|----------|--------------|---------|
| 机械代码翻译 (常量/阈值/F77工具) | 极高 — 模式明确，批量处理 | 5x |
| API 映射 (HDF5 Fortran→C) | 高 — 有文档，但需验证 | 2.5x |
| 算法翻译 (场景测试/RTM) | 高 — 逻辑明确 | 3x |
| 构建系统 (CMake/pybind11) | 极高 — 标准化 | 5x |
| 测试生成 | 极高 — 自动生成 | 5x |
| Python 编排层 | 极高 — 标准模式 | 5x |
| **精度验证** | **低 — 需要真实数据和人工判读** | **1.5x** |
| **性能调优** | **低 — 需要真实硬件** | **1.2x** |

### Agent 辅助时间线 (单人 + Agent)

#### 第 1 周: 基础设施 + 纯计算模块

**Day 1-2: 项目骨架**
- Agent 生成 CMake 构建系统 + pybind11 模块骨架
- Agent 翻译 constant.f90 → constants.hpp (物理常量、symbol_struct → enum class)
- Agent 翻译 names_module.f90 → config.hpp
- Agent 翻译 data_arrays_module.f90 → arrays.hpp
- Agent 翻译 cloudmask_data_arrays.f90 → scene_context.hpp
- Agent 翻译 platform_module.f90, message_module.f90
- Agent 翻译 25个 .inc 阈值文件 → constexpr headers
- Agent 翻译字符串工具 → std::string
- **开发者**: 审查生成的代码，确认常量值正确

**Day 3-4: F77 工具函数批量翻译**
- Agent 批量翻译 15 个 F77 函数:
  - conf_test.f, conf_test_2val.f
  - tview.f, trispc.f
  - set_bit.f, clear_bit.f, check_bits.f, check_qa_bits.f, set_qa_bit.f
  - get_sg_thresholds.f90, get_pn_thresholds.f, get_nl_thresholds.f90
  - get_regdif.f, get_regstd.f
  - pxinit.f, proc_path.f, set_unused_bits.f, set_confdnc.f, set_quality_A.f
- Agent 为每个函数生成单元测试
- Agent 包装 optmed.c / optmed_int1.c (已是 C 代码)
- **开发者**: 用 Fortran 参考数据验证每个函数的输出

**Day 5: 数值模块 + 编译验证**
- Agent 翻译 numerical.f90 → numerical.cpp
- Agent 翻译 thresholds_read_module.f90
- Agent 翻译 param_read_file.f
- Agent 翻译 fill_bit_pixel.f90
- **开发者**: 全量编译，修复编译错误

**第1周里程碑**: 所有纯计算函数迁移完成，C++ 库编译通过，pybind11 模块可加载。

#### 第 2 周: I/O 层

**Day 1-2: HDF5 读取**
- Agent 翻译 io_module.f90 L1b/GEO 读取部分
  - fylat_read_fy3_mersi_geo_data() → io::read_geo()
  - fylat_read_fy3_mersi_L1b_data() → io::read_l1b()
- Agent 处理 Fortran HDF5 API → C HDF5 API 映射
- **开发者**: 用真实 HDF5 文件验证读取结果

**Day 3: HDF5 写入 + 中间结果**
- Agent 翻译产品写入 (fylat_write_out_cloud_mask 等)
- Agent 翻译 io_module_intermediate.f90
- Agent 翻译 hdf4.f90
- **开发者**: 验证写出的 HDF5 可被其他工具正确读取

**Day 4-5: NWP + 辅助数据**
- Agent 翻译 read_nwp_data_module.f90
- Agent 翻译 nwp_utils_module.f90
- Agent 翻译 get_ancil_data_module.f90 (ecosystem, snow/ice, OISST, emissivity, albedo)
- **开发者**: 用真实 NWP/辅助数据验证插值结果

**第2周里程碑**: 完整 I/O 流水线可用。能读 L1b/GEO/NWP/辅助数据，能写产品。

#### 第 3 周: RTM + 核心算法 batch 1

**Day 1-2: 辐射传输**
- Agent 翻译 planck_module.f90 → rtm::planck_rad2tbb() / tbb2rad()
- Agent 翻译 rtm_tran_module.f90 → rtm::tran_vmodis_d101()
- Agent 翻译 rtm_utils_module.f90 → rtm::ir_rtm_driver()
- Agent 翻译 atmtran.f → rtm::atmtran()
- Agent 翻译 frontend_module.f90 → rtm::compute_angles() / extract_sattime()
- **开发者**: 验证 Planck 函数输出与 Fortran 一致

**Day 3-5: 场景测试 batch 1 (10个)**
- Agent 批量翻译 10 个场景测试子程序:
  - LandDay.f90 → scene::land_day()
  - LandNite.f90 → scene::land_nite()
  - LandDay_desert.f90 → scene::land_day_desert()
  - LandDay_desert_c.f90 → scene::land_day_desert_c()
  - LandDay_coast.f90 → scene::land_day_coast()
  - ocean_day.f90 → scene::ocean_day()
  - ocean_nite.f90 → scene::ocean_nite()
  - Day_snow.f90 → scene::day_snow()
  - Nite_snow.f90 → scene::nite_snow()
  - Antarctic_day.f90 → scene::antarctic_day()
- **开发者**: 逐个用像素级参考数据验证

**第3周里程碑**: RTM + 半数场景测试迁移完成。

#### 第 4 周: 核心算法 batch 2 + 集成

**Day 1-3: 场景测试 batch 2 (10个)**
- Agent 批量翻译:
  - PolarDay_land.f90 → scene::polar_day_land()
  - PolarDay_ocean.f90 → scene::polar_day_ocean()
  - PolarDay_snow.f90 → scene::polar_day_snow()
  - PolarDay_desert.f90 → scene::polar_day_desert()
  - PolarDay_desert_c.f90 → scene::polar_day_desert_c()
  - PolarDay_coast.f90 → scene::polar_day_coast()
  - PolarNite_land.f90 → scene::polar_nite_land()
  - PolarNite_ocean.f90 → scene::polar_nite_ocean()
  - PolarNite_snow.f90 → scene::polar_nite_snow()
  - chk_land.f90, chk_land_nite.f90, chk_coast.f90, chk_sunglint.f90, chk_shallow_water.f
- Agent 翻译空间测试: spatial_var.f, chk_spatial_var.f, chk_spatial2.f
- Agent 翻译辅助检验: shadows.f90, thin_ci_chk_ir.f90, noncld_obs_chk.f90

**Day 4-5: 像素主循环 + OpenMP**
- Agent 翻译 fylat_fy3mersi_cloud_mask.f90 (主循环引擎)
- Agent 集成所有 20 个场景测试
- Agent 添加 OpenMP parallel for
- Agent 翻译 cloudamount/fylat_fy3mersi_cloud_amount.f90
- **开发者**: 编译、链接、端到端测试

**第4周里程碑**: 端到端流水线可产出云检测结果。OpenMP 并行可用。

#### 第 5 周: 精度验证

**Day 1-3: 逐比特对比**
- Agent 编写精度验证脚本 (Python)
- Agent 编写 Fortran vs C++ 输出对比工具
- **开发者**: 运行 3 个基准场景，对比输出
- **开发者**: 修复发现的数值差异
- 目标: 差异率 <0.01%

**Day 4-5: 边缘情况测试**
- Agent 生成边缘情况测试用例
- **开发者**: 测试极地、沙漠、太阳耀斑、云边缘等场景
- **开发者**: 验证 QA 位与 Fortran 一致

**第5周里程碑**: 精度验证通过。

#### 第 6 周: Python 编排 + 可视化 + 打磨

**Day 1-2: Python 编排层**
- Agent 生成 python/fylat/ 包:
  - config.py (YAML 配置解析)
  - pipeline.py (主流水线)
  - cloudmask.py (云检测封装)
  - io/l1b_loader.py, io/product_writer.py
- Agent 生成 scripts/run_single.py, run_batch.py

**Day 3-4: 可视化**
- Agent 生成 visualize/ 模块 (复用现有 visualize/ 代码)
- Agent 生成 compare_results.py (Fortran vs C++ 对比图)

**Day 5: 打磨 + 文档**
- Agent 生成 README、构建文档
- **开发者**: 最终集成测试

**第6周里程碑**: 生产就绪系统。

---

### Agent 辅助时间线总结

| 周 | Agent 主要工作 | 开发者主要工作 |
|----|--------------|--------------|
| 1 | 常量/阈值/F77工具/数值模块翻译 | 审查代码，验证常量，编译调试 |
| 2 | I/O 层翻译 (HDF5/NWP/辅助数据) | 用真实数据验证 I/O |
| 3 | RTM + 场景测试 batch 1 | 验证 RTM + 场景测试 |
| 4 | 场景测试 batch 2 + 像素主循环 + OpenMP | 端到端集成测试 |
| 5 | 精度验证脚本 + 测试用例 | 跑基准场景，修复精度差异 |
| 6 | Python 编排 + 可视化 + 文档 | 最终集成测试 |

**总计: 6 周** (含精度验证)

### 与纯人工对比

| 场景 | 纯人工 | Agent 辅助 | 加速 |
|------|--------|-----------|------|
| 单人开发 | 19 周 | **6 周** | 3.2x |
| 最坏情况 | 22 周 | 8 周 | 2.8x |

### Agent 辅助的局限

Agent 加速的是**代码翻译**，但以下仍需开发者投入时间：

1. **提供测试数据** — Agent 没有你的 HDF5 卫星数据文件
2. **科学正确性验证** — 云检测结果对不对，只有领域专家能判断
3. **编译环境** — 需要在你的服务器上编译、修复环境问题
4. **精度对比** — 需要跑 Fortran 基准输出做对比
5. **性能调优** — 需要在真实硬件上测 OpenMP 效果

**最可能的节奏**: Agent 在几天内完成翻译，但验证和调试需要你持续投入。建议全职投入第 1-5 周，确保验证质量。
