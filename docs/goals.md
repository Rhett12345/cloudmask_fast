# FYLAT 项目迁移与算法升级总体规划

> 文档版本：v1.1（经代码核查修正）  
> 项目仓库：https://github.com/Rhett12345/cloudmask_fast  
> 编写时间：2026-06 | 修正：2026-06-29

---

## 一、背景与问题诊断

### 1.1 现状

FYLAT（FengYun Land-Atmosphere Toolkit）是基于 FY-3D/E MERSI-II 的云检测系统，核心算法框架从 MODIS MOD35 迁移而来，使用 Intel Fortran + C 编写，约 37,000 行。

当前架构：

```
Python（批处理调度）
└── Fortran/C 编译整体（ifort/icc 链接）
     ├── IO 层（HDF5 读写、GRIB2 解析）
     ├── NWP 插值（101 压力层，时空双插值）
     ├── PFAAST RTM（晴空亮温模拟）
     └── 云检测决策树（白天/夜间/极区各场景子模块）
```

### 1.2 两个核心问题

本规划要同时解决两件事，且将其**合并为一个工程过程**——替换哪个模块，就顺手升级该模块的算法，不做重复劳动。

| 问题 | 描述 |
|------|------|
| **工程问题** | Fortran 代码老旧，难以维护、难以集成 ML 模块，逐步全部替换为 C++ |
| **算法问题** | MERSI-II 缺少 MOD35 依赖的关键通道，夜间/极区检测效果差，阈值调整已到上限 |

---

## 二、算法缺陷根因分析

### 2.1 MERSI-II 相对 MODIS 的关键通道差异

MOD35 框架中部分测试依赖 MODIS 特有通道。经代码核查，实际存在 3 类差异：

**真实缺失 / 不匹配（影响大）：**

| MOD35 用途 | MODIS 通道 | 中心波长 | MERSI-II 实际 | 影响 |
|------------|-----------|---------|--------------|------|
| 11-3.9μm BTD 核心测试 | B22 | 3.959 μm | ch20=3.80μm / ch21=4.05μm，BT 偏差 +40~+90K | **夜间所有 BTD 逻辑失效的根因** |
| 极夜晴空恢复 | B33–36 | 13.3–14.2 μm（CO₂） | **无，MERSI 截止于 11.95μm** | 极区夜间 CO₂ 测试全部注释 |
| 高云 / 水汽检测 | B27 | 6.7 μm | ch22=7.233μm，权重函数峰值偏低（中对流层 vs 上对流层） | 6.7μm 测试被注释，高云检测能力减弱 |

**存在等效通道（影响小，无需特殊处理）：**

| MODIS 通道 | MERSI-II 通道 | 波长偏差 | 说明 |
|-----------|--------------|---------|------|
| B28 (7.3μm) | ch22 (7.233μm) | +0.09μm | 功能等效，夜间 7.3-11μm BTD 测试活跃 |
| B29 (8.6μm) | ch23 (8.560μm) | +0.01μm | 同一通道，有专用 RTM 剖面和 8.6-7.3 BTD 测试 |
| B17-18 (0.9μm) | ch16/17/18 (0.91/0.94/0.94μm) | ~0 | 3 个等效波段覆盖

### 2.2 MERSI-II 仪器 SRF 差异导致系统性 BTD 偏差

**SRF 校正系数差异（src/planck_module.f90）——代码确凿证据：**

非单色传感器需要 tci/tcs 系数将 Planck 亮温修正为实际通道亮温。MERSI-II 的系数与 MODIS 差异巨大：

| 通道 | Aqua MODIS tcs | FY-3D tcs | Aqua tci(K) | FY-3D tci(K) | Δtci(K) |
|------|---------------|-----------|-------------|-------------|---------|
| ch24 (10.8μm) | 0.9995 | 0.9980 | 0.1290 | **0.5763** | +0.447 |
| ch25 (12.0μm) | 0.9997 | 0.9984 | 0.0681 | **0.4317** | +0.364 |

BTD = BT(ch24) - BT(ch25)，通道间不对称的 tci/tcs 差异会产生 ~0.1K 量级的系统性 BTD 偏差。这意味着即使 MERSI-II 和 MODIS 观测完全相同的场景，算出的 BTD 也不同——直接沿用 MODIS 的 BTD 阈值在物理上就是不正确的。

注：文档记载 ch20 NEdT=0.25K vs MODIS B22 NEdT=0.07K，但代码中未引用 NEdT 值，该差异未被阈值方案考虑。

### 2.3 现有阈值调整的上限

**阈值沿用 MODIS 原始值——经代码核查确认：**

对比 `coeff/fylat_thresholds.mersi.ii3d.v8`（FY-3D）与 `coeff/fylat_thresholds.mersi.aqua.v8`（Aqua MODIS），夜间 BTD 阈值**字节级完全相同**：

- `nl_11_4l/h/m`（11-4μm BTD）、`nl4_12hi`（4-12μm）、`nl7_11s`（7.2-11μm）、`no11_4lo`（海洋 11-4μm）——全部未改
- 唯一调整的是 `ns4_12hi`（夜间雪 4-12μm，标注"revised by minmin 20190109"）
- `LandNite.f90:310` 有硬编码的 `-1.5K` 补丁修正 3.8μm BT

此外，夜间陆地模块并非文档早期版本所述的"退化为单一 BTD 测试"。经核查 `LandNite.f90` 实际有 5-6 个活跃测试（pfMFT、nfMFT、地表温度、BTD 11-12 薄卷云、BTD 11-4 雾/低云、BTD 7.3-11 厚中云、BTD 4-12 薄卷云），CO2 和 6.7μm 测试因通道原因被注释。**问题不是测试数量少，而是所有测试的阈值都未经 MERSI 重校准。**

在通道物理特性已有系统性偏差（见 2.2 节）且阈值未经调整的情况下，继续微调个别阈值的收益有限。特别是：

- **夜间陆地**：3.8μm 替代 3.959μm（BT 偏差 40-90K），所有 BTD 11-4 测试基线偏移
- **极区夜间**：缺少 CO₂ 通道，晴空恢复几乎无效

**结论：夜间/极区模块不是阈值问题，是信息源问题，需要算法层面的替换。**

---

## 三、整体迁移策略

### 3.1 核心原则

1. **迁移即升级**：重写哪个 Fortran 模块，该模块的算法同步升级，不做 1:1 翻译
2. **分阶段替换**：每个阶段结束时，系统可独立运行并与原 Fortran 输出做像素级对比验证
3. **最终目标**：彻底消除 Fortran 依赖，全栈为 C++（计算核心）+ Python（调度/ML 训练）

### 3.2 目标架构

```
Python 层
├── 批处理调度（现有，扩展）
├── ML 模型训练（CALIOP 匹配 + LightGBM/PyTorch）
├── 可视化与验证工具
└── pybind11 接口 ──────────────────────────────┐
                                                 ↓
C++ 核心层                                   [pybind11]
├── IO 模块（HighFive HDF5 + eccodes GRIB2）
├── NWP 插值模块（Eigen）
├── 白天决策树（阈值 YAML 配置，可热更新）
├── 夜间 / 极区 ML 模块（LightGBM C API）
└── RTM 模块（RTTOV C 接口 或 PFAAST C++重写）
```

---

## 四、分阶段执行计划

### 第一阶段：构建基础设施 + IO 迁移

**目标**：建立 CMake 构建体系，替换 IO 和格式转换层，Fortran 核心暂时保留

#### 4.1.1 构建系统迁移（Makefile → CMake）

现有 `build.sh` 调用 `ifort/icc` 编译整体，迁移到 CMake 分库管理：

```cmake
cmake_minimum_required(VERSION 3.20)
project(fylat LANGUAGES Fortran C CXX)

# 暂时保留的 Fortran 核心
add_library(fylat_fortran STATIC
    src/pfaast/pfaast.f90
    src/cloudmask/cloud_core.f90
    # ...其余 Fortran 源文件
)

# 新建 C++ IO 库
find_package(HighFive REQUIRED)
find_package(pybind11 REQUIRED)

add_library(fylat_io SHARED
    src/cpp/io/mersi_l1_reader.cpp
    src/cpp/io/nwp_grib_reader.cpp
    src/cpp/io/clm_writer.cpp
)
target_link_libraries(fylat_io HighFive::HighFive)

# pybind11 Python 模块
pybind11_add_module(fylat_py
    src/python_bindings/bindings.cpp
)
target_link_libraries(fylat_py PRIVATE fylat_io fylat_fortran)
```

#### 4.1.2 HDF5 IO 替换（Fortran → C++）

使用 **HighFive**（header-only C++ HDF5 封装）重写 L1 数据读取和 CLM 输出：

```cpp
// src/cpp/io/mersi_l1_reader.hpp
#include <highfive/H5File.hpp>
#include <Eigen/Dense>

struct MersiL1Data {
    // 25 个通道的反射率 / 亮温
    std::array<Eigen::ArrayXXf, 25> bands;
    // 几何参数
    Eigen::ArrayXXf solar_zenith;
    Eigen::ArrayXXf sensor_zenith;
    Eigen::ArrayXXf relative_azimuth;
    // 元数据
    int nlines, npixels;
    std::string sensor_id;
};

class MersiL1Reader {
public:
    explicit MersiL1Reader(const std::string& filepath);
    MersiL1Data read();
private:
    HighFive::File file_;
};
```

pybind11 绑定后，Python 侧直接拿到 NumPy 数组，零拷贝：

```python
import fylat_py as fylat
import numpy as np

data = fylat.read_mersi_l1("/data/FY3D_MERSI_20230101_0000.HDF")
bt_108 = data.bands[23]  # B24, 10.8μm, shape: (nlines, npixels)
```

#### 4.1.3 GRIB2 解析替换

将 `wgrib/` 脚本链替换为 Python `cfgrib` + `eccodes`：

```python
# python/fylat/nwp_reader.py
import cfgrib
import numpy as np

def read_gfs_grib2(path: str, valid_time: str) -> dict:
    ds = cfgrib.open_dataset(path, filter_by_keys={'typeOfLevel': 'isobaricInhPa'})
    return {
        't': ds['t'].values,      # 温度，shape: (level, lat, lon)
        'q': ds['q'].values,      # 比湿
        'rh': ds['r'].values,     # 相对湿度
    }
```

#### 4.1.4 阶段验证

每个替换模块必须通过像素级对比：

```python
# tests/test_io_regression.py
def test_l1_reader_regression():
    # 用同一景数据，对比 Fortran 输出和 C++ 输出
    fortran_out = load_fortran_reference("testdata/ref_bt108.npy")
    cpp_out = fylat.read_mersi_l1("testdata/FY3D_test.HDF").bands[23]
    np.testing.assert_allclose(cpp_out, fortran_out, rtol=1e-5)
```

**阶段完成标志**：Python 批处理脚本的数据读取全部走 C++ IO 接口，wgrib 脚本链退役。

---

### 第二阶段：NWP 插值 + 白天决策树迁移

**目标**：替换 NWP 插值模块，重写白天决策树（同步解耦硬编码阈值）

#### 4.2.1 NWP 插值（C++ + Eigen）

Fortran 的 101 压力层时空插值重写为 C++，使用 Eigen 做矩阵运算：

```cpp
// src/cpp/nwp/nwp_interpolator.hpp
#include <Eigen/Dense>

class NWPInterpolator {
public:
    struct NWPProfile {
        Eigen::VectorXf pressure_levels;   // 101层压力
        Eigen::ArrayXXXf temperature;      // (level, lat, lon)
        Eigen::ArrayXXXf specific_humidity;
        Eigen::ArrayXXXf relative_humidity;
    };

    // 双线性空间插值 + 线性时间插值（两个预报时次之间）
    NWPProfile interpolate_to_swath(
        const NWPProfile& nwp_t0,
        const NWPProfile& nwp_t1,
        float time_weight,
        const Eigen::ArrayXXf& lat,
        const Eigen::ArrayXXf& lon
    );
};
```

#### 4.2.2 白天决策树（阈值解耦）

白天模块的核心逻辑翻译质量较高（通道信息足够），但需要把所有硬编码阈值提取出来：

**现有问题**（Fortran 中的典型写法）：
```fortran
! 硬编码 MODIS 原始阈值，无法动态调整
if (ref_138 > 0.03) then
    cloud_flag = CLOUDY
end if
```

**C++ 重写后**：
```cpp
// coeff/thresholds_mersi_ii.yaml
// day_land:
//   cirrus_138_threshold: 0.025   # MERSI-II 校准后值，不同于 MODIS 的 0.03
//   vrat_threshold: 1.1
//   btd_11_12_threshold: 0.5

class DayLandDetector {
public:
    explicit DayLandDetector(const ThresholdConfig& cfg) : cfg_(cfg) {}

    CloudConfidence detect(const PixelFeatures& feat) {
        // 薄卷云测试
        if (feat.ref_138 > cfg_.cirrus_138_threshold)
            return CloudConfidence::CLOUDY;
        // 可见光比值测试
        if (feat.vrat > cfg_.vrat_threshold)
            return CloudConfidence::PROBABLY_CLOUDY;
        // BTD 测试
        if (std::abs(feat.btd_11_12) > cfg_.btd_11_12_threshold)
            return CloudConfidence::PROBABLY_CLOUDY;
        return CloudConfidence::CLEAR;
    }

private:
    ThresholdConfig cfg_;
};
```

阈值统一存放在 `coeff/thresholds_mersi_ii.yaml`，支持运行时热加载，**彻底和 MODIS 原始值解耦**。

#### 4.2.3 阶段完成标志

白天场景（陆地/海洋/沙漠/海岸线）全部走 C++ 决策树，Fortran 仅剩夜间/极区模块和 PFAAST RTM。

---

### 第三阶段：夜间 / 极区 ML 模块（算法质变）

**目标**：用 LightGBM 彻底替换夜间和极区 Fortran 模块，解决通道缺失导致的算法上限问题

#### 4.3.1 为什么选 LightGBM 而不是深度学习

| 考量 | LightGBM | 深度学习（CNN/Transformer） |
|------|----------|--------------------------|
| C++ 推理 | 官方 C API，直接加载 `.txt` 模型 | 需要 TensorRT / ONNXRuntime |
| 推理速度 | 极快，单景 <1s | 依赖 GPU，CPU 推理较慢 |
| 训练数据量要求 | 数十万样本即可 | 需要百万级以上 |
| 模型可解释性 | 特征重要性可分析 | 黑盒 |
| 集成复杂度 | 低，模型文件放 coeff/ 即可 | 高 |

对于逐像素的云检测（特征维度低、像素独立），LightGBM 完全够用，且工程集成成本低。

#### 4.3.2 训练数据准备（与代码迁移并行进行）

**数据源**：FY-3D MERSI-II L1 + CALIPSO CALIOP VFM（云垂直特征掩码）

匹配条件：
- 时间差 < 5 分钟
- 空间距离 < 1 km
- 仅选取夜间过境（太阳天顶角 > 90°）

```python
# python/fylat/data_prep/caliop_match.py

def match_mersi_caliop(
    mersi_l1_path: str,
    caliop_vfm_path: str,
    max_time_diff_min: float = 5.0,
    max_spatial_diff_km: float = 1.0
) -> pd.DataFrame:
    """
    返回匹配样本 DataFrame，列包含：
    - bt_38, bt_405, bt_108, bt_12   MERSI-II 红外通道亮温
    - btd_38_108, btd_405_108, btd_108_12  亮温差（含 3.8/4.05 双通道）
    - nwp_t2m, nwp_skt              NWP 温度场
    - nwp_rh500, nwp_rh700          NWP 湿度场（补充 ch22=7.23μm 水汽通道信息）
    - glcm_contrast, glcm_entropy   3×3 窗口纹理特征
    - surface_type                  IGBP 地表类型
    - lat, lon, sensor_zenith
    - label                         CALIOP 云标签（0=晴, 1=云）
    """
```

**特征设计说明**：
- `nwp_rh500` / `nwp_rh700`（对流层中上层湿度）补充 ch22=7.23μm 水汽通道，提供垂直水汽分布信息
- `nwp_skt - bt_108`（皮温与亮温差）辅助低云检测
- `bt_405` / `btd_405_108` 保留 4.05μm 通道信息，与 3.8μm 互补应对 3.959μm 缺失
- GLCM 纹理特征捕捉云场的空间结构，对均匀云场和破碎云场有区分能力

#### 4.3.3 模型训练

```python
# python/fylat/training/train_night_model.py
import lightgbm as lgb
from sklearn.model_selection import StratifiedGroupKFold

FEATURES = [
    'bt_38', 'bt_405', 'bt_108', 'bt_12',
    'btd_38_108', 'btd_405_108', 'btd_108_12',
    'nwp_t2m', 'nwp_skt', 'nwp_rh500', 'nwp_rh700',
    'glcm_contrast', 'glcm_entropy',
    'surface_type', 'sensor_zenith', 'lat'
]

# 按地表类型分层训练，6种场景各有独立模型
SURFACE_TYPES = ['land', 'ocean', 'snow_ice', 'desert', 'polar_land', 'polar_ocean']

params = {
    'objective': 'binary',
    'metric': 'binary_logloss',
    'num_leaves': 63,
    'learning_rate': 0.05,
    'feature_fraction': 0.8,
    'bagging_fraction': 0.8,
    'min_child_samples': 50,
}

for surface in SURFACE_TYPES:
    subset = df[df['surface_type'] == surface]
    train_data = lgb.Dataset(subset[FEATURES], label=subset['label'])
    model = lgb.train(params, train_data, num_boost_round=500)
    model.save_model(f'coeff/lgbm_night_{surface}.txt')
```

#### 4.3.4 C++ 推理集成

```cpp
// src/cpp/cloudmask/night_detector.hpp
#include <lightgbm/c_api.h>
#include <Eigen/Dense>

class NightCloudDetector {
public:
    struct Input {
        Eigen::ArrayXXf bt_38;
        Eigen::ArrayXXf bt_108;
        Eigen::ArrayXXf bt_12;
        Eigen::ArrayXXf nwp_t2m;
        Eigen::ArrayXXf nwp_skt;
        Eigen::ArrayXXf nwp_rh500;
        Eigen::ArrayXXf nwp_rh700;
        Eigen::ArrayXXi surface_type;
        Eigen::ArrayXXf sensor_zenith;
        Eigen::ArrayXXf lat;
    };

    // 加载各场景模型
    explicit NightCloudDetector(const std::string& coeff_dir);

    // 推理，返回云概率图（0~1）
    Eigen::ArrayXXf detect(const Input& in);

    // 按阈值转为置信度等级（0~3）
    Eigen::ArrayXXi to_confidence(
        const Eigen::ArrayXXf& prob,
        float cloudy_thresh = 0.5f,
        float prob_cloudy_thresh = 0.35f,
        float prob_clear_thresh = 0.2f
    );

private:
    std::array<BoosterHandle, 6> models_;  // 6种地表类型各一个模型
    std::string coeff_dir_;
};
```

**阈值配置**同样写入 `coeff/thresholds_mersi_ii.yaml`，可在不重编译的情况下调整云概率分割阈值。

#### 4.3.5 降级保底机制

夜间 ML 模型如果输入数据质量差（NWP 缺失、特殊天气事件），自动降级到保守阈值方法：

```cpp
Eigen::ArrayXXf NightCloudDetector::detect(const Input& in) {
    // 检查 NWP 数据有效性
    if (!nwp_is_valid(in)) {
        // 降级：仅用 BTD 做保守判断
        return threshold_fallback(in);
    }
    return ml_inference(in);
}
```

**阶段完成标志**：夜间/极区模块走 LightGBM 推理，Fortran 仅剩 PFAAST RTM。

---

### 第四阶段：RTM 替换（Fortran 彻底消除）

**目标**：替换 PFAAST，彻底消除 Fortran 依赖

#### 4.4.1 两条路径对比

**路径 A：切换到 RTTOV（推荐）**

RTTOV v13+ 提供 C 接口（`rttov_c_interface.h`），可以从 C++ 直接调用：

优点：
- ECMWF 长期维护，精度高于 PFAAST
- 支持 MERSI-II 的仪器系数文件（需向 ECMWF 申请）
- 有成熟的社区和文档

代价：
- 需要注册申请 license 和系数文件
- 初次集成有一定工作量

**路径 B：PFAAST Fortran → C++ 手工移植**

- 完全自主可控，无外部依赖
- 工作量大（PFAAST 约 10,000 行 Fortran）
- 适合不便使用外部 license 的场景

**建议**：优先尝试路径 A，评估 RTTOV 系数文件获取的可行性，如受阻则走路径 B。

#### 4.4.2 Fortran 消除验收标准

```bash
# 最终状态：src/ 下不再有 .f90 / .f / .f 文件
find src/ -name "*.f90" -o -name "*.f" | wc -l
# 预期输出：0

# CMake 不再需要 Fortran 编译器
# CMakeLists.txt 中 LANGUAGES 仅含 C CXX
```

---

## 五、并行工作安排

下表说明各工作流可以并行推进，不互相阻塞：

| 工作流 | 可并行阶段 | 说明 |
|--------|-----------|------|
| CALIOP 匹配数据积累 | 从现在开始 | 与所有代码阶段并行，等第三阶段用 |
| LightGBM 夜间模型训练 | 第一阶段结束后 | 有匹配数据即可训练，不依赖 C++ 重写 |
| RTTOV license 申请 | 从现在开始 | 申请周期可能较长，早启动 |
| 阈值统计优化 | 第一阶段结束后 | 用 CALIOP 数据对白天阈值做统计优化 |

---

## 六、目录结构演进

### 当前

```
cloudmask_fast/
├── src/                    # Fortran/C 混合，~37,000 行
├── coeff/                  # 阈值系数（混杂 MODIS 原始值）
├── python/fylat/           # 批处理驱动
├── wgrib/                  # NWP 格式转换脚本
└── scripts/
```

### 第一阶段结束

```
cloudmask_fast/
├── src/
│   ├── cpp/
│   │   └── io/             # 新增：C++ IO 模块
│   └── fortran/            # 原有 Fortran，整理归类
├── python/fylat/
│   ├── nwp_reader.py       # 新增：cfgrib 替换 wgrib
│   └── bindings/           # pybind11 绑定层
├── coeff/
└── CMakeLists.txt          # 新增：替换 Makefile
```

### 第三阶段结束

```
cloudmask_fast/
├── src/
│   ├── cpp/
│   │   ├── io/             # HDF5 + GRIB2
│   │   ├── nwp/            # 插值
│   │   ├── cloudmask/
│   │   │   ├── day/        # 白天决策树
│   │   │   └── night/      # 夜间 ML 模块
│   │   └── rtm/            # PFAAST 保留（iso_c_binding）
│   └── fortran/            # 仅剩 PFAAST
├── python/fylat/
│   ├── training/           # ML 模型训练
│   ├── data_prep/          # CALIOP 匹配
│   └── validation/         # 验证工具
└── coeff/
    ├── thresholds_mersi_ii.yaml    # 白天阈值（MERSI-II 专属）
    └── lgbm_night_*.txt            # 夜间 ML 模型（6 个场景）
```

### 第四阶段结束（最终状态）

```
cloudmask_fast/
├── src/cpp/                # 全部为 C++，无 Fortran
├── python/fylat/
└── coeff/
```

---

## 七、验证策略

### 7.1 回归测试（每个阶段必做）

每替换一个模块，必须与原 Fortran 输出做像素级对比，验收标准：

```python
# 白天决策树：允许极少量差异（浮点精度）
assert (cpp_result == fortran_result).mean() > 0.999

# 夜间 ML 模块：与 Fortran 不直接对比，改为与 CALIOP 真值对比
assert overall_accuracy > 0.85
assert pod > 0.80          # 云的命中率（Probability of Detection）
assert far < 0.15          # 晴空误判率（False Alarm Rate）
```

### 7.2 性能基准

```bash
# 对比单景处理时间（目标：C++ 版本不慢于 Fortran，最终应更快）
time ./fylat_fortran config.nml
time ./fylat_cpp config.yaml
```

---

## 八、风险与对策

| 风险 | 概率 | 对策 |
|------|------|------|
| PFAAST C++ 移植工作量超预期 | 高 | 优先申请 RTTOV，两条路并行评估 |
| CALIOP 夜间匹配样本不足（极区） | 中 | 补充 FY-3E 数据，扩大时间窗口 |
| LightGBM 模型在特殊天气泛化差 | 中 | 保留阈值降级兜底，持续更新训练数据 |
| RTTOV license 申请受阻 | 低 | 直接走 PFAAST 手工移植路径 B |

---

## 九、里程碑汇总

| 里程碑 | 内容 | 依赖 |
|--------|------|------|
| M1 | CMake 构建体系就绪，C++ HDF5 IO 可用 | — |
| M2 | wgrib 脚本链退役，GRIB2 走 Python cfgrib | M1 |
| M3 | 白天决策树 C++ 版本通过回归测试，阈值 YAML 化 | M2 |
| M4 | CALIOP 匹配数据集 > 100 万样本，夜间模型训练完成 | M1 并行 |
| M5 | 夜间/极区走 LightGBM 推理，OA > 85% | M3 + M4 |
| M6 | RTM 替换完成，`find src/ -name "*.f90" \| wc -l` 输出为 0 | M5 |

---

*文档维护：随各阶段进展更新，建议每个里程碑完成后同步修订。*