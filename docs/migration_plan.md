# FYLAT 迁移执行计划

> 更新时间：2026-06-30  
> 当前定位：本文件是当前推荐维护的迁移总计划。`goals.md` 和
> `improvement_plan.md` 暂作为历史背景保留，不再作为日常执行入口。

## 1. 当前状态

FYLAT 已经完成 M1/M2 的主体迁移工作，当前处在
“Python 调度 + C++ IO 起步 + Fortran 核心保留”的过渡阶段。

已完成或基本完成：

- `fylat_io` + `fylat_py` 已纳入 CMake，可从源码构建。
- `python/fylat/mersi_io.py` 已接入 C++ HDF5/pybind11 后端。
- `python/fylat/nwp_reader.py` 默认使用 ecCodes 解析 GRIB2。
- `wgrib2` 仅作为显式 fallback，通过 `FYLAT_NWP_BACKEND=wgrib2` 启用。
- 三轨、两套定标回归曾验证输出一致，说明 M1/M2 改动没有改变现有反演结果。
- `scripts/check_migration_build.sh` 已固化远端 `cloudmask` 环境的 CMake/CTest 入口。
- `scripts/compare_clm_regression.py` 已提供三轨、两套定标的 CLM/QA/SHA256 回归比较入口。
- `analyze_accuracy.py` 已升级为分层诊断入口，可输出 day/night、SZA、BT11、海岸/内陆水体、NWP age 与光谱特征统计。

仍需注意：

- 主反演核心仍主要由 Fortran 执行。
- C++ IO 目前偏底层 dataset 读写，还不是完整业务对象。
- `run_fylat.py` 的默认 NWP 日志已切到 Python GRIB2/ecCodes 语义，旧 wgrib2 仅作为 fallback 保留。
- Fortran `read_nwp_data_module.f90` 中仍保留旧脚本调用路径，直接运行 Fortran 时仍可能触发旧逻辑。

当前迁移完成度估计：

| 口径 | 估计完成度 | 说明 |
|------|------------|------|
| 迁移路线图整体 | 约 30% | M1/M2 主体完成，M3 尚未正式落地 |
| Fortran 消除比例 | 约 5-10% | 核心反演、RTM、场景决策树仍主要在 Fortran |

## 2. 下一步优先级

### P0：收尾 M1/M2 文档和入口一致性

- 保持 `run_fylat.py` 的 NWP 文案与 ecCodes 默认实现一致。
- 复用 `generate_fortran_nwp_binary()`，按 Fortran 当前读取的单时次命名生成 NWP binary。
- 明确测试 `FYLAT_NWP_BACKEND=eccodes` 和 `FYLAT_NWP_BACKEND=wgrib2` 两条路径。
- 在 `usage.md` 中继续维护 C++ IO 和 NWP backend 的推荐用法。
- 远端构建默认命令：

```bash
conda activate cloudmask
bash scripts/check_migration_build.sh /tmp/fylat_build_check
```

### P1：把 C++ IO 升级为业务级接口

当前 C++ IO 后端已经能读写基础 HDF5 dataset。下一步应从通用 dataset API
升级为面向 FYLAT 业务的 reader/writer：

- `MersiGeoReader`：读取纬经度、太阳/卫星角、DEM、LSM 等 GEO 数据。
- `MersiL1Reader`：读取并组织 1 km L1B 反射率、亮温和通道元数据。
- `ClmWriter`：封装 CLM/QA 输出，减少 Python 与 Fortran 输出结构差异。

目标是让 Python 侧尽量面对业务对象，而不是手动拼 HDF5 dataset 路径。

### P2：启动 M3，优先迁移 ocean_day

`python/fylat/ocean_day.py` 已有一版接近 Fortran 语义的 Python replica，
因此最适合作为第一块 C++ 决策树迁移对象。

执行顺序：

1. 以 `python/fylat/ocean_day.py` 作为行为基线。
2. 新增 C++ `ocean_day` detector，先追求逐像素一致，而不是马上调优阈值。
3. 阈值统一读取 `coeff/thresholds_mersi_ii.yaml`。
4. 通过 pybind11 暴露实验接口。
5. Python 侧先增加实验入口，不直接替换生产主流程。
6. 回归稳定后，再引入批处理开关。

预留接口：

```bash
FYLAT_CLOUDMASK_BACKEND=fortran        # 默认，继续走现有稳定生产路径
FYLAT_CLOUDMASK_BACKEND=cpp_ocean_day  # 仅 ocean_day 走 C++ 实验路径
FYLAT_CLOUDMASK_BACKEND=auto           # 后续按场景自动选择 C++/Fortran
```

当前 `run_fylat.py` 已校验该环境变量，但非 Fortran 后端仍只作为 M3 占位提示，不改变默认生产结果。

### P3：建立三轨、两套定标自动回归

M1/M2 已经做过三轨、两套定标对比，但应沉淀成可重复测试。

建议固定：

- 三个代表轨道：例如 `0715`、`0720`、`0740`。
- 两套定标：`business` 与 `recali`。
- 对比对象：`Cloud_Mask`、`Quality_Assurance`、文件级 SHA256。
- 验收标准：M1/M2 收尾阶段要求完全一致；M3 算法迁移阶段先要求目标场景逐像素一致。

当前可重复入口：

```bash
python scripts/compare_clm_regression.py --date 20220803 --times 0715,0720,0740
```

若要验证迁移前后完全一致，使用 `--baseline-root` 指向旧输出目录，并加
`--require-identical`。

### P4：夜间、极区与 ML 模块后置

CALIOP/LightGBM 夜间与极区模块依赖训练数据、验证协议和推理接口，不应早于
M3 白天决策树稳定前全面展开。

当前建议只做准备工作：

- 梳理 CALIOP 匹配数据需求。
- 固化特征列表。
- 明确模型验收指标。
- 保留 Fortran 阈值路径作为 fallback。

## 3. 必须保留的文档

以下文档是长期保留文档，不应在迁移清理中删除：

- `docs/usage.md`：运行、配置、定标、数据路径和 backend 使用入口。
- `docs/mersi_ii_channels.md`：FY-3D/E MERSI-II 通道参数。
- `docs/modis_channels.md`：Aqua MODIS / MYD35 对照通道参数。
- `docs/validation_20220803_batch.md`：20220803 批处理验证历史报告。

以下文档暂时保留，但不再作为当前执行入口：

- `docs/goals.md`：早期总体路线和背景分析。
- `docs/improvement_plan.md`：早期算法改进规划和验证记录。

后续如需精简 docs，可新建 `docs/archive/`，再将历史规划移动进去。

## 4. 接口和兼容性约定

当前已经存在并应继续维护：

```bash
FYLAT_IO_BACKEND=auto   # 默认：有 C++ 后端则优先使用，否则回退 h5py
FYLAT_IO_BACKEND=cpp    # 强制使用 C++ HDF5 后端
FYLAT_IO_BACKEND=h5py   # 强制使用 Python h5py 后端

FYLAT_NWP_BACKEND=eccodes  # 默认：Python ecCodes 解析 GRIB2
FYLAT_NWP_BACKEND=wgrib2   # 显式兼容旧 wgrib2 路径
```

未来预留：

```bash
FYLAT_CLOUDMASK_BACKEND=fortran|cpp_ocean_day|auto
```

兼容性原则：

- 默认运行结果必须保持与当前 Fortran 主流程一致。
- 新 C++ 模块先通过显式开关进入实验路径。
- 每个模块替换都必须有像素级回归测试。
- 生产批处理默认路径只能在回归稳定后切换。

## 5. 验收标准

M1/M2 收尾完成标准：

- CMake 可构建 `fylat_py`。
- `ctest --output-on-failure` 全部通过。
- `run_fylat.py` 的 NWP 文案、入口和默认 backend 与 ecCodes 实现一致。
- `usage.md` 中记录当前推荐运行方式。

M3 第一阶段完成标准：

- C++ `ocean_day` detector 可通过 pybind11 调用。
- 同一输入下，C++ 输出与 Python/Fortran 基线逐像素一致。
- 阈值来自 `coeff/thresholds_mersi_ii.yaml`，不在 C++ 中硬编码业务阈值。
- 默认生产路径仍走 Fortran。

长期完成标准：

- 白天场景决策树逐步迁移到 C++。
- 夜间/极区模块在 CALIOP 验证后再切换到 ML 推理。
- RTM 替换完成后，Fortran 依赖最终清零。
