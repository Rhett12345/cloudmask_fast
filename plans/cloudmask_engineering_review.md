# FY-3D MERSI-II Cloud Mask 工程化与规范性审查报告

**审查日期**: 2026-06-24
**审查角度**: 工程化、代码规范、可维护性、项目管理
**基于提交**: `78d1d6f` fix: memory leak, file handle leak, error handling (v0.2.1)

---

## 总体评价

项目从卫星遥感算法原型向可维护的工程化软件过渡中，存在较多历史遗留问题。核心算法逻辑（MOD35 决策树移植）基本正确，但软件工程实践薄弱，主要体现在：二进制文件入库、巨型单体模块、缺乏测试、调试代码残留、命名混乱等方面。

---

## 一、版本控制问题

### 1. 二进制文件提交到仓库 (HIGH)

**问题**: 以下二进制文件被 git 跟踪：

| 文件 | 大小 | 类型 |
|------|------|------|
| `wgrib/wgrib2` | 2.6 MB | ELF 64-bit 可执行文件 |
| `wgrib/wgrib` | 152 KB | ELF 32-bit 可执行文件 |
| `fylat_FY3_MERSI_II_PGS` | 6.5 MB | 编译产物 |

总计约 9.3 MB 二进制文件。

**影响**:
- 仓库体积膨胀
- 每次二进制变更都完整存储，历史不可压缩
- 不同平台的二进制不兼容

**建议**: 从 git 历史中清除，加入 `.gitignore`，提供构建说明替代。

### 2. 备份文件被跟踪 (MEDIUM)

**问题**: `visualize/io_mersi.py.bak` 是备份文件，不应进入版本控制。

**建议**: 加入 `.gitignore`: `*.bak`

### 3. .gitignore 不完整 (MEDIUM)

**当前 `.gitignore`**:
```
*.o
*.mod
fylat_FY3_MERSI_II_PGS
.codegraph/
coeff/fylat_ecosystem.hdf
...
```

**缺失的规则**:
```
*.bak
*.pyc
__pycache__/
*.log
*.nc
*.hdf
*.h5
wgrib/wgrib
wgrib/wgrib2
cal_mode.txt
plans/
```

### 4. 提交历史质量 (LOW)

- 只有 4 次提交，缺乏渐进式开发历史
- 提交信息格式不统一（有的用中文，有的用英文）
- 缺乏 conventional commits 规范

---

## 二、构建系统问题

### 5. clean.sh 引用错误的 Makefile (MEDIUM)

**文件**: `clean.sh`

```bash
cd src
make clean -f fylat_makefile    # 错误！应该是 fylat_makefile_cldmask
cd ..
```

### 6. build.sh 缺少错误处理 (MEDIUM)

**文件**: `build.sh`

```bash
#bshell                    # shebang 被注释掉了？
cd src
make -f fylat_makefile_cldmask
cd ..
```

**问题**:
- 没有 `#!/bin/bash` shebang
- 没有 `set -e` 错误退出
- 没有检查 make 是否成功
- 没有环境检查（ifort/h5fc 是否存在）

**建议**:
```bash
#!/bin/bash
set -euo pipefail

command -v ifort >/dev/null 2>&1 || { echo "ERROR: ifort not found"; exit 1; }
command -v h5fc >/dev/null 2>&1 || { echo "ERROR: h5fc not found"; exit 1; }

cd src
make -f fylat_makefile_cldmask || { echo "BUILD FAILED"; exit 1; }
echo "BUILD SUCCESS"
```

### 7. Makefile 硬编码路径 (MEDIUM)

**文件**: `src/fylat_makefile_cldmask`

```makefile
root_path = /home/liusy2020/yuq/cloudmask/retrieval_system_V3.1_cldmask/src/
hdf5= /opt/intel_lib/bin/h5fc
```

**问题**: 路径只在特定机器上有效，换机器需要手动修改。

**建议**: 使用环境变量或相对路径：
```makefile
root_path = $(shell pwd)/
h5fc ?= h5fc
```

### 8. 缺少依赖管理 (MEDIUM)

**问题**: 没有记录外部依赖的版本和安装方式：
- Intel ifort 编译器版本
- HDF5 库版本和路径
- wgrib2 版本
- Python 包依赖（无 requirements.txt）

**建议**: 添加 `README.md` 中的依赖说明或 `INSTALL.md`。

---

## 三、代码结构问题

### 9. 巨型单体模块 (HIGH)

| 文件 | 行数 | 子程序数 | 评估 |
|------|------|----------|------|
| `io_module.f90` | 3,223 | 48 | 过大，应拆分 |
| `nwp_utils_module.f90` | 3,730 | 110 | 过大，应拆分 |
| `read_nwp_data_module.f90` | 2,770 | 39 | 过大 |
| `get_ancil_data_module.f90` | 2,260 | 117 | 过大 |
| `fylat_fy3mersi_cloud_mask.f90` | 1,322 | 16 | 可接受 |

**问题**:
- 单个文件超过 3000 行，难以理解和维护
- 职责不单一，违反单一职责原则
- 修改一个功能可能影响其他功能

**建议**: 按功能拆分，例如 `io_module.f90` 可拆为：
- `io_geo_module.f90` - GEO 数据读取
- `io_l1b_module.f90` - L1b 数据读取
- `io_cloud_mask_module.f90` - 云掩码输出
- `io_utils_module.f90` - 通用 I/O 工具

### 10. 文件命名不规范 (MEDIUM)

**问题**: 文件名风格混乱：

| 风格 | 示例 |
|------|------|
| 小写下划线 | `land_module.f90`, `water_module.f90` |
| 驼峰混合 | `LandDay.f90`, `PolarNite_ocean.f90` |
| 大写开头 | `Antarctic_day.f90`, `Day_snow.f90` |
| 缩写 | `chk_land.f90`, `tview.f` |

**建议**: 统一为小写下划线风格，如 `land_day.f90`, `polar_nite_ocean.f90`。

### 11. F77/F90 混合 (MEDIUM)

**问题**: 同一功能模块中混合使用 F77（`.f`）和 F90（`.f90`）：
- `.f` 文件使用固定格式、COMMON 块、GOTO
- `.f90` 文件使用自由格式、模块、显式接口

**影响**: 难以维护，风格不一致。

**建议**: 逐步将 `.f` 文件迁移到 `.f90`，使用模块替代 COMMON 块。

---

## 四、代码质量问题

### 12. 调试代码大量残留 (MEDIUM)

**统计**:
- `fylat_fy3mersi_cloud_mask.f90`: 62 处 `jincheng`/`lyj` 标记
- `io_module.f90`: 70 处注释掉的 print/write 语句
- `io_module_intermediate.f90`: 10 处调试标记

**示例**:
```fortran
!print*,'ij',iline,ielem         !jincheng
!open(1,file='test.txt')         !jincheng
!if (iline ==1800 .and. ielem == 1100) then
!   print*,'desert,hi_elev',desert,hi_elev
!endif
```

**影响**: 降低可读性，增加维护负担。

**建议**: 全部删除，需要调试时使用 `-DDEBUG` 预处理宏。

### 13. 大量注释掉的代码 (MEDIUM)

**问题**: 整个子程序、模块引用、COMMON 块被注释掉而非删除：

```fortran
!use io_module
!integer(kind=1) :: lsf
!common / bug / debug, h_output
```

**建议**: 使用 git 历史保留旧代码，当前代码只保留活跃代码。

### 14. SAVE 语句滥用 (LOW)

**问题**: 34 个子程序使用了 `save` 语句。在 F90 模块中，模块变量自动具有 SAVE 属性，子程序中的局部变量默认不保留。过度使用 `save` 可能导致意外的状态保留。

**建议**: 审查每个 `save` 的必要性，移除不必要的。

### 15. INCLUDE vs USE 混用 (MEDIUM)

**问题**: 同一文件中混合使用 `include` 和 `use`：

```fortran
use names_module
use data_arrays_module
include 'global.inc'           ! include 文件定义参数
include 'LandDay_thr.inc'      ! include 文件定义阈值变量
```

**影响**:
- `include` 文件中的变量没有命名空间，容易冲突
- `global.inc` 定义的 `bad_data` 等常量应该放在模块中

**建议**: 将 `global.inc` 中的常量移入 `constant.f90` 模块，阈值变量移入 `thresholds_read_module`。

### 16. 全局状态过多 (MEDIUM)

**问题**: `data_arrays_module.f90` 定义了大量全局变量：

```fortran
type(fylat_fy3_mersi_geo), public  :: geo
type(fylat_fy3_mersi_L1b), public  :: sat
type(nwpdata), save, public :: nwpo
```

所有模块通过 `use data_arrays_module` 共享这些全局状态，形成紧耦合。

**建议**: 通过参数传递数据，减少全局状态。

---

## 五、测试与验证问题

### 17. 缺乏自动化测试 (HIGH)

**问题**:
- 没有单元测试
- 没有集成测试
- 没有回归测试
- 没有 CI/CD 配置

**影响**: 修改代码后无法验证正确性，依赖人工检查。

**建议**:
- 为核心数值函数（`conf_test`, `tview`, `planck` 等）编写单元测试
- 为端到端流程编写集成测试（使用小规模测试数据）
- 使用 CMake + CTest 或自定义测试框架

### 18. 缺乏输入验证 (MEDIUM)

**问题**: 程序对输入数据的有效性检查不足：
- Namelist 参数没有范围检查
- HDF5 文件路径没有预先验证
- 通道索引没有边界检查

**建议**: 在程序入口添加参数验证。

---

## 六、配置管理问题

### 19. 配置文件散乱 (MEDIUM)

**问题**: 项目根目录下有 12 个 `.nml` 配置文件和 4 个 `.xcfg` 文件：

```
cfg_20200308_1345_business.nml
cfg_20200308_1345_recali.nml
cfg_20200308_1435_business.nml
cfg_20200308_1435_recali.nml
...
VIS_Cal_Coeff_business.xcfg
VIS_Cal_Coeff_recali.xcfg
VIS_Cal_Coeff.xcfg
VIS_Cal_Coeff_old.xcfg
```

**建议**:
- 配置文件移入 `config/` 目录
- 使用模板配置 + 参数覆盖机制
- 删除重复/过时的配置

### 20. cal_mode.txt 硬编码文件名 (LOW)

**文件**: `src/io_module.f90:433`

```fortran
inquire(file='cal_mode.txt', exist=vis_flag)
```

**问题**: 使用相对路径，依赖工作目录。且 `cal_mode.txt` 被 git 跟踪（含 `recali`），不应进入版本控制。

**建议**: 通过 namelist 参数传入，或使用 `trim(code_root_path)//'cal_mode.txt'`。

---

## 七、文档问题

### 21. README 不完整 (MEDIUM)

**问题**: README 文件存在两个版本（`README.md` 英文，`README_zh.md` 中文），但：
- 缺少安装步骤
- 缺少依赖说明
- 缺少使用示例
- 缺少算法原理简述

### 22. 代码注释风格不统一 (LOW)

**问题**: 注释风格混乱：
- F77 风格: `C`, `!`, `!!`
- Doxygen 风格: `!!F90`, `!!Description:`
- 中文注释: `!lyj`, `!jincheng test`
- 英文注释: `! revised by minmin`

**建议**: 统一使用 `!` 注释，关键子程序使用 Doxygen 风格文档。

---

## 八、建议的改进路线图

### 阶段 1: 立即清理（1-2 天）
1. 从 git 中移除二进制文件（`wgrib/wgrib`, `wgrib/wgrib2`, `fylat_FY3_MERSI_II_PGS`）
2. 完善 `.gitignore`
3. 修复 `clean.sh`
4. 删除 `.bak` 文件

### 阶段 2: 构建系统改进（2-3 天）
5. `build.sh` 添加 shebang、错误处理、环境检查
6. Makefile 改用相对路径或环境变量
7. 添加 `INSTALL.md` 依赖说明

### 阶段 3: 代码清理（3-5 天）
8. 删除所有调试代码残留（`jincheng`, `lyj`）
9. 删除注释掉的代码块
10. 将 `global.inc` 常量移入模块

### 阶段 4: 结构重构（1-2 周）
11. 拆分巨型模块（`io_module`, `nwp_utils_module`）
12. 统一文件命名风格
13. 将 `.f` 文件迁移到 `.f90`
14. 减少全局状态

### 阶段 5: 质量保障（持续）
15. 添加单元测试框架
16. 为核心函数编写测试
17. 设置 CI/CD 流程
18. 代码审查流程

---

## 附录：问题清单

| # | 类别 | 级别 | 问题 | 文件 |
|---|------|------|------|------|
| 1 | 版本控制 | HIGH | 二进制文件提交到仓库 | wgrib/*, fylat_* |
| 2 | 版本控制 | MEDIUM | .bak 文件被跟踪 | visualize/io_mersi.py.bak |
| 3 | 版本控制 | MEDIUM | .gitignore 不完整 | .gitignore |
| 4 | 版本控制 | LOW | 提交历史质量 | - |
| 5 | 构建系统 | MEDIUM | clean.sh 引用错误 Makefile | clean.sh |
| 6 | 构建系统 | MEDIUM | build.sh 缺少错误处理 | build.sh |
| 7 | 构建系统 | MEDIUM | Makefile 硬编码路径 | fylat_makefile_cldmask |
| 8 | 构建系统 | MEDIUM | 缺少依赖管理 | - |
| 9 | 代码结构 | HIGH | 巨型单体模块 | io_module, nwp_utils_module |
| 10 | 代码结构 | MEDIUM | 文件命名不规范 | 多个文件 |
| 11 | 代码结构 | MEDIUM | F77/F90 混合 | .f 文件 |
| 12 | 代码质量 | MEDIUM | 调试代码残留 | 62+ 处 |
| 13 | 代码质量 | MEDIUM | 注释掉的代码 | 多处 |
| 14 | 代码质量 | LOW | SAVE 语句滥用 | 34 个子程序 |
| 15 | 代码质量 | MEDIUM | INCLUDE vs USE 混用 | 多个文件 |
| 16 | 代码质量 | MEDIUM | 全局状态过多 | data_arrays_module |
| 17 | 测试 | HIGH | 缺乏自动化测试 | - |
| 18 | 测试 | MEDIUM | 缺乏输入验证 | Driver |
| 19 | 配置 | MEDIUM | 配置文件散乱 | 根目录 |
| 20 | 配置 | LOW | cal_mode.txt 硬编码 | io_module.f90 |
| 21 | 文档 | MEDIUM | README 不完整 | README.md |
| 22 | 文档 | LOW | 注释风格不统一 | 多个文件 |
