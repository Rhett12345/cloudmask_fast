# CLM 诊断与改进执行记录

> 基线：20220803，40 条 MYD35 匹配轨道，约 1.18 亿像元。

## 当前基线

- Overall HSS 约 0.607，云/晴二分类已经可用。
- Overall agreement 约 68.5%，四级 CLM 置信等级仍弱。
- BUSINESS 与 RECALI 的 HSS 基本一致，说明白天效果一般不主要由可见光定标系数造成。
- 重点问题层：
  - `Mid 30-60 Land`：FAR 约 24%，是陆地白天/中纬误报主战场。
  - `High 60-90 Land`：HSS 约 0.19，极区陆地/雪冰背景仍不可靠。
  - `Coast/InlandWater`：混合像元、湖岸和海岸线应独立看，不再混入普通水体或陆地。

## 诊断命令

```bash
conda activate cloudmask
python analyze_accuracy.py --date 20220803 --output accuracy_20220803.json
```

快速模式：

```bash
python analyze_accuracy.py --date 20220803 --no-spectral-diagnostics
```

输出会同时包含：

- `recal` / `onboard`：二分类 HSS、POD、FAR、CSI 与四级 class agreement。
- `feature_stats`：按 MYD35 cloud/clear 标签聚合的 BT11、BTD11-12、BTD8-11、0.65/0.86/1.38 μm 特征均值与标准差。
- day/night、SZA、BT11 温度段、海岸/内陆水体、NWP forecast age 等诊断层。

## 改进优先级

1. 先用 `feature_stats.DayMidLand` 和 `feature_stats.MidLand` 重标定陆地白天阈值，目标是压低中纬陆地 FAR，且 Overall HSS 不退化。
2. 单独处理 `Coast/InlandWater`，必要时给湖岸、海岸和小水体更保守的 clear 恢复逻辑。
3. 对 `High 60-90 Land` 引入雪冰/低温专项 QA、空间纹理和 NWP 温度背景一致性诊断，不强行只靠 BTD。
4. 单独校准 0/1/2/3 置信等级映射；二分类 HSS 与四分类 agreement 分开验收。

## 验收

每次阈值或决策树调整后至少报告：

- Overall：HSS、POD、FAR、class agreement。
- Day SZA<80。
- Mid 30-60 Land。
- High 60-90 Land。
- Coast/InlandWater。
- BUSINESS 与 RECALI 差异。

默认生产路径必须保持可回退：`FYLAT_CLOUDMASK_BACKEND=fortran`。
