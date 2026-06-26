"""基于实际BTD分布计算最优阈值，输出建议的阈值调整方案"""
import h5py, numpy as np

C1 = 1.191042e-5; C2 = 1.4387752
BT_RANGE = {20:(200,350), 21:(200,380), 22:(180,280), 23:(180,300), 24:(180,330), 25:(180,330)}

# ---- 读数据 ----
f = h5py.File('/data/Data_yuq/mersi/20220803/FY3D_MERSI_GBAL_L1_20220803_0740_1000M_MS.HDF', 'r')
rad0_1km = f['Data/EV_1KM_Emissive'][:]
rad0_agg = f['Data/EV_250_Aggr.1KM_Emissive'][:]
slope_1km = f['Data/EV_1KM_Emissive'].attrs['Slope'][:]
slope_agg = f['Data/EV_250_Aggr.1KM_Emissive'].attrs['Slope'][:]
wavelengths = f.attrs['Effect_Center_WaveLength'][:]
A = f.attrs['TBB_Trans_Coefficient_A'][:]
B = f.attrs['TBB_Trans_Coefficient_B'][:]
f.close()

# 标定
bt = {}
for ch in [20,21,22,23,24,25]:
    if ch <= 23:
        rad0 = rad0_1km[ch-20]; slope = slope_1km[ch-20]
    else:
        rad0 = rad0_agg[ch-24]; slope = slope_agg[ch-24]
    vmask = (rad0 > 0) & (rad0 < 65533)
    rad = rad0.astype(np.float64) * slope
    rad = np.maximum(rad, 1e-10)
    wl = wavelengths[ch-1]; nu = 10000.0 / wl
    Te = C2 * nu / np.log(C1 * nu**3 / rad + 1.0)
    bt[ch] = A[ch-20] * Te + B[ch-20]

# 读云掩码
clm = h5py.File('/data/Data_yuq/fy3_cloud/20220803/FY3D_MERSI_ORBT_L2_CLM_MLT_NUL_20220803_0740_1000M_MS_BUSINESS.HDF', 'r')['Cloud_Mask'][:]
byte0 = clm[0].astype(np.uint8)
determined = (byte0 & 0x01) != 0
clear_bits = (byte0 >> 1) & 0x03
is_cloudy = (clear_bits == 0) & determined
is_clear = ((clear_bits == 2) | (clear_bits == 3)) & determined

# ---- 计算最优阈值 ----
def find_optimal_threshold(btd_clear, btd_cloudy, op):
    """扫描阈值，找最大化 Youden's J = TPR - FPR 的点"""
    all_vals = np.concatenate([btd_clear, btd_cloudy])
    candidates = np.linspace(np.percentile(all_vals, 1), np.percentile(all_vals, 99), 200)

    best_j, best_thr = -1, 0
    for thr in candidates:
        if op == '>=':
            tpr = np.mean(btd_clear >= thr)    # clear passes
            fpr = np.mean(btd_cloudy >= thr)   # cloudy passes
        elif op == '<=':
            tpr = np.mean(btd_clear <= thr)
            fpr = np.mean(btd_cloudy <= thr)
        elif op == '>':
            tpr = np.mean(btd_clear > thr)
            fpr = np.mean(btd_cloudy > thr)
        j = tpr - fpr
        if j > best_j:
            best_j = j
            best_thr = thr

    return best_thr, best_j

def eval_threshold(btd_clear, btd_cloudy, thr, op):
    """Evaluate a specific threshold"""
    if op == '>=':
        tpr = np.mean(btd_clear >= thr) * 100
        fpr = np.mean(btd_cloudy >= thr) * 100
    elif op == '<=':
        tpr = np.mean(btd_clear <= thr) * 100
        fpr = np.mean(btd_cloudy <= thr) * 100
    elif op == '>':
        tpr = np.mean(btd_clear > thr) * 100
        fpr = np.mean(btd_cloudy > thr) * 100
    return tpr, fpr, tpr - fpr

# 定义测试及其 BTD 计算
# (name, ch1, ch2, current_thr, op, scene_filter_description)
tests = [
    # 11-4 类测试
    ("do11_4lo",   24, 21, -8.0,  ">=", "海洋白天 11-4.05"),
    ("do11_4lo",   24, 20, -8.0,  ">=", "海洋白天 11-3.8"),
    ("no11_4lo",   24, 21,  1.0,  "<=", "海洋夜间 11-4.05"),
    ("no11_4lo",   24, 20,  1.0,  "<=", "海洋夜间 11-3.8"),
    ("dl11_4lo",   24, 20, -12.0, ">=", "陆地白天 11-3.8"),
    ("dl11_4lo",   24, 21, -12.0, ">=", "陆地白天 11-4.05"),
    ("nl_11_4l",   24, 20, -2.5,  "<=", "陆地夜间 11-3.8"),
    ("lds11_4lo",  24, 20, -23.0, ">=", "沙漠 11-3.8 low"),
    ("lds11_4hi",  24, 20, -5.0,  "<=", "沙漠 11-3.8 high"),
    ("ns11_4lo",   24, 20,  0.6,  "<=", "雪面夜间 11-3.8"),
    # 水汽测试
    ("no86_73",    23, 22, 17.0,  ">",  "海洋夜间 8.6-7.2"),
    ("nl7_11s",    22, 24, -10.0, "<=", "陆地夜间 7.2-11"),
    # 4-12 测试
    ("ns4_12hi",   20, 25,  5.5,  "<=", "雪面夜间 4-12"),
]

print("=" * 85)
print(f"{'参数':<14} {'通道':<10} {'当前阈值':>8} {'最优阈值':>9} {'当前区分':>8} {'最优区分':>8}")
print(f"{'':14} {'':10} {'':8} {'(midpt)':>9} {'(TPR-FPR)':>8} {'(TPR-FPR)':>8}")
print("-" * 85)

suggestions = {}
for name, ch1, ch2, cur_thr, op, desc in tests:
    v = (bt[ch1] > 0) & (bt[ch2] > 0) & ~np.isnan(bt[ch1]) & ~np.isnan(bt[ch2])
    btd = bt[ch1][v] - bt[ch2][v]
    btd_c = btd[is_clear[v]]
    btd_ld = btd[is_cloudy[v]]

    cur_tpr, cur_fpr, cur_sep = eval_threshold(btd_c, btd_ld, cur_thr, op)
    opt_thr, opt_sep = find_optimal_threshold(btd_c, btd_ld, op)

    # 置信区间参数: lo-cut = midpt + range/2, hi-cut = midpt - range/2
    # 取 clear p10 和 cloudy p90 来定置信区间宽度
    if op in ('>=', '>'):
        lo_cut = np.percentile(btd_c, 5)   # clear 低分位 → 0% confidence
        hi_cut = np.percentile(btd_c, 95)  # clear 高分位 → 100% confidence
    else:
        lo_cut = np.percentile(btd_c, 95)
        hi_cut = np.percentile(btd_c, 5)
    power = 1.0

    key = (name, f"ch{ch1}-ch{ch2}")
    suggestions[key] = {
        'desc': desc, 'cur_thr': cur_thr, 'opt_thr': round(opt_thr, 1),
        'cur_sep': round(cur_sep, 1), 'opt_sep': round(opt_sep, 1),
        'lo_cut': round(lo_cut, 1), 'hi_cut': round(hi_cut, 1),
        'op': op
    }

    marker = "←" if abs(opt_thr - cur_thr) > 1.5 else " "
    print(f"{name:<14} {'ch'+str(ch1)+'-ch'+str(ch2):<10} {cur_thr:>8.1f} {opt_thr:>9.1f} {cur_sep:>8.1f}% {opt_sep:>8.1f}% {marker}")

print("\n" + "=" * 85)
print("建议的阈值调整方案")
print("=" * 85)
print()
for key, s in suggestions.items():
    name, chpair = key
    if abs(s['opt_thr'] - s['cur_thr']) > 1.5:
        print(f"  {name} ({s['desc']}, {chpair}):")
        print(f"    {s['cur_thr']:.1f} → {s['opt_thr']:.1f}  (midpt), 区分度 {s['cur_sep']:.1f}% → {s['opt_sep']:.1f}%")
        print(f"    置信区间: [{s['lo_cut']:.1f}, {s['hi_cut']:.1f}], power=1.0")
        print()

# --- 按系数文件格式打印建议 ---
print("=" * 85)
print("coeff/fylat_thresholds.mersi.ii3d.v8 中需修改的行")
print("=" * 85)
threshold_map = {
    ('do11_4lo', 'ch24-ch21'): 'do11_4lo',
    ('no11_4lo', 'ch24-ch21'): 'no11_4lo',
    ('dl11_4lo', 'ch24-ch20'): 'dl11_4lo',
    ('lds11_4lo', 'ch24-ch20'): 'lds11_4lo',
    ('lds11_4hi', 'ch24-ch20'): 'lds11_4hi',
    ('ns11_4lo', 'ch24-ch20'): 'ns11_4lo',
    ('no86_73', 'ch23-ch22'): 'no86_73',
    ('nl7_11s', 'ch22-ch24'): 'nl7_11s',
    ('ns4_12hi', 'ch20-ch25'): 'ns4_12hi',
}
for key, coeff_name in threshold_map.items():
    s = suggestions[key]
    fmt = f"{coeff_name:<15}: {s['lo_cut']:.1f}, {s['opt_thr']:.1f}, {s['hi_cut']:.1f}, 1.0"
    cur = f"{coeff_name:<15}: {s['lo_cut']:.1f}, {s['cur_thr']:.1f}, {s['hi_cut']:.1f}, 1.0"
    if abs(s['opt_thr'] - s['cur_thr']) > 1.5:
        print(f"  OLD: {cur}")
        print(f"  NEW: {fmt}")
        print()
