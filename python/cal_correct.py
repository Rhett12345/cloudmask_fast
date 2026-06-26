"""按正确的L1标定链路重新计算所有IR通道BT"""
import h5py, numpy as np, sys

C1 = 1.191042e-5
C2 = 1.4387752
BT_RANGE = {20:(200,350), 21:(200,380), 22:(180,280), 23:(180,300), 24:(180,330), 25:(180,330)}

fname = sys.argv[1] if len(sys.argv) > 1 else \
    '/data/Data_yuq/mersi/20220803/FY3D_MERSI_GBAL_L1_20220803_0740_1000M_MS.HDF'

f = h5py.File(fname, 'r')

# 数据 (RAD0, 已标定辐亮度放大值)
rad0_1km = f['Data/EV_1KM_Emissive'][:]      # [4,2000,2048] B20-23
rad0_agg = f['Data/EV_250_Aggr.1KM_Emissive'][:]  # [2,2000,2048] B24-25

# SDS属性: Slope, Intercept
slope_1km = f['Data/EV_1KM_Emissive'].attrs['Slope'][:]        # [4]
slope_agg = f['Data/EV_250_Aggr.1KM_Emissive'].attrs['Slope'][:]  # [2]
intercept_1km = f['Data/EV_1KM_Emissive'].attrs['Intercept'][:]
intercept_agg = f['Data/EV_250_Aggr.1KM_Emissive'].attrs['Intercept'][:]

# 文件全局属性
wavelengths = f.attrs['Effect_Center_WaveLength'][:]   # [25] all channels
A_coeff = f.attrs['TBB_Trans_Coefficient_A'][:]        # [6] B20-25
B_coeff = f.attrs['TBB_Trans_Coefficient_B'][:]        # [6] B20-25
f.close()

print("Effect_Center_WaveLength[19:25]:", wavelengths[19:25])
print("TBB_Trans_Coefficient_A:", A_coeff)
print("TBB_Trans_Coefficient_B:", B_coeff)

# 标定各通道
results = {}
for ch in [20, 21, 22, 23, 24, 25]:
    if ch <= 23:
        idx = ch - 20
        rad0 = rad0_1km[idx]
        slope = slope_1km[idx]
        intercept = intercept_1km[idx]
    else:
        idx = ch - 24
        rad0 = rad0_agg[idx]
        slope = slope_agg[idx]
        intercept = intercept_agg[idx]

    # Step1: RAD = RAD0 * Slope + Intercept
    rad = rad0.astype(np.float64) * slope + intercept
    rad = np.maximum(rad, 1e-10)

    # Step2: 波数
    wl = wavelengths[ch - 1]  # 0-indexed
    nu = 10000.0 / wl  # cm⁻¹

    # Step3: Planck逆变换 → Te
    Te = C2 * nu / np.log(C1 * nu**3 / rad + 1.0)

    # Step4: 线性修正 → Tbb
    ai = ch - 20  # 0-5
    A = A_coeff[ai]
    B = B_coeff[ai]
    Tbb = A * Te + B

    results[ch] = Tbb.astype(np.float32)

# 统计
clm_file = '/data/Data_yuq/fy3_cloud/20220803/FY3D_MERSI_ORBT_L2_CLM_MLT_NUL_20220803_0740_1000M_MS_BUSINESS.HDF'
with h5py.File(clm_file, 'r') as f:
    clm = f['Cloud_Mask'][:]
byte0 = clm[0].astype(np.uint8)
determined = (byte0 & 0x01) != 0
clear_bits = (byte0 >> 1) & 0x03
is_cloudy = (clear_bits == 0) & determined
is_clear = ((clear_bits == 2) | (clear_bits == 3)) & determined

print("\n" + "=" * 80)
print("CORRECTED BT STATISTICS (RAD0 * Slope → Te → Tbb correction)")
print("=" * 80)

for ch in [20, 21, 22, 23, 24, 25]:
    bt = results[ch]
    valid = ~np.isnan(bt)
    bv = bt[valid]
    bt_lo, bt_hi = BT_RANGE[ch]
    n_over = np.sum(bv > bt_hi)
    n_under = np.sum(bv < bt_lo)
    n_in = np.sum((bv >= bt_lo) & (bv <= bt_hi))
    total = len(bv)

    print(f"\nCH{ch} (wl={wavelengths[ch-1]:.3f}um, nu={10000/wavelengths[ch-1]:.1f}cm-1, A={A_coeff[ch-20]:.5f}, B={B_coeff[ch-20]:.4f})")
    print(f"  BT: [{bv.min():.1f}, {bv.max():.1f}] K, median={np.median(bv):.1f}K")
    print(f"  p5={np.percentile(bv,5):.1f}, p25={np.percentile(bv,25):.1f}, p75={np.percentile(bv,75):.1f}, p95={np.percentile(bv,95):.1f}")
    print(f"  Dynamic range [{bt_lo}, {bt_hi}] K: in={n_in}({100*n_in/total:.1f}%), over={n_over}({100*n_over/total:.1f}%)")

    # 按晴空/云天分
    for lbl, mask in [('clear', is_clear), ('cloudy', is_cloudy)]:
        subset = bv[mask[valid]]
        if len(subset) > 0:
            print(f"  {lbl}: median={np.median(subset):.1f}K, p25={np.percentile(subset,25):.1f}, p75={np.percentile(subset,75):.1f}")

# Key BTD statistics
print("\n" + "=" * 80)
print("KEY BTD STATISTICS (corrected)")
print("=" * 80)

# BT11 - BT3.8 (ch24-ch20)
# BT11 - BT4.05 (ch24-ch21)
# BT8.6 - BT7.2 (ch23-ch22)
# BT7.2 - BT11 (ch22-ch24)
for (ch1, ch2, label) in [(24, 20, 'BT11-BT3.8'), (24, 21, 'BT11-BT4.05'),
                            (23, 22, 'BT8.6-BT7.2'), (22, 24, 'BT7.2-BT11'),
                            (20, 21, 'BT3.8-BT4.05')]:
    v = ~np.isnan(results[ch1]) & ~np.isnan(results[ch2])
    btd = results[ch1][v] - results[ch2][v]
    print(f"\n{label}:")
    for lbl, mask in [('All', slice(None)), ('Clear', is_clear[v]), ('Cloudy', is_cloudy[v])]:
        if lbl == 'All':
            d = btd
        else:
            d = btd[mask]
        if len(d) > 0:
            print(f"  {lbl}: mean={np.mean(d):.1f}K, std={np.std(d):.1f}K, "
                  f"p25={np.percentile(d,25):.1f}, p50={np.median(d):.1f}, p75={np.percentile(d,75):.1f}")

# Threshold pass rates
print("\n" + "=" * 80)
print("THRESHOLD PASS RATES (corrected BT)")
print("=" * 80)
tests = [
    ("Ocean Day 11-4.05  do11_4lo=-8.0 >=", 24, 21, -8.0, ">="),
    ("Ocean Day 11-3.8   do11_4lo=-8.0 >=", 24, 20, -8.0, ">="),
    ("Ocean Nite 11-4.05 no11_4lo=1.0 <=",  24, 21, 1.0,  "<="),
    ("Ocean Nite 11-3.8  no11_4lo=1.0 <=",  24, 20, 1.0,  "<="),
    ("Day Land 11-3.8   dl11_4lo=-12.0 >=", 24, 20, -12.0, ">="),
    ("Day Land 11-4.05  dl11_4lo=-12.0 >=", 24, 21, -12.0, ">="),
    ("Nite Land 11-3.8  nl_11_4l=-2.5 <=",  24, 20, -2.5,  "<="),
    ("Desert 11-3.8 lo  lds11_4lo=-23.0 >=", 24, 20, -23.0, ">="),
    ("Desert 11-3.8 hi  lds11_4hi=-5.0 <=",  24, 20, -5.0,  "<="),
    ("Nite Ocean 8.6-7.2 no86_73=17.0 >",   23, 22, 17.0,  ">"),
    ("Nite Land 7.2-11  nl7_11s=-10.0 <=",   22, 24, -10.0, "<="),
    ("Nite Snow 11-3.8  ns11_4lo=0.6 <=",    24, 20, 0.6,   "<="),
    ("Nite Snow 4-12    ns4_12hi=5.5 <=",    20, 25, 5.5,   "<="),
]
print(f"{'Test':<35} {'Op':>4} {'Thr':>8} {'ClearPass':>10} {'CloudPass':>10}")
print("-" * 75)
for name, ch1, ch2, thr, op in tests:
    v = ~np.isnan(results[ch1]) & ~np.isnan(results[ch2])
    btd = results[ch1][v] - results[ch2][v]
    if op == ">=": passes = btd >= thr
    elif op == "<=": passes = btd <= thr
    elif op == ">": passes = btd > thr
    clear_pass = np.sum(passes & is_clear[v]) / max(np.sum(is_clear[v]), 1) * 100
    cloud_pass = np.sum(passes & is_cloudy[v]) / max(np.sum(is_cloudy[v]), 1) * 100
    print(f"{name:<35} {op:>4} {thr:>8.1f} {clear_pass:>9.1f}% {cloud_pass:>9.1f}%")
