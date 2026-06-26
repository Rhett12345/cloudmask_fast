"""诊断所有IR通道BT标定，对比Planck预期值"""
import h5py
import numpy as np
import sys

C1 = 1.191042e-5
C2 = 1.4387752

WNUM = {20:2643.44, 21:2471.65, 22:1382.62, 23:1168.18, 24:933.36, 25:836.94}
BT_RANGE = {20:(200,350), 21:(200,380), 22:(180,280), 23:(180,300), 24:(180,330), 25:(180,330)}
NAME = {20:'3.8um', 21:'4.05um', 22:'7.2um', 23:'8.55um', 24:'10.7um', 25:'12.0um'}

fname = sys.argv[1] if len(sys.argv) > 1 else \
    '/data/Data_yuq/mersi/20220803/FY3D_MERSI_GBAL_L1_20220803_0740_1000M_MS.HDF'

f = h5py.File(fname, 'r')
dn_1km = f['Data/EV_1KM_Emissive'][:]
dn_agg = f['Data/EV_250_Aggr.1KM_Emissive'][:]
cal = f['Calibration/IR_Cal_Coeff'][:]
f.close()

print("=" * 75)
print("PLANCK EXPECTED RADIANCE AT TYPICAL BT")
print("=" * 75)
for ch in [20, 21, 22, 23, 24, 25]:
    w = WNUM[ch]
    T = {20:300, 21:300, 22:270, 23:270, 24:300, 25:300}[ch]
    B = C1 * w**3 / (np.exp(C2 * w / T) - 1.0)
    print(f"  CH{ch}({NAME[ch]}): Planck at {T}K = {B:.4f} mW/(m2 sr cm-1)")

print("\n" + "=" * 75)
print("ACTUAL CALIBRATION RESULTS (before dynamic range clipping)")
print("=" * 75)

for ch in [20, 21, 22, 23, 24, 25]:
    ci = ch - 20
    dn = dn_1km[ci] if ch <= 23 else dn_agg[ch - 24]
    k0 = float(np.mean(cal[ci, 0, :]))
    k1 = float(np.mean(cal[ci, 1, :]))
    k2 = float(np.mean(cal[ci, 2, :]))
    k3 = float(np.mean(cal[ci, 3, :]))

    dn_f = dn.astype(np.float64)
    rad = k0 + k1 * dn_f + k2 * dn_f**2 + k3 * dn_f**3
    rad = np.maximum(rad, 1e-10)

    # Also compute using simple slope attribute
    slope_attr = 0.0
    if ch <= 23:
        slope_attr = 0.0002 if ch <= 21 else 0.01
    else:
        slope_attr = 0.01

    rad_simple = slope_attr * dn_f
    rad_simple = np.maximum(rad_simple, 1e-10)

    w = WNUM[ch]
    bt = C2 * w / np.log(1.0 + C1 * w**3 / rad)
    bt_simple = C2 * w / np.log(1.0 + C1 * w**3 / rad_simple)

    v = (dn > 0) & (dn < 65535)
    bv = bt[v]
    bv_s = bt_simple[v]

    bt_lo, bt_hi = BT_RANGE[ch]
    n_over = np.sum(bv > bt_hi)
    n_under = np.sum(bv < bt_lo)
    n_ok = np.sum((bv >= bt_lo) & (bv <= bt_hi))

    dn_med = np.median(dn[v])
    rad_med = k0 + k1 * float(dn_med) + k2 * float(dn_med)**2

    print(f"\n  CH{ch} ({NAME[ch]}): range [{bt_lo}, {bt_hi}] K")
    print(f"    cal coef: k0={k0:.6f}, k1={k1:.6f}, k2={k2:.2e}")
    print(f"    DN: median={dn_med:.0f}, min={dn[v].min()}, max={dn[v].max()}")
    print(f"    rad at median DN: {rad_med:.3f}")
    print(f"    BT (IR_Cal_Coeff): median={np.median(bv):.1f}K, p5={np.percentile(bv,5):.1f}, p95={np.percentile(bv,95):.1f}")
    print(f"    BT (Slope={slope_attr}): median={np.median(bv_s):.1f}K, p5={np.percentile(bv_s,5):.1f}, p95={np.percentile(bv_s,95):.1f}")
    print(f"    in range: {n_ok} ({100*n_ok/len(bv):.1f}%), >upper: {n_over} ({100*n_over/len(bv):.1f}%), <lower: {n_under}")

# Key check: BT23
print("\n" + "=" * 75)
print("CH23 (8.55um) DIAGNOSTIC - IS CALIBRATION TOO WARM?")
print("=" * 75)

ci = 23 - 20
dn23 = dn_1km[ci]
k0 = float(np.mean(cal[ci, 0, :]))
k1 = float(np.mean(cal[ci, 1, :]))
k2 = float(np.mean(cal[ci, 2, :]))

v = (dn23 > 0) & (dn23 < 65535)

# What radiance does Planck predict at 290K, 300K, 310K for 8.55um?
w23 = WNUM[23]
for T in [270, 280, 290, 300, 310, 320, 330]:
    B = C1 * w23**3 / (np.exp(C2 * w23 / T) - 1.0)
    # What DN would produce this radiance?
    # rad = k0 + k1*DN + k2*DN^2 -> solve quadratic for DN
    # k2*DN^2 + k1*DN + (k0 - B) = 0
    a = k2
    b = k1
    c = k0 - B
    if abs(a) > 1e-12:
        disc = b**2 - 4*a*c
        if disc >= 0:
            dn_est = (-b + np.sqrt(disc)) / (2*a)
        else:
            dn_est = -1
    else:
        dn_est = (B - k0) / k1
    print(f"  T={T}K -> Planck B={B:.3f} -> DN(est)={dn_est:.0f}")

print(f"\n  Actual CH23 DN: median={np.median(dn23[v]):.0f}, mean={dn23[v].mean():.1f}")

# DN should be roughly proportional to radiance. If BT23 is ~320K,
# Planck B(320K) should match the radiance at median DN
B320 = C1 * w23**3 / (np.exp(C2 * w23 / 320) - 1.0)
B300 = C1 * w23**3 / (np.exp(C2 * w23 / 300) - 1.0)
B290 = C1 * w23**3 / (np.exp(C2 * w23 / 290) - 1.0)

dn_med = float(np.median(dn23[v]))
rad_med = k0 + k1 * dn_med + k2 * dn_med**2
print(f"  Actual rad at median DN({dn_med:.0f}): {rad_med:.3f}")
print(f"  Planck B(290K)={B290:.3f}, B(300K)={B300:.3f}, B(320K)={B320:.3f}")
print(f"  Ratio rad_med/B(290K)={rad_med/B290:.2f}")
print(f"  This suggests BT23 is ~{C2*w23/np.log(1+C1*w23**3/max(rad_med,1e-10)):.0f}K")

# Compare with BT22
print("\n" + "=" * 75)
print("BT22 vs BT23 CROSS-VALIDATION")
print("=" * 75)
for ch in [22, 23]:
    ci = ch - 20
    dn = dn_1km[ci]
    k0 = float(np.mean(cal[ci, 0, :]))
    k1 = float(np.mean(cal[ci, 1, :]))
    k2 = float(np.mean(cal[ci, 2, :]))
    dn_f = dn.astype(np.float64)
    rad = np.maximum(k0 + k1*dn_f + k2*dn_f**2, 1e-10)
    bt = C2 * WNUM[ch] / np.log(1.0 + C1 * WNUM[ch]**3 / rad)
    v = (dn > 0) & (dn < 65535)
    bv = bt[v]
    print(f"  CH{ch}({NAME[ch]}): BT median={np.median(bv):.1f}K, p25={np.percentile(bv,25):.1f}, p75={np.percentile(bv,75):.1f}")
    print(f"    DN: median={np.median(dn[v]):.0f}, min={dn[v].min()}, max={dn[v].max()}")

# BTD compute
v_both = (dn_1km[0] > 0) & (dn_1km[0] < 65535) & (dn_1km[2] > 0) & (dn_1km[2] < 65535) & (dn_1km[3] > 0) & (dn_1km[3] < 65535)
for ch in [22, 23]:
    ci = ch - 20
    dn = dn_1km[ci]
    k0 = float(np.mean(cal[ci, 0, :]))
    k1 = float(np.mean(cal[ci, 1, :]))
    k2 = float(np.mean(cal[ci, 2, :]))
    dn_f = dn.astype(np.float64)
    rad = np.maximum(k0 + k1*dn_f + k2*dn_f**2, 1e-10)
    bt = C2 * WNUM[ch] / np.log(1.0 + C1 * WNUM[ch]**3 / rad)
    if ch == 22: bt22_arr = bt
    if ch == 23: bt23_arr = bt

btd_23_22 = bt23_arr[v_both] - bt22_arr[v_both]
print(f"\n  BT8.6-BT7.2 (from IR_Cal_Coeff):")
print(f"    mean={np.mean(btd_23_22):.1f}K, std={np.std(btd_23_22):.1f}K")
print(f"    p5={np.percentile(btd_23_22,5):.1f}, p25={np.percentile(btd_23_22,25):.1f}")
print(f"    p50={np.median(btd_23_22):.1f}, p75={np.percentile(btd_23_22,75):.1f}, p95={np.percentile(btd_23_22,95):.1f}")
print(f"  MODIS expected: BT29(8.55um)~290-300K, BT28(7.3um)~260-275K, BTD~15-25K clear, ~25-40K cloudy")
print(f"  Conclusion: CH23 cal likely has +20~30K systematic warm bias")
PYEOF
