import h5py, numpy as np

C1 = 1.191042e-5; C2 = 1.4387752
WNUM = {20:2643.44, 21:2471.65, 22:1382.62, 23:1168.18, 24:933.36, 25:836.94}

f = h5py.File('/data/Data_yuq/mersi/20220803/FY3D_MERSI_GBAL_L1_20220803_0740_1000M_MS.HDF', 'r')
dn1 = f['Data/EV_1KM_Emissive'][:]
dn2 = f['Data/EV_250_Aggr.1KM_Emissive'][:]
cal = f['Calibration/IR_Cal_Coeff'][:]

print("Slope attributes from L1:")
for name in ['Data/EV_1KM_Emissive', 'Data/EV_250_Aggr.1KM_Emissive']:
    ds = f[name]
    s = ds.attrs['Slope'][:]
    i = ds.attrs['Intercept'][:]
    bn = ds.attrs['band_name']
    print(f"  {name}: slope={s}, intercept={i}, bands={bn}")

f.close()

print("\nIR_Cal_Coeff k1 vs Slope_attr ratio:")
slope_attr = {20:0.0002, 21:0.0002, 22:0.01, 23:0.01, 24:0.01, 25:0.01}
for ch in [20,21,22,23,24,25]:
    k1 = float(np.mean(cal[ch-20, 1, :]))
    k0 = float(np.mean(cal[ch-20, 0, :]))
    k2 = float(np.mean(cal[ch-20, 2, :]))
    ratio = k1 / slope_attr[ch]
    print(f"  CH{ch}: k1={k1:.6f}, slope_attr={slope_attr[ch]:.4f}, k1/k1_simple={ratio:.2f}x, k0={k0:.4f}")

print("\nRadiance at median DN vs Planck at expected T:")
T_exp = {20:300, 21:300, 22:265, 23:295, 24:295, 25:295}
for ch in [20,21,22,23,24,25]:
    ci = ch - 20
    dn_arr = dn1[ci] if ch <= 23 else dn2[ch-24]
    v = (dn_arr > 0) & (dn_arr < 65535)
    dn_med = float(np.median(dn_arr[v]))
    k0 = float(np.mean(cal[ci, 0, :]))
    k1 = float(np.mean(cal[ci, 1, :]))
    k2 = float(np.mean(cal[ci, 2, :]))
    rad = k0 + k1*dn_med + k2*dn_med**2
    w = WNUM[ch]
    B = C1 * w**3 / (np.exp(C2 * w / T_exp[ch]) - 1.0)
    bt = C2 * w / np.log(1.0 + C1 * w**3 / max(rad, 1e-10))
    print(f"  CH{ch}: rad_cal={rad:.2f}, B({T_exp[ch]}K)={B:.2f}, ratio={rad/B:.2f}x → BT={bt:.1f}K")
