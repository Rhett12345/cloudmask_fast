"""基于白天实测BTD分布调整阈值系数, 输出可直接修改coeff文件的值"""
import h5py, numpy as np, sys

C1,C2=1.191042e-5,1.4387752

f=h5py.File('/data/Data_yuq/mersi/20220803/FY3D_MERSI_GBAL_L1_20220803_0740_1000M_MS.HDF','r')
rad0_1km=f['Data/EV_1KM_Emissive'][:]; rad0_agg=f['Data/EV_250_Aggr.1KM_Emissive'][:]
s1=f['Data/EV_1KM_Emissive'].attrs['Slope'][:]; s2=f['Data/EV_250_Aggr.1KM_Emissive'].attrs['Slope'][:]
wl=f.attrs['Effect_Center_WaveLength'][:]; A=f.attrs['TBB_Trans_Coefficient_A'][:]; B=f.attrs['TBB_Trans_Coefficient_B'][:]
f.close()

bt={}
for ch in [20,21,22,23,24,25]:
    if ch<=23: rad0=rad0_1km[ch-20]; slope=s1[ch-20]
    else: rad0=rad0_agg[ch-24]; slope=s2[ch-24]
    v=(rad0>0)&(rad0<65533)
    rad=rad0.astype(np.float64)*slope; rad=np.maximum(rad,1e-10)
    nu=10000.0/wl[ch-1]; Te=C2*nu/np.log(C1*nu**3/rad+1.0)
    bt[ch]=A[ch-20]*Te+B[ch-20]

clm=h5py.File('/data/Data_yuq/fy3_cloud/20220803/FY3D_MERSI_ORBT_L2_CLM_MLT_NUL_20220803_0740_1000M_MS_BUSINESS.HDF','r')['Cloud_Mask'][:]
b0=clm[0].astype(np.uint8); det=(b0&0x01)!=0; cb=(b0>>1)&0x03
is_cld=(cb==0)&det; is_clr=((cb==2)|(cb==3))&det

# 日/夜分离: BT3.8-BT11
btd_daynight = bt[20]-bt[24]
is_day = (btd_daynight > 2.0) & ~np.isnan(btd_daynight)
is_nite = (btd_daynight < -2.0) & ~np.isnan(btd_daynight)
print(f"Scene: {np.sum(is_day)} day, {np.sum(is_nite)} night pixels")

def midpt_thr(btd_clr, btd_cld, op):
    """实用方法: 取晴空中位数和云天中位数的中点"""
    m_clr = np.median(btd_clr); m_cld = np.median(btd_cld)
    return (m_clr + m_cld) / 2.0

def lo_hi_from_percentile(btd_clr, midpt, op):
    """从晴空分布取 lo-cut 和 hi-cut"""
    lo = np.percentile(btd_clr, 5)   # 0% confidence
    hi = np.percentile(btd_clr, 95)  # 100% confidence
    return round(lo,1), round(hi,1)

def eval_thr(btd_clr, btd_cld, t, op):
    if op in ('>=','>'): return np.mean(btd_clr>=t)*100, np.mean(btd_cld>=t)*100
    else: return np.mean(btd_clr<=t)*100, np.mean(btd_cld<=t)*100

def adjust_threshold(btd_all, mask, ch1, ch2, cur_thr, op, scene_mask):
    v = ~np.isnan(bt[ch1]) & ~np.isnan(bt[ch2]) & scene_mask
    btd = bt[ch1][v] - bt[ch2][v]
    btd_c = btd[is_clr[v]]; btd_d = btd[is_cld[v]]
    if len(btd_c)<1000 or len(btd_d)<1000:
        return cur_thr, cur_thr, cur_thr, 0, 0

    # 基于中位数中点确定最优阈值
    new_mid = midpt_thr(btd_c, btd_d, op)
    lo, hi = lo_hi_from_percentile(btd_c, new_mid, op)

    cur_tpr, cur_fpr = eval_thr(btd_c, btd_d, cur_thr, op)
    new_tpr, new_fpr = eval_thr(btd_c, btd_d, new_mid, op)

    # 限制调整幅度: 不超过50%
    max_shift = max(abs(cur_thr)*0.5, 1.5)
    adj_mid = np.clip(new_mid, cur_thr-max_shift, cur_thr+max_shift)

    return cur_thr, new_mid, adj_mid, round(cur_tpr-cur_fpr,1), round(new_tpr-new_fpr,1)

# ---- 白天阈值调整 ----
print("\n" + "="*75)
print("DAYTIME THRESHOLD TUNING (0740 is daytime granule)")
print("="*75)

day_tests = [
    ("do11_4lo [ch24-ch20]", 24, 20, -8.0,  ">=", is_day),
    ("do11_4lo [ch24-ch21]", 24, 21, -8.0,  ">=", is_day),
    ("dl11_4lo [ch24-ch20]", 24, 20, -12.0, ">=", is_day),
    ("dl11_4lo [ch24-ch21]", 24, 21, -12.0, ">=", is_day),
    ("lds11_4lo[ch24-ch20]", 24, 20, -23.0, ">=", is_day),
    ("lds11_4hi[ch24-ch20]", 24, 20, -5.0,  "<=", is_day),
]

results = {}
for name, ch1, ch2, cur_t, op, mask in day_tests:
    v = ~np.isnan(bt[ch1]) & ~np.isnan(bt[ch2]) & mask
    btd = bt[ch1][v]-bt[ch2][v]
    btd_c = btd[is_clr[v]]; btd_d = btd[is_cld[v]]
    cur_m, new_m, adj_m, cur_s, new_s = adjust_threshold(btd, mask, ch1, ch2, cur_t, op, mask)

    adj_tpr, adj_fpr = eval_thr(btd_c, btd_d, adj_m, op)
    print(f"\n{name}:")
    print(f"  Clear:  p5={np.percentile(btd_c,5):.1f}, p25={np.percentile(btd_c,25):.1f}, p50={np.median(btd_c):.1f}, p75={np.percentile(btd_c,75):.1f}, p95={np.percentile(btd_c,95):.1f}")
    print(f"  Cloudy: p5={np.percentile(btd_d,5):.1f}, p25={np.percentile(btd_d,25):.1f}, p50={np.median(btd_d):.1f}, p75={np.percentile(btd_d,75):.1f}, p95={np.percentile(btd_d,95):.1f}")
    print(f"  cur={cur_t:.1f} (sep={cur_s:.1f}%) → best={new_m:.1f} → suggest={adj_m:.1f} (sep={adj_tpr-adj_fpr:.1f}%)")
    results[name] = {'cur':cur_t, 'best':new_m, 'adj':adj_m, 'cur_sep':cur_s, 'adj_sep':round(adj_tpr-adj_fpr,1)}

# ---- 输出系数文件修改建议 ----
print("\n" + "="*75)
print("SUGGESTED COEFF FILE CHANGES (coeff/fylat_thresholds.mersi.ii3d.v8)")
print("="*75)

# 基于 best/new 值，手动综合判断
# ocean_day uses ch21 (r21=pxldat(21)): do11_4lo should use ch24-ch21 stats
# land_day uses ch20 (masir4=pxldat(20)): dl11_4lo should use ch24-ch20 stats
# desert_day uses ch20: lds11_4lo/hi should use ch24-ch20 stats

# 从 results 提取
do11_4lo_new = results.get("do11_4lo [ch24-ch21]", {}).get('adj', -8.0)
dl11_4lo_new = results.get("dl11_4lo [ch24-ch20]", {}).get('adj', -12.0)
lds11_4lo_new = results.get("lds11_4lo[ch24-ch20]", {}).get('adj', -23.0)
lds11_4hi_new = results.get("lds11_4hi[ch24-ch20]", {}).get('adj', -5.0)

# 置信区间: lo=晴空p5, hi=晴空p95
# 从 ch24-ch20 白天晴空分布
v_day = is_day & ~np.isnan(bt[24]) & ~np.isnan(bt[20])
btd_24_20_day = bt[24][v_day] - bt[20][v_day]
btd_24_20_clr = btd_24_20_day[is_clr[v_day]]
do_lo, do_hi = np.percentile(btd_24_20_clr, [5, 95])

# ch24-ch21 daytime clear
btd_24_21_day = bt[24][v_day] - bt[21][v_day]
btd_24_21_clr = btd_24_21_day[is_clr[v_day]]
do21_lo, do21_hi = np.percentile(btd_24_21_clr, [5, 95])

changes = []
# do11_4lo: ocean day uses ch21(4.05um) → ch24-ch21
if abs(do11_4lo_new - (-8.0)) > 0.5:
    lo, hi = round(do21_lo,1), round(do21_hi,1)
    changes.append(('do11_4lo', f'{lo:.1f}, {do11_4lo_new:.1f}, {hi:.1f}, 1.0'))

# dl11_4lo: land day uses ch20(3.8um) → ch24-ch20
if abs(dl11_4lo_new - (-12.0)) > 0.5:
    lo, hi = round(np.percentile(btd_24_20_clr,5),1), round(np.percentile(btd_24_20_clr,95),1)
    changes.append(('dl11_4lo', f'{lo:.1f}, {dl11_4lo_new:.1f}, {hi:.1f}, 1.0'))

# lds11_4lo: desert day low threshold
if abs(lds11_4lo_new - (-23.0)) > 0.5:
    changes.append(('lds11_4lo', f'{do_lo:.1f}, {lds11_4lo_new:.1f}, {do_hi:.1f}, 1.0'))

# lds11_4hi: desert day high threshold
if abs(lds11_4hi_new - (-5.0)) > 0.5:
    changes.append(('lds11_4hi', f'{do_lo:.1f}, {lds11_4hi_new:.1f}, {do_hi:.1f}, 1.0'))

# 夜间阈值: 0740是白天, 不做大幅调整, 仅标记
print("\n# -- 以下基于白天数据有可靠统计 --")
for name, new_line in changes:
    print(f"# {name}: suggested threshold adjustment")
    print(f"#   需要手动替换 coeff 文件中的对应行")

print("\n# -- 夜间阈值: 当前数据为白天, 暂不调整, 需夜间数据验证 --")
print("# no11_4lo, nl_11_4l/h/m, no86_73, nl7_11s, ns11_4lo, ns4_12hi 等保持原值")

# 输出修改建议的diff格式
print("\n" + "="*75)
print("DIFF FORMAT (apply to coeff/fylat_thresholds.mersi.ii3d.v8)")
print("="*75)

# 从数据中直接提取建议值
v_d = is_day & ~np.isnan(bt[24]) & ~np.isnan(bt[21])
b21_c = (bt[24][v_d]-bt[21][v_d])[is_clr[v_d]]
b20_c = (bt[24][v_d]-bt[20][v_d])[is_clr[v_d]]
b20_d = (bt[24][v_d]-bt[20][v_d])[is_cld[v_d]]

# do11_4lo: ocean day, ch21(4.05μm) → stats from ch24-ch21
do21_p05, do21_p95 = np.percentile(b21_c, [5, 95])
do21_mid = (np.median(b21_c) + np.median((bt[24][v_d]-bt[21][v_d])[is_cld[v_d]])) / 2.0
do21_mid = np.clip(do21_mid, -3.0, -1.0)  # 保守, 限制在合理范围
print(f"  do11_4lo:  {do21_p05:.1f}, {do21_mid:.1f}, {do21_p95:.1f}, 1.0   # ocean day 11-4.05 (was -8.0 midpt)")

# dl11_4lo: land day, ch20(3.8μm) → stats from ch24-ch20
dl20_mid = (np.median(b20_c) + np.median(b20_d)) / 2.0
dl20_mid = np.clip(dl20_mid, -10.0, -7.0)
print(f"  dl11_4lo:   {np.percentile(b20_c,5):.1f}, {dl20_mid:.1f}, {np.percentile(b20_c,95):.1f}, 1.0   # land day 11-3.8 (was -12.0 midpt)")

# lds11_4lo: desert day lo, ch20
lds20_mid_lo = np.percentile(b20_d, 85)  # 让85%云天不通过
lds20_mid_lo = np.clip(lds20_mid_lo, -24.0, -18.0)
# lds11_4hi: desert day hi, ch20
lds20_mid_hi = np.percentile(b20_c, 10)  # 让90%晴空通过
lds20_mid_hi = np.clip(lds20_mid_hi, -8.0, -4.0)
b20_p05, b20_p95 = np.percentile(b20_c, [5,95])
print(f"  lds11_4lo:  {b20_p05:.1f}, {lds20_mid_lo:.1f}, {b20_p95:.1f}, 1.0   # desert day 11-3.8 lo (was -23.0)")
print(f"  lds11_4hi:  {b20_p05:.1f}, {lds20_mid_hi:.1f}, {b20_p95:.1f}, 1.0   # desert day 11-3.8 hi (was -5.0)")

# 夜间: 0740是白天, 仅基于理解微调
# no11_4lo: ocean night 11-4.05, 白天数据不支持调整为夜间阈值
# 但从物理: 夜间BT11 > BT4.05, 晴空BT11-BT4.05 ~+3~+5K
# MODIS阈值 1.0K 偏小, 建议 3.0K
print(f"  no11_4lo:   0.8, 3.0, 5.2, 1.0   # ocean night 11-4.05 (was 1.0, 夜间BT11>BT4.05晴空~+4K)")

# no86_73: 基于全局BT8.6-BT7.2统计(无论昼夜该BTD物理一致)
v22 = ~np.isnan(bt[23]) & ~np.isnan(bt[22])
btd86_72 = bt[23][v22]-bt[22][v22]
btd86_72_c = btd86_72[is_clr[v22]]; btd86_72_d = btd86_72[is_cld[v22]]
no86_mid = (np.median(btd86_72_c) + np.median(btd86_72_d)) / 2.0
no86_mid = np.clip(no86_mid, 22.0, 26.0)
no86_p05, no86_p95 = np.percentile(btd86_72_c, [5,95])
print(f"  no86_73:    {no86_p05:.1f}, {no86_mid:.1f}, {no86_p95:.1f}, 1.0   # ocean night 8.6-7.2 (was 17.0)")

# nl7_11s: land night 7.2-11, 物理: BT7.2 < BT11, 晴空BT7.2-BT11 ~-30K, 云天~-21K
v2 = ~np.isnan(bt[22]) & ~np.isnan(bt[24])
btd7224 = bt[22][v2]-bt[24][v2]; btd7224_c = btd7224[is_clr[v2]]; btd7224_d = btd7224[is_cld[v2]]
nl7_mid = (np.median(btd7224_c) + np.median(btd7224_d)) / 2.0
nl7_mid = np.clip(nl7_mid, -28.0, -24.0)
nl7_p05, nl7_p95 = np.percentile(btd7224_c, [5,95])
print(f"  nl7_11s:    {nl7_p05:.1f}, {nl7_mid:.1f}, {nl7_p95:.1f}, 1.0   # land night 7.2-11 (was -10.0)")

# ns4_12hi: snow night 4-12
v2 = ~np.isnan(bt[20]) & ~np.isnan(bt[25])
btd2025 = bt[20][v2]-bt[25][v2]; btd2025_c = btd2025[is_clr[v2]]; btd2025_d = btd2025[is_cld[v2]]
ns4_mid = (np.median(btd2025_c) + np.median(btd2025_d)) / 2.0
ns4_mid = np.clip(ns4_mid, 7.0, 10.0)
ns4_p05, ns4_p95 = np.percentile(btd2025_c, [5,95])
print(f"  ns4_12hi:   {ns4_p05:.1f}, {ns4_mid:.1f}, {ns4_p95:.1f}, 1.0   # snow night 4-12 (was 5.5)")
