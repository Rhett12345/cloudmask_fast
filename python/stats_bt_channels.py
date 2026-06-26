"""统计 MERSI-II ch20(3.8μm), ch21(4.05μm), ch22(7.2μm) BT分布并与MYD35对比.

Usage:
  python stats_bt_channels.py --date 20220803 --time 0740
"""

import argparse
import os
import sys
from collections import defaultdict
import h5py
import numpy as np

# Planck constants
C1 = 1.191042e-5  # mW / (m^2 sr cm^-4)
C2 = 1.4387752     # K / cm^-1

# MERSI-II IR channel central wavenumbers (cm^-1) from platform_module.f90 (sensor_id=21/22)
WNUM = {
    20: 2643.44,  # 3.79599um
    21: 2471.65,  # 4.04587um
    22: 1382.62,  # 7.23264um
    23: 1168.18,  # 8.56031um
    24: 933.36,   # 10.7139um
    25: 836.94,   # 11.94827um
}

# MERSI-II IR channel valid BT dynamic range (from instrument specs)
BT_VALID_RANGE = {
    20: (200.0, 350.0),   # 3.8um
    21: (200.0, 380.0),   # 4.05um
    22: (180.0, 280.0),   # 7.2um
    23: (180.0, 330.0),   # 8.55um (estimated)
    24: (180.0, 340.0),   # 10.8um (estimated)
    25: (180.0, 340.0),   # 12.0um (estimated)
}

# MODIS equivalent central wavenumbers
# MODIS 3.959um -> 2525.9, MODIS 7.325um -> 1365.2
WNUM_MODIS_39 = 2525.9

BAND_NAMES = {20: "3.8um(ch20)", 21: "4.05um(ch21)", 22: "7.2um(ch22)"}


def read_mersi_l1(filename):
    """Read MERSI-II L1 IR channels and calibrate to BT.

    According to HDF5 band_name attributes:
      EV_1KM_Emissive:  bands 20, 21, 22, 23
      EV_250_Aggr.1KM_Emissive: bands 24, 25

    IR_Cal_Coeff shape: (6, 4, 200) for bands 20-25
      coeff[b, 0, :] = k0 (intercept)
      coeff[b, 1, :] = k1 (linear slope)
      coeff[b, 2, :] = k2 (quadratic, near zero)
      coeff[b, 3, :] = k3 (cubic, zero)
      Radiance = k0 + k1*DN + k2*DN^2 + k3*DN^3

    Returns:
        dict with ch20..ch25 BT arrays (shape: 2000x2048)
    """
    with h5py.File(filename, 'r') as f:
        dn_1km = f['Data/EV_1KM_Emissive'][:]  # (4, 2000, 2048): bands 20,21,22,23
        dn_ch20 = dn_1km[0]
        dn_ch21 = dn_1km[1]
        dn_ch22 = dn_1km[2]
        dn_ch23 = dn_1km[3]

        dn_agg = f['Data/EV_250_Aggr.1KM_Emissive'][:]  # (2, 2000, 2048): bands 24,25
        dn_ch24 = dn_agg[0]
        dn_ch25 = dn_agg[1]

        cal = f['Calibration/IR_Cal_Coeff'][:]  # (6, 4, 200)

    # IR_Cal_Coeff[0]=ch20, [1]=ch21, ..., [5]=ch25
    dns = {20: dn_ch20, 21: dn_ch21, 22: dn_ch22, 23: dn_ch23, 24: dn_ch24, 25: dn_ch25}
    results = {}

    for ch in [20, 21, 22, 23, 24, 25]:
        cal_idx = ch - 20
        # Radiance = k0 + k1*DN + k2*DN^2 + k3*DN^3; k2/k3 are essentially zero
        k0 = float(np.mean(cal[cal_idx, 0, :]))
        k1 = float(np.mean(cal[cal_idx, 1, :]))
        k2 = float(np.mean(cal[cal_idx, 2, :]))
        k3 = float(np.mean(cal[cal_idx, 3, :]))

        dn = dns[ch]
        rad = k0 + k1 * dn.astype(np.float64) + k2 * dn.astype(np.float64)**2 + k3 * dn.astype(np.float64)**3
        rad = np.maximum(rad, 1e-10)
        wnum = WNUM[ch]
        bt = C2 * wnum / np.log(1.0 + C1 * wnum**3 / rad)

        # Clip to instrument valid dynamic range
        bt_min, bt_max = BT_VALID_RANGE[ch]
        bt = np.clip(bt, bt_min, bt_max)

        results[ch] = bt.astype(np.float32)

    return results


def read_fy3d_cloudmask(filename):
    """Read FY-3D cloud mask L2 HDF5 output.

    Returns:
        cloud_mask: (6, 2000, 2048) uint8
        qa: (10, 2000, 2048) uint8
    """
    with h5py.File(filename, 'r') as f:
        clm = f['Cloud_Mask'][:]  # (6, 2000, 2048)
        qa = f['Quality_Assurance'][:]  # (10, 2000, 2048)
    return clm, qa


def get_cloud_flag(clm):
    """Extract cloudy/clear flag from cloud mask byte array.

    Cloud_Mask[0] bit layout (MODIS MOD35 convention):
      bit0: Cloud Mask Flag (0=undetermined, 1=determined)
      bit1-2: Unobstructed FOV Quality (00=cloudy, 01=uncertain, 10=probClear, 11=confClear)

    Returns:
        is_cloudy: (2000, 2048) bool, True if cloudy
        is_clear: (2000, 2048) bool, True if clear (probClear or confClear)
    """
    byte0 = clm[0].astype(np.uint8)

    # Cloud mask determined?
    determined = (byte0 & 0x01) != 0

    # Bits 1-2: clear sky confidence
    clear_bits = (byte0 >> 1) & 0x03
    # 00=cloudy, 01=uncertain, 10=probably clear, 11=confident clear
    is_cloudy_conf = (clear_bits == 0) & determined
    is_uncertain = (clear_bits == 1) & determined
    is_prob_clear = (clear_bits == 2) & determined
    is_conf_clear = (clear_bits == 3) & determined

    is_clear = is_prob_clear | is_conf_clear

    return is_cloudy_conf, is_clear, is_uncertain


def read_myd35_cloudmask(filename):
    """Read MYD35 HDF4 cloud mask.

    Returns:
        clm: (6, 2030, 1354) uint8
    """
    from pyhdf.SD import SD
    f = SD(filename)
    clm = f.select('Cloud_Mask')[:]  # (6, 2030, 1354)
    f.end()
    return clm


def compute_stats(bt_data, clm, ch, label):
    """Compute statistics for a channel's BT under clear/cloudy conditions."""
    is_cloudy, is_clear, is_uncertain = get_cloud_flag(clm)
    bt = bt_data[ch]

    # Filter valid BT
    valid = ~np.isnan(bt)

    stats = {}
    for mask, name in [(is_cloudy, 'cloudy'), (is_clear, 'clear'),
                        (is_uncertain, 'uncertain')]:
        subset = bt[valid & mask]
        if len(subset) > 0:
            stats[name] = {
                'count': len(subset),
                'mean': float(np.mean(subset)),
                'std': float(np.std(subset)),
                'p5': float(np.percentile(subset, 5)),
                'p25': float(np.percentile(subset, 25)),
                'p50': float(np.percentile(subset, 50)),
                'p75': float(np.percentile(subset, 75)),
                'p95': float(np.percentile(subset, 95)),
                'min': float(np.min(subset)),
                'max': float(np.max(subset)),
            }
        else:
            stats[name] = {'count': 0}

    return stats


def compute_btd_stats(bt_data, ch1, ch2, clm, label):
    """Compute BT difference statistics (ch1-ch2)."""
    is_cloudy, is_clear, is_uncertain = get_cloud_flag(clm)

    bt1 = bt_data[ch1]
    bt2 = bt_data[ch2]
    valid = ~np.isnan(bt1) & ~np.isnan(bt2)
    btd = bt1 - bt2

    stats = {}
    for mask, name in [(is_cloudy, 'cloudy'), (is_clear, 'clear'),
                        (is_uncertain, 'uncertain')]:
        subset = btd[valid & mask]
        if len(subset) > 0:
            stats[name] = {
                'count': len(subset),
                'mean': float(np.mean(subset)),
                'std': float(np.std(subset)),
                'p5': float(np.percentile(subset, 5)),
                'p25': float(np.percentile(subset, 25)),
                'p50': float(np.percentile(subset, 50)),
                'p75': float(np.percentile(subset, 75)),
                'p95': float(np.percentile(subset, 95)),
            }
        else:
            stats[name] = {'count': 0}

    return stats


def compare_myd35_same_region(fy3d_clm, myd35_clm, bt_data, ch):
    """Compare FY3D clear/cloudy BT distributions with MYD35 over overlapping region."""
    # This requires spatial matching; simplified version:
    # Just report BT statistics from both FY3D and MYD35 classifications
    # separately for the same approximate region

    fy_is_cloudy, fy_is_clear, _ = get_cloud_flag(fy3d_clm)
    myd_is_cloudy, myd_is_clear, _ = get_cloud_flag(myd35_clm)

    bt = bt_data[ch]
    valid = ~np.isnan(bt)

    # Where they agree
    agree_clear = fy_is_clear & myd_is_clear
    agree_cloudy = fy_is_cloudy & myd_is_cloudy
    # Where they disagree
    fy_clear_myd_cloudy = fy_is_clear & myd_is_cloudy
    fy_cloudy_myd_clear = fy_is_cloudy & myd_is_clear

    result = {}
    for name, mask in [('agree_clear', agree_clear), ('agree_cloudy', agree_cloudy),
                        ('FYclear_MYDcloudy', fy_clear_myd_cloudy),
                        ('FYcloudy_MYDclear', fy_cloudy_myd_clear)]:
        subset = bt[valid & mask]
        if len(subset) > 0:
            result[name] = {
                'count': len(subset),
                'mean': float(np.mean(subset)),
                'std': float(np.std(subset)),
                'p50': float(np.percentile(subset, 50)),
            }
        else:
            result[name] = {'count': 0}

    return result


def print_stats_table(stats, title):
    """Pretty-print statistics table."""
    print(f"\n{'='*80}")
    print(f"  {title}")
    print(f"{'='*80}")

    has_range = 'min' in next(iter(stats.values()))
    if has_range:
        header = f"{'Category':<15} {'Count':>10} {'Mean':>8} {'Std':>8} {'P5':>8} {'P25':>8} {'P50':>8} {'P75':>8} {'P95':>8} {'Min':>8} {'Max':>8}"
    else:
        header = f"{'Category':<15} {'Count':>10} {'Mean':>8} {'Std':>8} {'P5':>8} {'P25':>8} {'P50':>8} {'P75':>8} {'P95':>8}"
    print(header)
    print("-" * len(header))

    for name, s in stats.items():
        if s.get('count', 0) > 0:
            if has_range:
                print(f"{name:<15} {s['count']:>10d} {s['mean']:>8.2f} {s['std']:>8.2f} "
                      f"{s['p5']:>8.2f} {s['p25']:>8.2f} {s['p50']:>8.2f} "
                      f"{s['p75']:>8.2f} {s['p95']:>8.2f} {s['min']:>8.2f} {s['max']:>8.2f}")
            else:
                print(f"{name:<15} {s['count']:>10d} {s['mean']:>8.2f} {s['std']:>8.2f} "
                      f"{s['p5']:>8.2f} {s['p25']:>8.2f} {s['p50']:>8.2f} "
                      f"{s['p75']:>8.2f} {s['p95']:>8.2f}")
        else:
            print(f"{name:<15} {'N/A':>10}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--date', default='20220803')
    parser.add_argument('--time', default='0740')
    parser.add_argument('--l1-dir', default='/data/Data_yuq/mersi')
    parser.add_argument('--clm-dir', default='/data/Data_yuq/fy3_cloud')
    parser.add_argument('--myd35-dir', default='/data/Data_yuq/aqua_modis/MYD35_L2')
    args = parser.parse_args()

    l1_file = os.path.join(args.l1_dir, args.date,
                           f'FY3D_MERSI_GBAL_L1_{args.date}_{args.time}_1000M_MS.HDF')
    clm_file = os.path.join(args.clm_dir, args.date,
                            f'FY3D_MERSI_ORBT_L2_CLM_MLT_NUL_{args.date}_{args.time}_1000M_MS_BUSINESS.HDF')
    clm_file_recali = os.path.join(args.clm_dir, args.date,
                                   f'FY3D_MERSI_ORBT_L2_CLM_MLT_NUL_{args.date}_{args.time}_1000M_MS_RECALI.HDF')

    print(f"L1 file: {l1_file}")
    print(f"CLM (business): {clm_file}")
    print(f"CLM (recali):   {clm_file_recali}")

    if not os.path.exists(l1_file):
        print(f"ERROR: L1 file not found: {l1_file}")
        sys.exit(1)

    # Read L1 and calibrate
    print("\nReading MERSI L1 and calibrating IR channels...")
    bt_data = read_mersi_l1(l1_file)

    # ============ Part 1: FY-3D self statistics ============
    for cal_name, clm_path in [('business', clm_file), ('recali', clm_file_recali)]:
        if not os.path.exists(clm_path):
            print(f"WARNING: {clm_path} not found, skipping")
            continue

        print(f"\n{'#'*80}")
        print(f"#  CALIBRATION: {cal_name}")
        print(f"{'#'*80}")

        clm, qa = read_fy3d_cloudmask(clm_path)

        # Basic cloud mask stats
        is_cloudy, is_clear, is_uncertain = get_cloud_flag(clm)
        total_valid = np.sum(is_cloudy | is_clear | is_uncertain)
        print(f"\nTotal valid pixels: {total_valid}")
        print(f"  Cloudy:    {np.sum(is_cloudy):>10d} ({100*np.sum(is_cloudy)/max(total_valid,1):.1f}%)")
        print(f"  Uncertain: {np.sum(is_uncertain):>10d} ({100*np.sum(is_uncertain)/max(total_valid,1):.1f}%)")
        print(f"  Clear:     {np.sum(is_clear):>10d} ({100*np.sum(is_clear)/max(total_valid,1):.1f}%)")

        # Stats for each channel
        for ch in [20, 21, 22]:
            stats = compute_stats(bt_data, clm, ch, BAND_NAMES[ch])
            print_stats_table(stats, f"BT {BAND_NAMES[ch]} [{cal_name}]")

        # Key BTD statistics
        btd_stats_11_20 = compute_btd_stats(bt_data, 24, 20, clm, "BT11-BT3.8")  # 24=ch24=10.8um
        print_stats_table(btd_stats_11_20, f"BT11 - BT3.8 [ch24-ch20] [{cal_name}]")

        btd_stats_11_21 = compute_btd_stats(bt_data, 24, 21, clm, "BT11-BT4.05")  # 24=ch24=10.8um, 21=4.05um
        print_stats_table(btd_stats_11_21, f"BT11 - BT4.05 [ch24-ch21] [{cal_name}]")

        btd_stats_8_22 = compute_btd_stats(bt_data, 23, 22, clm, "BT8.6-BT7.2")  # 23=8.55um, 22=7.2um
        print_stats_table(btd_stats_8_22, f"BT8.6 - BT7.2 [ch23-ch22] [{cal_name}]")

        btd_stats_22_11 = compute_btd_stats(bt_data, 22, 24, clm, "BT7.2-BT11")  # 22=7.2um, 24=10.8um
        print_stats_table(btd_stats_22_11, f"BT7.2 - BT11 [ch22-ch24] [{cal_name}]")

        btd_stats_20_21 = compute_btd_stats(bt_data, 20, 21, clm, "BT3.8-BT4.05")  # ch20-ch21
        print_stats_table(btd_stats_20_21, f"BT3.8 - BT4.05 [ch20-ch21] [{cal_name}]")

    # ============ Part 2: BT range validation ============
    print(f"\n{'#'*80}")
    print(f"#  BT RANGE VALIDATION (dynamic range clipping applied)")
    print(f"{'#'*80}")

    for ch in [20, 21, 22]:
        bt = bt_data[ch]
        bt_min, bt_max = BT_VALID_RANGE[ch]
        nan_mask = np.isnan(bt)
        valid = ~nan_mask
        vals = bt[valid]
        total = len(vals)

        # Count pixels at clip boundaries (saturated)
        at_max = np.sum(np.abs(vals - bt_max) < 0.01)
        at_min = np.sum(np.abs(vals - bt_min) < 0.01)
        in_range = total - at_min - at_max

        if total > 0:
            pct_at_max = 100 * at_max / total
            pct_at_min = 100 * at_min / total
            print(f"\n  {BAND_NAMES[ch]}: dynamic range [{bt_min}, {bt_max}] K")
            print(f"    Total valid: {total}")
            print(f"    At upper bound (saturated): {at_max} ({pct_at_max:.1f}%)")
            print(f"    At lower bound (clipped):    {at_min} ({pct_at_min:.1f}%)")
            print(f"    Within range:                {in_range} ({100*pct_at_min-pct_at_max:.1f}% unused)")
            print(f"    Range of in-range values:    [{np.min(vals):.1f}, {np.max(vals):.1f}] K")
            print(f"    Median: {np.median(vals):.1f} K,  P5: {np.percentile(vals,5):.1f},  P95: {np.percentile(vals,95):.1f}")
            # Separate clear/cloudy at boundaries
            is_cloudy, is_clear, _ = get_cloud_flag(read_fy3d_cloudmask(clm_file)[0])
            for lbl, mask in [('clear', is_clear), ('cloudy', is_cloudy)]:
                sat = np.sum((np.abs(bt[mask & valid] - bt_max) < 0.01))
                print(f"    {lbl} at upper bound: {sat} ({100*sat/max(np.sum(mask&valid),1):.1f}%)")

    # ============ Part 3: Threshold pass rate check ============
    print(f"\n{'#'*80}")
    print(f"#  THRESHOLD PASS RATE CHECKS (using business calibration)")
    print(f"{'#'*80}")

    clm_biz, _ = read_fy3d_cloudmask(clm_file)
    is_cloudy, is_clear, _ = get_cloud_flag(clm_biz)

    bt11 = bt_data[24]  # ch24 = 10.8um
    bt20 = bt_data[20]
    bt21 = bt_data[21]
    bt22 = bt_data[22]
    bt23 = bt_data[23]  # 8.55um

    valid_all = ~np.isnan(bt11) & ~np.isnan(bt20) & ~np.isnan(bt21) & ~np.isnan(bt22)

    # Test thresholds from mersi.ii3d.v8 coefficient file
    tests = [
        # (name, condition_expr, expected_clear, expected_cloudy)
        ("Day Ocean 11-4.05  do11_4lo(2)=-8.0", bt11 - bt21, -8.0, ">="),
        ("Day Ocean 11-3.8   do11_4lo(2)=-8.0", bt11 - bt20, -8.0, ">="),
        ("Nite Ocean 11-4.05 no11_4lo(2)=1.0", bt11 - bt21, 1.0, "<="),
        ("Nite Ocean 11-3.8  no11_4lo(2)=1.0", bt11 - bt20, 1.0, "<="),
        ("Day Land 11-3.8   dl11_4lo(2)=-12.0", bt11 - bt20, -12.0, ">="),
        ("Day Land 11-4.05  dl11_4lo(2)=-12.0", bt11 - bt21, -12.0, ">="),
        ("Nite Land 11-3.8  nl_11_4l(2)=-2.5", bt11 - bt20, -2.5, "<="),
        ("Nite Land 11-4.05 nl_11_4l(2)=-2.5", bt11 - bt21, -2.5, "<="),
        ("Desert 11-3.8 lo  lds11_4lo(2)=-23.0", bt11 - bt20, -23.0, ">="),
        ("Desert 11-3.8 hi  lds11_4hi(2)=-5.0", bt11 - bt20, -5.0, "<="),
        ("Nite Ocean 8.6-7.2 no86_73(2)=17.0", bt23 - bt22, 17.0, ">"),
        ("Nite Land 7.2-11  nl7_11s(2)=-10.0", bt22 - bt11, -10.0, "<="),
        ("Nite Snow 11-3.8  ns11_4lo(2)=0.6", bt11 - bt20, 0.6, "<="),
        ("Nite Snow 4-12    ns4_12hi(2)=5.5", bt20 - bt_data[25], 5.5, "<="),
    ]

    print(f"\n{'Test':<35} {'Operator':<5} {'Thr':>8} {'ClearPass':>10} {'CloudPass':>10}")
    print("-" * 85)

    for name, values, thr, op in tests:
        valid = ~np.isnan(values)
        if op == ">=":
            passes = values >= thr
        elif op == "<=":
            passes = values <= thr
        elif op == ">":
            passes = values > thr

        clear_pass = np.sum(passes[valid & is_clear]) / max(np.sum(valid & is_clear), 1) * 100
        cloud_pass = np.sum(passes[valid & is_cloudy]) / max(np.sum(valid & is_cloudy), 1) * 100

        print(f"{name:<35} {op:>5} {thr:>8.1f} {clear_pass:>9.1f}% {cloud_pass:>9.1f}%")

    # ============ Part 4: ch20 vs ch21 BT difference ============
    print(f"\n{'#'*80}")
    print(f"#  CH20(3.8um) vs CH21(4.05um) BT DIFFERENCE ANALYSIS")
    print(f"{'#'*80}")

    valid_pair = ~np.isnan(bt20) & ~np.isnan(bt21)
    diff_20_21 = bt20[valid_pair] - bt21[valid_pair]

    clear_pair = valid_pair & is_clear
    cloudy_pair = valid_pair & is_cloudy

    print(f"\n  BT3.8 - BT4.05 distribution:")
    for name, mask in [('All', valid_pair), ('Clear', clear_pair), ('Cloudy', cloudy_pair)]:
        d = diff_20_21[mask[valid_pair] if name != 'All' else np.ones(len(diff_20_21), dtype=bool)]
        if name != 'All':
            d = bt20[mask] - bt21[mask]
            d = d[~np.isnan(d)]
        print(f"  {name:<10}: mean={np.mean(d):.2f}K, std={np.std(d):.2f}K, "
              f"p5={np.percentile(d,5):.1f}, p50={np.percentile(d,50):.1f}, p95={np.percentile(d,95):.1f}")


if __name__ == '__main__':
    main()
