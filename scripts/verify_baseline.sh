#!/bin/bash
# Verify current code produces identical results to v0.2.1 baseline
set -e
cd "$(dirname "$0")/.."

OUTDIR=/tmp/fylat_verify_$$
mkdir -p "$OUTDIR"
NML="$OUTDIR/test.nml"

cat > "$NML" << 'NML'
&config
  fylat_sensor_id     = 21,
  code_root_path      = "CODEPATH",
  L1b_data_path       = "/data/Data_yuq/mersi/20220803/",
  nwp_data_path       = "/data/nwp/20220803/",
  oisst_data_path     = "/data/Data_minmin/oisst/",
  fy3_mersi_GEO_data  = "/data/Data_yuq/mersi/20220803/FY3D_MERSI_GBAL_L1_20220803_0740_GEO1K_MS.HDF",
  fy3_mersi_L1b_data  = "/data/Data_yuq/mersi/20220803/FY3D_MERSI_GBAL_L1_20220803_0740_1000M_MS.HDF",
  fy3_mersi_CLM_data  = "OUTDIR/v021_CLM.HDF",
  fy3_mersi_CLA_data  = "OUTDIR/v021_TMP.HDF", fy3_mersi_CLP_data  = "OUTDIR/v021_TMP.HDF",
  fy3_mersi_CTP_data  = "OUTDIR/v021_TMP.HDF", fy3_mersi_COT_data  = "OUTDIR/v021_TMP.HDF",
  fy3_mersi_CON_data  = "OUTDIR/v021_TMP.HDF", fy3_mersi_SST_data  = "OUTDIR/v021_TMP.HDF",
  fy3_intermediate    = "OUTDIR/v021_TMP.HDF",
  fylat_nwp_opt       = 10, fylat_rtm_opt = 1,
  nwp_grib_data1      = "/data/nwp/20220803/ORG/gfs.t06z.pgrb2.0p25.f018",
  nwp_grib_data2      = "/data/nwp/20220803/ORG/gfs.t06z.pgrb2.0p25.f021",
  oisst_data          = "/data/Data_minmin/oisst/sst.day.mean.20200401.hdf5",
  cloudmask_id        = 1, cloudamount_id = 0, cloudphase_id = 0, cloudtopz_id = 0,
  cloudtau_day_id     = 0, cloudtau_night_id = 0, cloudtypeII_id = 0,
  surface_sst_id      = 0, write_inter_id = 0/
NML

CURDIR=$(pwd)
sed -i "s|CODEPATH|$CURDIR/|g" "$NML"
sed -i "s|OUTDIR|$OUTDIR|g" "$NML"

echo "=== Building v0.2.1 baseline ==="
git stash -q 2>/dev/null
CUR_BRANCH=$(git branch --show-current)
git checkout v0.2.1 -q
cd src && make -f fylat_makefile_cldmask -j$(nproc) 2>&1 | tail -1 && cd ..
cp VIS_Cal_Coeff_business.xcfg VIS_Cal_Coeff.xcfg 2>/dev/null; rm -f cal_mode.txt
echo "Running v0.2.1..."
./fylat_FY3_MERSI_II_PGS "$NML" 2>&1 | grep "FYLAT Over"

echo "=== Building current ==="
git checkout "$CUR_BRANCH" -q
git stash pop -q 2>/dev/null
cd src && make -f fylat_makefile_cldmask -j$(nproc) 2>&1 | tail -1 && cd ..
cp VIS_Cal_Coeff_business.xcfg VIS_Cal_Coeff.xcfg 2>/dev/null; rm -f cal_mode.txt
sed "s|v021_CLM|v09_CLM|g" "$NML" > "${NML}_curr"
echo "Running current..."
./fylat_FY3_MERSI_II_PGS "${NML}_curr" 2>&1 | grep "FYLAT Over"

echo ""
echo "=== Comparing datasets ==="
python3 -c "
import h5py, numpy as np
f1 = h5py.File('$OUTDIR/v021_CLM.HDF','r')
f2 = h5py.File('$OUTDIR/v09_CLM.HDF','r')
ok = True
for k in sorted(f1.keys()):
    d1, d2 = f1[k][:], f2[k][:]
    ndiff = int(np.sum(d1 != d2))
    if ndiff > 0:
        print(f'  FAIL: {k}: {ndiff}/{d1.size} diffs')
        ok = False
f1.close(); f2.close()
if ok:
    print('  PASS: ALL datasets identical to v0.2.1 baseline')
else:
    print('  FAIL: Differences detected!')
    exit(1)
" 2>&1

rm -rf "$OUTDIR"
echo "Verification complete."
