#!/bin/bash

# NCEP文件名
NCEPFile=$1
# 输出路径
OutFile=$2


#调用样例   
# ./NCEP_to_bin.sh /export/home/fy4operation/DISK1/SourceData/NWP/20150201/ORG/fnl_20150201_00_00.grib2  outfile
NameList=("PRES:surface"  "PRMSL:mean"  "TMP:surface"  "HGT:surface"  "PRES:surface"  \
          "TMP:0.995 sigma level" "RH:0.995 sigma level" "UGRD:0.995 sigma level" "VGRD:0.995 sigma level" \
          "PWAT"          "WEASD:surface"         "TOZNE"  "TMP:tropopause"  \
          "TMP:1 mb" "TMP:2 mb" "TMP:3 mb" "TMP:5 mb" "TMP:7 mb"  \
          "TMP:10 mb" "TMP:20 mb" "TMP:30 mb" "TMP:50 mb" "TMP:70 mb" "TMP:100 mb" "TMP:150 mb"  "TMP:200 mb" "TMP:250 mb"  \
          "TMP:300 mb"  "TMP:350 mb"  "TMP:400 mb"  "TMP:450 mb"  "TMP:500 mb"  "TMP:550 mb"  "TMP:600 mb"  "TMP:650 mb"  \
          "TMP:700 mb"  "TMP:750 mb"  "TMP:800 mb"  "TMP:850 mb"  "TMP:900 mb"  "TMP:925 mb"  "TMP:950 mb"  "TMP:975 mb" "TMP:1000 mb" \
          "HGT:1 mb"     "HGT:2 mb"     "HGT:3 mb"   "HGT:5 mb"   "HGT:7 mb"   \
          "HGT:10 mb"     "HGT:20 mb"     "HGT:30 mb"   "HGT:50 mb"   "HGT:70 mb"   \
          "HGT:100 mb"  "HGT:150 mb"  "HGT:200 mb"  "HGT:250 mb"   "HGT:300 mb"  "HGT:350 mb"  "HGT:400 mb"  "HGT:450 mb" \
          "HGT:500 mb"  "HGT:550 mb"  "HGT:600 mb"  "HGT:650 mb"   "HGT:700 mb"  "HGT:750 mb"  "HGT:800 mb"  "HGT:850 mb"  \
          "HGT:900 mb"  "HGT:925 mb"    "HGT:950 mb"    "HGT:975 mb"    "HGT:1000 mb"  \
          "O3MR:1 mb"   "O3MR:2 mb"     "O3MR:3 mb"    "O3MR:5 mb"    "O3MR:7 mb"  \
          "O3MR:10 mb"  "O3MR:20 mb"    "O3MR:30 mb"    "O3MR:50 mb"    "O3MR:70 mb"    "O3MR:100 mb"   \
          "O3MR:150 mb"  "O3MR:200 mb"  "O3MR:250 mb"  "O3MR:300 mb"  "O3MR:350 mb"  "O3MR:400 mb"  \
          "RH:1 mb"   "RH:2 mb"   "RH:3 mb"   "RH:5 mb"   "RH:7 mb" \
          "RH:10 mb"  "RH:20 mb"  "RH:30 mb"  "RH:50 mb"  "RH:70 mb" \
          "RH:100 mb"  "RH:150 mb"  "RH:200 mb"  "RH:250 mb"  "RH:300 mb"  "RH:350 mb"  "RH:400 mb"  "RH:450 mb" \
          "RH:500 mb"  "RH:550 mb"  "RH:600 mb"  "RH:650 mb"  "RH:700 mb"  "RH:750 mb"  "RH:800 mb"  "RH:850 mb"  "RH:900 mb" \
          "RH:925 mb"  "RH:950 mb"  "RH:975 mb"  "RH:1000 mb"  \
          "CLWMR:100 mb"  "CLWMR:150 mb"  "CLWMR:200 mb"  "CLWMR:250 mb"  "CLWMR:300 mb"  "CLWMR:350 mb"  "CLWMR:400 mb"  \
          "CLWMR:450 mb"  "CLWMR:500 mb"  "CLWMR:550 mb"  "CLWMR:600 mb"  "CLWMR:650 mb"  "CLWMR:700 mb"  "CLWMR:750 mb"  \
          "CLWMR:800 mb"  "CLWMR:850 mb"  "CLWMR:900 mb"  "CLWMR:925 mb"  "CLWMR:950 mb"  "CLWMR:975 mb"  "CLWMR:1000 mb"  \
          "UGRD:1 mb"     "UGRD:2 mb"     "UGRD:3 mb"   "UGRD:5 mb"   "UGRD:7 mb"   \
          "UGRD:10 mb"     "UGRD:20 mb"     "UGRD:30 mb"   "UGRD:50 mb"   "UGRD:70 mb"   \
          "UGRD:100 mb"  "UGRD:150 mb"  "UGRD:200 mb"  "UGRD:250 mb"   "UGRD:300 mb"  "UGRD:350 mb"  "UGRD:400 mb"  "UGRD:450 mb" \
          "UGRD:500 mb"  "UGRD:550 mb"  "UGRD:600 mb"  "UGRD:650 mb"   "UGRD:700 mb"  "UGRD:750 mb"  "UGRD:800 mb"  "UGRD:850 mb"  \
          "UGRD:900 mb"  "UGRD:925 mb"    "UGRD:950 mb"    "UGRD:975 mb"    "UGRD:1000 mb"  \
          "VGRD:1 mb"     "VGRD:2 mb"     "VGRD:3 mb"   "VGRD:5 mb"   "VGRD:7 mb"   \
          "VGRD:10 mb"     "VGRD:20 mb"     "VGRD:30 mb"   "VGRD:50 mb"   "VGRD:70 mb"   \
          "VGRD:100 mb"  "VGRD:150 mb"  "VGRD:200 mb"  "VGRD:250 mb"   "VGRD:300 mb"  "VGRD:350 mb"  "VGRD:400 mb"  "VGRD:450 mb" \
          "VGRD:500 mb"  "VGRD:550 mb"  "VGRD:600 mb"  "VGRD:650 mb"   "VGRD:700 mb"  "VGRD:750 mb"  "VGRD:800 mb"  "VGRD:850 mb"  \
          "VGRD:900 mb"  "VGRD:925 mb"   "VGRD:950 mb"    "VGRD:975 mb"    "VGRD:1000 mb"  "UGRD:10 m above ground"  "VGRD:10 m above ground" )

#NameList2=("PRES:surface"  "PRMSL:mean")  

rm -rf $OutFile
for name in "${NameList[@]}"
#for name in "${NameList2[@]}"
do
    #echo $name
    /opt/software/grib2/wgrib2/wgrib2 -set local_table 1 -s $NCEPFile | grep "$name" |/opt/software/grib2/wgrib2/wgrib2 -i $NCEPFile -no_header -append -bin $OutFile >/dev/null
done
#/opt/software/grib2/wgrib2/wgrib2 -set local_table 1 -s $NCEPFile | grep "TMP:2 m above ground" |/opt/software/grib2/wgrib2/wgrib2 -i $NCEPFile -no_header -append -text $OutFile
#/opt/software/grib2/wgrib2/wgrib2 -set local_table 1 -s $NCEPFile | grep "$name" |/opt/software/grib2/wgrib2/wgrib2 -i $NCEPFile -no_header -append -bin $OutFile
