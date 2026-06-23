#!/bin/bash

# T639文件名
T639File=$1
# 输出路径
OutDir=$2


#调用样例   
# ./T639_to_bin.sh /export/home/fy4operation/DISK1/SourceData/NWP/20141119/ORG/gmf.639.2014111806042.grb2  /outdir
NameList=("PRES:surface"  "PRMSL:mean sea level"  "TMP:surface"  "HGT:surface"  "PRES:surface"  \
          "TMP:2 m above ground" "RH:2 m aboveground" "UGRD:10 m above ground" "VGRD:10 m above ground" \
          "PWAT"          "WEASD"         "TOZNE"  "TMP:250 mb"  \
          "TMP:0.1 mb" "TMP:0.2 mb" "TMP:0.5 mb" "TMP:1 mb" "TMP:1.5 mb" "TMP:2 mb" "TMP:3 mb" "TMP:4 mb" "TMP:5 mb" "TMP:7mb"\
          "TMP:10 mb" "TMP:20 mb" "TMP:30 mb" "TMP:50 mb" "TMP:70 mb" "TMP:100 mb" "TMP:150 mb"  "TMP:200 mb" "TMP:250 mb"  \
          "TMP:300 mb"  "TMP:350 mb"  "TMP:400 mb"  "TMP:450 mb"  "TMP:500 mb"  "TMP:550 mb"  "TMP:600 mb"  "TMP:650 mb"  \
          "TMP:700 mb"  "TMP:750 mb" "TMP:800 mb" "TMP:850 mb" "TMP:900 mb" "TMP:925 mb" "TMP:950 mb"  "TMP:975 mb" "TMP:1000 mb" \
          "HGT:0.1 mb"     "HGT:0.2 mb"    "HGT:0.5 mb"   "HGT:1 mb"     "HGT:1.5 mb"     "HGT:2 mb"   "HGT:3 mb"   "HGT:4 mb"   \
          "HGT:5 mb"   "HGT:7 mb"  "HGT:10 mb"     "HGT:20 mb"     "HGT:30 mb"   "HGT:50 mb"   "HGT:70 mb"   \
          "HGT:100 mb"  "HGT:150 mb"  "HGT:200 mb"  "HGT:250 mb"   "HGT:300 mb"  "HGT:350 mb"  "HGT:400 mb"  "HGT:450 mb" \
          "HGT:500 mb" "HGT:550 mb"  "HGT:600 mb"  "HGT:650 mb"  "HGT:700 mb"  "HGT:750 mb"  "HGT:800 mb"  "HGT:850 mb"  \
          "HGT:900 mb"  "HGT:925 mb"    "HGT:950 mb"    "HGT:975 mb"    "HGT:1000 mb"  \
          "O3MR:10 mb"  "O3MR:20 mb"    "O3MR:30 mb"    "O3MR:50 mb"    "O3MR:70 mb"    "O3MR:100 mb"   \
          "RH:0.1 mb"   "RH:0.2 mb"  "RH:0.5 mb"   "RH:1 mb"   "RH:1.5 mb"  "RH:2 mb"   "RH:3 mb"   "RH:4 mb"   "RH:5 mb"   \
          "RH:7 mb" "RH:10 mb" "RH:20 mb" "RH:30 mb" "RH:50 mb" "RH:70 mb" "RH:100 mb"  "RH:150 mb"  "RH:200 mb"  \
          "RH:250 mb"  "RH:300 mb"  "RH:350 mb"  "RH:400 mb"  "RH:450 mb" \
          "RH:500 mb"  "RH:550 mb"  "RH:600 mb"  "RH:650 mb"  "RH:700 mb"  "RH:750 mb"  "RH:800 mb"  "RH:850 mb"  "RH:900 mb" \
          "RH:925 mb"     "RH:950 mb"     "RH:975 mb"     "RH:1000 mb"  \
          "CLWMR:100 mb"  "CLWMR:150 mb"  "CLWMR:200 mb"  "CLWMR:250 mb"  "CLWMR:300 mb"  "CLWMR:350 mb"  "CLWMR:400 mb"  \
          "CLWMR:450 mb"  "CLWMR:500 mb"  "CLWMR:550 mb"  "CLWMR:600 mb"  "CLWMR:650 mb"  "CLWMR:700 mb"  "CLWMR:750 mb"  \
          "CLWMR:800 mb"  "CLWMR:850 mb"  "CLWMR:900 mb"  "CLWMR:925 mb"  "CLWMR:950 mb"  "CLWMR:975 mb"  "CLWMR:1000 mb")
  
#获取T639文件名上的时间
start_time=`echo $T639File | grep -Eo [0-9]{13}`
ymd=${start_time:0:8}
hh=${start_time:8:2}
#预报时间,小时计数
hours=`echo ${start_time:10:3} | bc`


#转换预报标准格式的预报时间
new_time=`date +"%Y%m%d_%H_00" --date="$ymd $hh:00:00 $hours hours"`

BinFile=$OutDir/T639_$new_time


rm -rf $BinFile
for name in "${NameList[@]}"
do
    #echo $name
    wgrib2 -set local_table 1 -s $T639File | grep "$name" |wgrib2 -i $T639File -no_header -append -bin $BinFile >/dev/null
done
#wgrib2 -set local_table 1 -s $T639File | grep "TMP:2 m above ground" |wgrib2 -i $T639File -no_header -append -text $BinFile
#wgrib2 -set local_table 1 -s $T639File | grep "$name" |wgrib2 -i $T639File -no_header -append -bin $BinFile
