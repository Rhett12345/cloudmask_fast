#!/bin/bash
#===============================================================
# Name:        main_run_cloud.sh
# author:      wangpeng
# Version:     0.0.1
# Date:        2017-03-16
# Description: 运行云检测程序
# Input:       
#     1、日期(YYYYMMDD)
#     2、时分(HHMM)
# OutPut:      
#     1、转换后的海面温度产品(nc年产品转换到hdf5的日产品)
#     2、云检测产品
#===============================================================

stime=$1
etime=$2


# 进入主程序根位置
ROOT_PATH=/home/huxq/fy3cloud/fy3_mersi2_cloud/retrieval
cd $ROOT_PATH

#海面温度转换进程
OISST_CONVERT_EXE=$ROOT_PATH/convert_oisst/oisst_daily_nc2hdf5
#海面温度产品路径
OISST_IPATH=/PUBLICDATAFY3PGS/FTPDATA/SAT/Data_fy3cloud/oisst/daily/ORG/
OISST_OPATH=/PUBLICDATAFY3PGS/FTPDATA/SAT/Data_fy3cloud/oisst/daily/

# 开始转换
# sst.day.mean.20050101.hdf5
while :
do
    if [[ $stime -gt $etime ]]; then
        break;
    fi 
    echo "start convert oisst year -> day $stime"
    otime=$(date -d "$stime 2day ago"  +%Y%m%d) 
    oisst_ofile=$OISST_OPATH/sst.day.mean.${otime}.hdf5
    if [ -f $oisst_ofile ];then
       echo "$oisst_ofile exist" 
       stime=$(date -d "$stime 1day"  +%Y%m%d)
       continue
    fi
    
    # 执行转换命令
    $OISST_CONVERT_EXE $OISST_IPATH $OISST_OPATH $stime > /dev/null
    if [ $? -eq 0 ];then
        echo "convert oisst Success"
    else
	echo "conver oisst Failed"
    fi 
    stime=$(date -d "$stime 1day"  +%Y%m%d)
done
