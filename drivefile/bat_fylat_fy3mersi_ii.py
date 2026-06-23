#!/usr/bin/env python3
# coding:utf8
#--------------------
__author__ = 'Min Min'
#--------------------


def julianday(IY, IM, ID, IH, MIT): 


# 2. begin program
	if IM <= 2:  # january & february
		IY1 = IY-1
		IM1 = IM+12
		JD = float(int( 365.25*(IY1 + 4716.0)) + int( 30.6001*( IM1 + 1.0)) + 2.0 - int( IY1/100.0 ) + int( int( IY1/100.0 )/4.0 ) + ID - 1524.5) + float((IH + MIT/60.+0./3600.)/24.)
     
	else:
  		JD = float( int( 365.25*(IY + 4716.0)) + int( 30.6001*( IM + 1.0)) + 2.0 - int( IY/100.0 ) + int( int( IY/100.0 )/4.0 ) + ID - 1524.5) + float((IH + MIT/60.+0./3600.)/24.)
     
	return JD



def julian_to_date(JD):

	I = int(JD + 0.5)
	Fr = abs( I - ( JD + 0.5) )	
	
	if I >= 2299160.0:
		A  = int( ( I- 1867216.25 ) / 36524.25 )
		a4 = int( A / 4 )
		B  = I + 1. + float(A - a4)
	else:
		B = I

	C = B + 1524.
	D = int( ( C - 122.1 ) / 365.25 )
	E = int( 365.25 * D )
	G = int( ( C - E ) / 30.6001 )
	day = int( C - E + Fr - int( 30.6001 * G ) )

	if G <= 13.5:
		month = int(G - 1)
	else:
		month = int(G - 13)


	if month > 2.5:
		year = int(D - 4716)
	else:
		year = int(D - 4715)


	hour = int( Fr * 24. )
	mint = int( abs( hour -( Fr * 24. ) ) * 60. )
	return (year, month, day)


def int_to_str(a):
	
	if a<10:
		c = '0'+str(a)
	else:
		c = str(a)

	return c

def find_nwp_name(hour):
	
	nwp_hour = [0,3,6,9,12,15,18,21,24]
	for i in range(0,8):
		if (hour>=nwp_hour[i] and hour<nwp_hour[i+1]):
			n1 = int_to_str(nwp_hour[i]+18)
			n2 = int_to_str(nwp_hour[i+1]+18)

	return (n1,n2)



import sys
sys.dont_write_bytecode = True
import os 
import glob
import time
#import numpy as np 
#reload(sys)
#sys.setdefaultencoding('utf8')


#if __name__ == '__main__':
# ---------- MAIN ----------
# step 0: set path, time and sensor id
print ('   ')
print ('  --- step 0: input data saving path and result output path')
L1_path   = '/ARSDISK/DATA/FY3/FY3D/MERSI/L1/1000M/'  #'/PUBLICDATAFY3PGS/FTPDATA/SAT/Data_fy3cloud/FY3D_MERSI_II_test/L1/'
GEO_path  = '/ARSDISK/DATA/FY3/FY3D/MERSI/L1/GEO1K/' #'/PUBLICDATAFY3PGS/FTPDATA/SAT/Data_fy3cloud/FY3D_MERSI_II_test/GEO1K/'
nwp_path  = '/PUBLICDATAFY3PGS/FTPDATA/SAT/Data_fy3cloud/FY3D_MERSI_II_test/nwp/'
sst_path  = '/PUBLICDATAFY3PGS/FTPDATA/SAT/Data_fy3cloud/oisst/daily/'
L2_path   = '/PUBLICDATAFY3PGS/FTPDATA/SAT/Data_fy3cloud/FY3D_MERSI_II_test/L2/'
code_path = '/home/huxq/fy3cloud/fy3mersi_ii/retrieval/'
print ('  L1_path   = '+L1_path)
print ('  GEO_path  = '+GEO_path)
print ('  nwp_path  = '+nwp_path)
print ('  sst_path  = '+sst_path)
print ('  L2_path   = '+L2_path)
print ('  code_path = '+code_path)
sensor_id = 21  # 41=h8; 31=fy4a 
nwp_id    = 4  # 41=h8; 31=fy4a 
#  !  fylat_nwp_opt	          
#  !   1 = ncep reanalysis 1*1 (grib1)
#  !   2 = gfs1p00 1*1 (grib2)
#  !   3 = T639 0.125*0.125 (grib2)
#  !   4 = ncep reanalysis 1*1 (grib2) 
#  !   5 = gfs0p50 0.5*0.5 (grib2) 
#  !   6 = EC prediction 1*1 (grib2) 
#  !   7 = grapes 0.25*0.25 (grib2)

#        start  end
year  = [2019, 2019]
month = [   4,  4]
day   = [ 19,   19]
hour_lmt  = [   0,   24]
mint_lmt  = [  0,  59]

# step 1: calculate year month day
print('   ')
print('  --- step 1: calculate year month day')
jd1  = julianday(year[0],month[0],day[0],0,0)
jd2  = julianday(year[1],month[1],day[1],0,0)
nday = int(jd2-jd1+1)


for iday in range(0,nday):

	print ('----------------------------------')
	print ('day number:', iday+1)
	jd = jd1 + iday
	#convert date
	[y1,m1,d1] = julian_to_date(jd)
	print ('Date = ', y1, m1, d1)
	yn = int_to_str(y1)
	mn = int_to_str(m1)
	dn = int_to_str(d1)
	
	date_fn = yn+mn+dn
	year_fn = yn
	#directory_name = L1_path+date_fn+'/'
	directory_name = L1_path+year_fn+'/'+date_fn+'/'
	print ('directory_name = ',directory_name)
	if os.path.isdir(directory_name):
	
		if sensor_id == 21:  # fy3d mersi_ii  ;FY3D_MERSI_GBAL_L1_20180220_2355_1000M_MS.HDF
			file_list = glob.glob(directory_name+'FY3D_MERSI_GBAL_L1_'+date_fn+'_*_1000M_MS.HDF')
			file_list = sorted(file_list)
		if sensor_id == 11:  # fy3d mersi_ii using modis
			file_list = glob.glob(directory_name+'FY3D_MERSI_GBAL_L1_'+date_fn+'_*_1000M_MS.HDF')
			file_list = sorted(file_list)
		
		#print(file_list)
			
		i=1
		for filename in file_list:
		
			retr_id = 0
			print ('   ')
			print ('file number:', i,'/',len(file_list))
			temp_fn = os.path.split(filename)
			fy3_mersi_L1b_data = filename
			fy3d_L1_fname = temp_fn[1]
			hm = fy3d_L1_fname[28:32]
			hour = int(fy3d_L1_fname[28:30])
			mint = int(fy3d_L1_fname[30:32])
			print ('fy3d L1 file  = ',fy3d_L1_fname)
			print ('hour minute = ', hm, hour, mint)
			# find geo data
			print (GEO_path+year_fn+'/'+date_fn+'/'+'FY3D_MERSI_GBAL_L1_'+date_fn+'_'+hm+'_GEO1K_MS.HDF')
			file_geo = glob.glob(GEO_path+year_fn+'/'+date_fn+'/'+'FY3D_MERSI_GBAL_L1_'+date_fn+'_'+hm+'_GEO1K_MS.HDF')
			print (len(file_geo))
			if (len(file_geo)>0 and (hour*60.+mint)>=(hour_lmt[0]*60.+mint_lmt[0]) and (hour*60.+mint)<=(hour_lmt[1]*60.+mint_lmt[1])):
				retr_id = 1
				print ('fy3 geo file = ', file_geo[0])
				fy3_mersi_GEO_data = file_geo[0]
				print ('FYLAT Retrieval is OK = ')
				# FY3D_MERSI_ORBT_L2_CLM_MLT_NUL_20180721_2355_1000M_MS.HDF
				fy3_mersi_CLM_data  = L2_path+date_fn+'/'+'FY3D_MERSI_ORBT_L2_CLM_MLT_NUL_'+date_fn+'_'+hm+'_1000M_MS.HDF'
				fy3_mersi_CLA_data  = L2_path+date_fn+'/'+'FY3D_MERSI_ORBT_L2_CLA_MLT_NUL_'+date_fn+'_'+hm+'_5000M_MS.HDF'
				fy3_mersi_CLP_data  = L2_path+date_fn+'/'+'FY3D_MERSI_ORBT_L2_CLP_MLT_NUL_'+date_fn+'_'+hm+'_1000M_MS.HDF'
				fy3_mersi_CTP_data  = L2_path+date_fn+'/'+'FY3D_MERSI_ORBT_L2_CTP_MLT_NUL_'+date_fn+'_'+hm+'_1000M_MS.HDF'
				fy3_mersi_COT_data  = L2_path+date_fn+'/'+'FY3D_MERSI_ORBT_L2_COT_MLT_NUL_'+date_fn+'_'+hm+'_1000M_MS.HDF'
				fy3_mersi_CON_data  = L2_path+date_fn+'/'+'FY3D_MERSI_ORBT_L2_CON_MLT_NUL_'+date_fn+'_'+hm+'_1000M_MS.HDF'
				fy3_intermediate    = L2_path+date_fn+'/'+'FY3D_MERSI_ORBT_L2_XXX_MLT_NUL_'+date_fn+'_'+hm+'_INTERMED.HDF'
				
				# oisst date   ;sst.day.mean.20180201.hdf5
				[y2,m2,d2] = julian_to_date(jd-3)
				print ('Date oisst = ', y2, m2, d2)
				yn2 = int_to_str(y2)
				mn2 = int_to_str(m2)
				dn2 = int_to_str(d2)
				oisst_data =sst_path+'sst.day.mean.'+yn2+mn2+dn2+'.hdf5'
				
				# nwp   ;gfs.t06z.pgrb2.0p50.f021
				[nwp1,nwp2] = find_nwp_name(hour)
				print ('nwp time = ',nwp1,' ',nwp2)
				nwp_grib_data1 = nwp_path+date_fn+'/ORG/gfs.t06z.pgrb2.0p50.f0'+nwp1
				nwp_grib_data2 = nwp_path+date_fn+'/ORG/gfs.t06z.pgrb2.0p50.f0'+nwp2
				
			else:
				print ('ERROR: No GEO data!')
				continue
			
			i = i+1
			#++++++++++
			if (retr_id == 1):
				print ('Write fy3d config control file')
				f = open('temp_fy3d_config.nml','w+')
				f.write('&config'+'\n') 
				f.write('  fylat_sensor_id     = '+str(sensor_id)+',\n')
				f.write('  code_root_path      = "'+code_path+'",\n')
				f.write('  L1b_data_path       = "'+L1_path+year_fn+'/'+date_fn+'/'+'",\n')
				f.write('  nwp_data_path       = "'+nwp_path+date_fn+'/'+'",\n')
				f.write('  oisst_data_path     = "'+sst_path+'",\n')
				f.write('  fy3_mersi_GEO_data  = "'+fy3_mersi_GEO_data+'",\n')
				f.write('  fy3_mersi_L1b_data  = "'+fy3_mersi_L1b_data+'",\n')
				f.write('  fy3_mersi_CLM_data  = "'+fy3_mersi_CLM_data+'",\n')
				f.write('  fy3_mersi_CLA_data  = "'+fy3_mersi_CLA_data+'",\n')
				f.write('  fy3_mersi_CLP_data  = "'+fy3_mersi_CLP_data+'",\n')
				f.write('  fy3_mersi_CTP_data  = "'+fy3_mersi_CTP_data+'",\n')
				f.write('  fy3_mersi_COT_data  = "'+fy3_mersi_COT_data+'",\n')
				f.write('  fy3_mersi_CON_data  = "'+fy3_mersi_CON_data+'",\n')
				f.write('  fy3_intermediate    = "'+fy3_intermediate+'",\n')
				f.write('  fylat_nwp_opt       = '+str(nwp_id)+',\n')
				f.write('  fylat_rtm_opt       = 1'+',\n')
				f.write('  nwp_grib_data1      = "'+nwp_grib_data1+'",\n')
				f.write('  nwp_grib_data2      = "'+nwp_grib_data2+'",\n')
				f.write('  oisst_data          = "'+oisst_data+'",\n')
				f.write('  cloudmask_id        = 1'+',\n')
				f.write('  cloudamount_id      = 0'+',\n')
				f.write('  cloudphase_id       = 0'+',\n')
				f.write('  cloudtopz_id        = 0'+',\n')
				f.write('  cloudtau_day_id     = 0'+',\n')
				f.write('  cloudtau_night_id   = 0'+',\n')
				f.write('  cloudtypeII_id      = 0'+',\n')
				f.write('  write_inter_id      = 0/'+'\n')
				f.close()
			
				print ('Drive fy3 mersi_ii FYLAT')
				os.chdir(code_path)
				os.system('./fylat_FY3_MERSI_II_PGS temp_fy3d_config.nml')
				os.chdir(code_path)
	else:
		print ('find next directory!!!')

else:
    print ('day cycle over!!!!')
	
