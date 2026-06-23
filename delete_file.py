#!/usr/bin/env python3
# coding:utf8
#--------------------
__author__ = 'Min Min'
#--------------------

def del_file(path):
	ls = os.listdir(path)
	for i in ls:
		c_path = os.path.join(path, i)
		if os.path.isdir(c_path):
			del_file(c_path)
		else:
			os.remove(c_path)
			

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

print ('  L1_path   = '+L1_path)

#        start  end
year  = [2018, 2018]
month = [   5,    5]
day   = [ 22,   28]
hour_lmt  = [   0,   24]

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
	directory_name = L1_path+date_fn+'/'
	
	print ('directory_name = ',directory_name)
	
	print ('delete Cloud file!')
	if os.path.isdir(directory_name+'Cloud'):
		del_file(directory_name+'Cloud')
		
	print ('delete Aviation file!')
	if os.path.isdir(directory_name+'Aviation'):
		del_file(directory_name+'Aviation')
	
	print ('delete Aerosol file!')
	if os.path.isdir(directory_name+'Aerosol'):
		del_file(directory_name+'Aerosol')
	
	print ('delete Surface file!')
	if os.path.isdir(directory_name+'Surface'):
		del_file(directory_name+'Surface')
	
	print ('delete Weather file!')
	if os.path.isdir(directory_name+'Weather'):
		del_file(directory_name+'Weather')
	
	print ('delete Radiation Figure file!')
	if os.path.isdir(directory_name+'Radiation/Figure'):
		del_file(directory_name+'Radiation/Figure')	
	
	
print ('day cycle over!!!!')
	
