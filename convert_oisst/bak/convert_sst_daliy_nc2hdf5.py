#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2020/6/3 16:30
# @Author  : Min Min
# @File    : Predict.py
import sys
import os
import numpy as np
import h5py
import glob
import netCDF4

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
	
	
	

nc_path   = '/wind1/home/minm/Data_fy4test/oisst/ORG/'
hdf5_path = '/wind1/home/minm/Data_fy4test/oisst/'
#        start  end
date_name = sys.argv[1]


print(nc_path)
#oisst-avhrr-v02r01.20200330.nc ; sst.day.mean.30000330.hdf5
nc_fname = nc_path + 'oisst-avhrr-v02r01.'+date_name+'.nc'
hdf5_fname = hdf5_path + 'sst.day.mean.'+date_name+'.hdf5'

print ('read and deal with nc sst data')
f = netCDF4.Dataset(nc_fname,'r')
sst0 = f.variables['sst']   # skt
sst0 = np.array(sst0)
print (sst0.shape)

sst = np.zeros([sst0.shape[2],sst0.shape[3]],dtype=np.float32)
sst[:,:] = sst0[0,0,:,:]

loc = np.where(sst < -1)
sst[loc] = -9999.0
print (sst.shape)

f.close()

print ('write hdf5 sst data')
f1 = h5py.File(hdf5_fname,'w')
a1 = f1.create_dataset('sst', data = np.float32(sst))
f1.close

