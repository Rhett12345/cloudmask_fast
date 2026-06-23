#!/usr/bin/env python3
# coding:utf8
#--------------------
__author__ = 'Min Min'
#--------------------


import sys
#sys.dont_write_bytecode = True
import os 
import glob
import h5py
#import time
#import multiprocessing 
#from multiprocessing import Pool
import numpy as np 
import netCDF4
#reload(sys)
#sys.setdefaultencoding('utf8')


#if __name__ == '__main__':
# ---------- MAIN ----------
# step 0: set path, time and sensor id
print ('   ')
print ('  --- input data saving path and result output path')
#example for using this code: 
#     python3 convert_sst_nc2hdf5_code.py /500T/AVHRR_sst/2022/oisst-avhrr-v02r01.20220909_preliminary.nc /500T/AVHRR_sst/hdf5/

argv_array1 = sys.argv[1]
argv_array2 = sys.argv[2]

nc_file = argv_array1
hdf5_path = argv_array2
fff = os.path.basename(nc_file)
date_name = fff[19:27]
print (' date name = ', date_name)

hdf5_fname = hdf5_path + 'sst.day.mean.'+date_name+'.hdf5'

code_path = sys.path[0]+'/'

print ('read and deal with nc sst data')
f = netCDF4.Dataset(nc_file,'r')
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

print ('  --- program is over !!!')
	
