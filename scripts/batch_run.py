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

def find_nwp_name2(yn, date_fn, path, hour): 

	nwp_hour = [0,3,6,9,12,15,18,21,24]
	for i in range(0,8):
		if (hour>=nwp_hour[i] and hour<nwp_hour[i+1]):
			n1 = nwp_hour[i]
			n2 = nwp_hour[i+1]
			
	print ('convert t639 nwp grib2 to binary data',n1,n2)
	file_list = glob.glob(path+date_fn+'/ORG/gmf.639.'+yn+'*grb2')
	file_list = sorted(file_list)
	y = int(yn)
	m = int(date_fn[4:6])
	d = int(date_fn[6:8])
	#print ('y , m, d',y , m, d)
	
	for filename in file_list:
		temp_fn = os.path.split(filename)
		print (filename)
		fname = temp_fn[1]
		path_nwp0 = temp_fn[0]
		len1 = len(path_nwp0)
		path_nwp = path_nwp0[0:len1-4]
		print ('path_nwp = ',path_nwp)
		#gmf.639.2018093012024.grb2
		tnn = fname[8:16]
		ps = int(fname[16:18])
		tl = int(fname[18:21])
		print ('tnn ps tl = ', tnn, ps, tl)
		h1 = ps + tl - 24
		
		if (h1 == n1) :
			nwp1 = fname
		if (h1 == n2) :
			nwp2 = fname
	print ('find t639 nwp data over!!!',nwp1, nwp2)
	
	return (nwp1,nwp2)


def find_nwp_name3(yn, date_fn, path, hour): 

	nwp_hour = [0,3,6,9,12,15,18,21,24]
	for i in range(0,8):
		if (hour>=nwp_hour[i] and hour<nwp_hour[i+1]):
			n1 = nwp_hour[i]
			n2 = nwp_hour[i+1]
			
	print ('convert grapes nwp grib2 to binary data',n1,n2)
	file_list = glob.glob(path+date_fn+'/ORG/gmf.gra.'+yn+'*grb2')
	file_list = sorted(file_list)
	y = int(yn)
	m = int(date_fn[4:6])
	d = int(date_fn[6:8])
	#print ('y , m, d',y , m, d)
	
	for filename in file_list:
		temp_fn = os.path.split(filename)
		print (filename)
		fname = temp_fn[1]
		path_nwp0 = temp_fn[0]
		len1 = len(path_nwp0)
		path_nwp = path_nwp0[0:len1-4]
		print ('path_nwp = ',path_nwp)
		#gmf.639.2018093012024.grb2
		tnn = fname[8:16]
		ps = int(fname[16:18])
		tl = int(fname[18:21])
		print ('tnn ps tl = ', tnn, ps, tl)
		h1 = ps + tl - 24
		
		if (h1 == n1) :
			nwp1 = fname
		if (h1 == n2) :
			nwp2 = fname
	print ('find grapes nwp data over!!!',nwp1, nwp2)
	
	return (nwp1,nwp2)


def find_nwp_name4(hour):
	
	nwp_hour = [0,6,12,18,24]
	for i in range(0,3):
		if (hour>=nwp_hour[i] and hour<nwp_hour[i+1]):
			n1 = int_to_str(nwp_hour[i])
			n2 = int_to_str(nwp_hour[i+1])

	return (n1,n2)


def driver_fylat(filename, year_fn, date_fn, date_fn_p1, jd, sensor_id, nwp_id, L1_path, GEO_path, L2_path, code_path, sst_path, nwp_path, hour_lmt ):
	
	retr_id = 0
	print ('   ')
#	print ('file number:', i,'/',len(file_list))
	temp_fn = os.path.split(filename)
	fy3_mersi_L1b_data = filename
	fy3d_L1_fname = temp_fn[1]
	hm = fy3d_L1_fname[28:32]
	hour = int(fy3d_L1_fname[28:30])
	mint = int(fy3d_L1_fname[30:32])
	print ('fy3d L1 file  = ',fy3d_L1_fname)
	print ('hour minute = ', hm, hour, mint)
	# find geo data
	file_geo = glob.glob(GEO_path+year_fn+'/'+date_fn+'/'+'FY3D_MERSI_GBAL_L1_'+date_fn+'_'+hm+'_GEO1K_MS.HDF')
	#file_geo = glob.glob(GEO_path+date_fn+'/'+'FY3D_MERSI_GBAL_L1_'+date_fn+'_'+hm+'_GEO1K_MS.HDF')

	if os.path.isdir(L2_path+date_fn+'/'):
		print ('output dir is ok!')
	else:
		os.system('mkdir '+L2_path+date_fn)
	
	if (len(file_geo)>0 and hour>=hour_lmt[0] and hour<=hour_lmt[1]):
		retr_id = 1
		print ('fy3 geo file = ', file_geo[0])
		fy3_mersi_GEO_data = file_geo[0]
		print ('FYLAT Retrieval is OK ')
		fy3_mersi_CLM_data  = L2_path+date_fn+'/'+'FY3D_MERSI_ORBT_L2_CLM_MLT_NUL_'+date_fn+'_'+hm+'_1000M_MS.HDF'
		fy3_mersi_CLA_data  = L2_path+date_fn+'/'+'FY3D_MERSI_ORBT_L2_CLA_MLT_NUL_'+date_fn+'_'+hm+'_5000M_MS.HDF'
		fy3_mersi_CLP_data  = L2_path+date_fn+'/'+'FY3D_MERSI_ORBT_L2_CLP_MLT_NUL_'+date_fn+'_'+hm+'_1000M_MS.HDF'
		fy3_mersi_CTP_data  = L2_path+date_fn+'/'+'FY3D_MERSI_ORBT_L2_CTP_MLT_NUL_'+date_fn+'_'+hm+'_1000M_MS.HDF'
		fy3_mersi_COT_data  = L2_path+date_fn+'/'+'FY3D_MERSI_ORBT_L2_COT_MLT_NUL_'+date_fn+'_'+hm+'_1000M_MS.HDF'
		fy3_mersi_CON_data  = L2_path+date_fn+'/'+'FY3D_MERSI_ORBT_L2_CON_MLT_NUL_'+date_fn+'_'+hm+'_1000M_MS.HDF'
		fy3_mersi_SST_data  = L2_path+date_fn+'/'+'FY3D_MERSI_ORBT_L2_SST_MLT_NUL_'+date_fn+'_'+hm+'_1000M_MS.HDF'
		fy3_intermediate    = L2_path+date_fn+'/'+'FY3D_MERSI_ORBT_L2_XXX_MLT_NUL_'+date_fn+'_'+hm+'_INTERMED.HDF'

		# oisst date   ;sst.day.mean.20180201.hdf5
		[y2,m2,d2] = julian_to_date(jd-0) #(jd-2)
		print ('Date oisst = ', y2, m2, d2)
		yn2 = int_to_str(y2)
		mn2 = int_to_str(m2)
		dn2 = int_to_str(d2)
		oisst_data = sst_path+'sst.day.mean.'+yn2+mn2+dn2+'.hdf5'
		sst_id = 1
		file_id = os.path.exists(oisst_data)
		if file_id == True:
			s1 = os.path.getsize(oisst_data)
			s1 = s1/1024.0/1024.0
			print ('s1 size = ',s1)
			if (s1>0.3):
				sst_id = 1
				
		if (sst_id ==0):
			[y2,m2,d2] = julian_to_date(jd)
			print ('Date oisst = ', y2, m2, d2)
			yn2 = '3000'
			mn2 = int_to_str(m2)
			dn2 = int_to_str(d2)
			oisst_data = sst_path+'sst.day.mean.'+yn2+mn2+dn2+'.hdf5'
		
		if nwp_id == 5:
			# nwp   ;gfs.t06z.pgrb2.0p50.f021
			print ('find gfs =', nwp_id)
			[nwp1,nwp2] = find_nwp_name(hour)
			print ('nwp time = ',nwp1,' ',nwp2)
			nwp_grib_data1 = nwp_path+date_fn+'/ORG/gfs.t06z.pgrb2.0p50.f0'+nwp1
			nwp_grib_data2 = nwp_path+date_fn+'/ORG/gfs.t06z.pgrb2.0p50.f0'+nwp2
			
		if nwp_id == 3: 
			print ('find t639 =', nwp_id,year_fn, date_fn, nwp_path, hour)
			# nwp   ;'gmf.639.2018093012012.grb2'
			[nwp1,nwp2] = find_nwp_name2(year_fn, date_fn, nwp_path, hour)
			print ('nwp time = ',nwp1,' ',nwp2)
			nwp_grib_data1 = nwp_path+date_fn+'/ORG/'+nwp1
			nwp_grib_data2 = nwp_path+date_fn+'/ORG/'+nwp2

		if nwp_id == 4: 
			print ('find ncep grib2 =', nwp_id)
			#date_fn_p1 fnl_20190419_18_00.grib2
			[nwp1,nwp2] = find_nwp_name4(hour)
			print ('nwp time = ',nwp1,' ',nwp2)
			nwp_grib_data1 = nwp_path+date_fn+'/ORG/fnl_'+date_fn+'_'+nwp1+'_00.grib2'
			nwp_grib_data2 = nwp_path+date_fn+'/ORG/fnl_'+date_fn+'_'+nwp2+'_00.grib2'
			if nwp2 == '24':
				nwp_grib_data2 = nwp_path+date_fn_p1+'/ORG/fnl_'+date_fn_p1+'_00_00.grib2'
			print (nwp_grib_data1)
			print (nwp_grib_data2)

		if nwp_id == 6: 
			print ('find grapes =', nwp_id,year_fn, date_fn, nwp_path, hour)
			# nwp   ;'gmf.gra.2018093012012.grb2'
			[nwp1,nwp2] = find_nwp_name3(year_fn, date_fn, nwp_path, hour)
			print ('nwp time = ',nwp1,' ',nwp2)
			nwp_grib_data1 = nwp_path+date_fn+'/ORG/'+nwp1
			nwp_grib_data2 = nwp_path+date_fn+'/ORG/'+nwp2

		if nwp_id == 7:
    			# nwp   ;gdas.t06z.pgrb2.0p50.f021
			print ('find gfs =', nwp_id)
			[nwp1,nwp2] = find_nwp_name(hour)
			print ('nwp time = ',nwp1,' ',nwp2)
			nwp_grib_data1 = nwp_path+date_fn+'/ORG/gdas.t06z.pgrb2.0p25.f0'+nwp1
			nwp_grib_data2 = nwp_path+date_fn+'/ORG/gdas.t06z.pgrb2.0p25.f0'+nwp2

		if nwp_id == 8 or nwp_id == 10:
    			# nwp   ;gfs.t06z.pgrb2.0p50.f021
			print ('find gfs =', nwp_id)
			[nwp1,nwp2] = find_nwp_name(hour)
			print ('nwp time = ',nwp1,' ',nwp2)
			nwp_grib_data1 = nwp_path+date_fn+'/ORG/gfs.t06z.pgrb2.0p25.f0'+nwp1
			nwp_grib_data2 = nwp_path+date_fn+'/ORG/gfs.t06z.pgrb2.0p25.f0'+nwp2

		if nwp_id == 9:
    			# nwp   ;gfs.t06z.pgrb2.0p50.f021
			print ('find gfs =', nwp_id)
			[nwp1,nwp2] = find_nwp_name(hour)
			print ('nwp time = ',nwp1,' ',nwp2)
			nwp_grib_data1 = nwp_path+date_fn+'/ORG/gfs.t06z.pgrb2.0p50.f0'+nwp1
			nwp_grib_data2 = nwp_path+date_fn+'/ORG/gfs.t06z.pgrb2.0p50.f0'+nwp2

	else:
		print ('ERROR: No GEO data!')
		#continue
		
	#i = i+1
	#++++++++++
	if (retr_id == 1):
		print ('Write fy3d config control file')
		print ('temp_fy3d_config_'+date_fn+'_'+hm+'.nml')
		f = open('temp_fy3d_config_'+date_fn+'_'+hm+'.nml','w+')
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
		f.write('  fy3_mersi_SST_data  = "'+fy3_mersi_SST_data+'",\n')
		f.write('  fy3_intermediate    = "'+fy3_intermediate+'",\n')
		f.write('  fylat_nwp_opt       = '+str(nwp_id)+',\n')
		f.write('  fylat_rtm_opt       = 0'+',\n')
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
		f.write('  surface_sst_id      = 1'+',\n')
		f.write('  write_inter_id      = 0/'+'\n')
		f.close()
		print ('Drive fy3 mersi_ii FYLAT')
		os.chdir(code_path)
		os.system('./fylat_FY3_MERSI_II_PGS temp_fy3d_config_'+date_fn+'_'+hm+'.nml')
		#os.system('rm temp_fy3d_config_'+date_fn+'_'+hm+'.nml')
		os.chdir(code_path)

	
def convert_nwp(yn,mn,dn,path,bin_path,ww): 

	print ('convert nwp grib2 to binary data')
	file_list = glob.glob(path+'gfs.t06z.pgrb2.0p50*')
	file_list = sorted(file_list)
	y = int(yn)
	m = int(mn)
	d = int(dn)
	#TRIM(TRIM(code_root_path)//'/wgrib/NCEP_grib2_to_bin.sh')
	#gfs0p50_20180408_03_00
	for filename in file_list:
		temp_fn = os.path.split(filename)
		fname = temp_fn[1]
		#gfs.t06z.pgrb2.0p50.f018
		ps = int(fname[5:7])
		tl = int(fname[21:24])
		h1 = ps + tl - 24
		if (h1 >=0) :
			hn = int_to_str(h1)
			bin_name = bin_path+'gfs0p50_'+yn+mn+dn+'_'+hn+'_00'
			print (filename)
			print (bin_name)
			file_bin = glob.glob(bin_name)
			print ('file id = ',len(file_bin))
			if (len(file_bin) < 1):
				os.system(ww+' '+filename+' '+bin_name)
		#call system(wgrib2_exe//' '//TRIM(grib_name)//' '//TRIM(bin_name))
	print ('convert nwp data over!!!')


def convert_nwp2(yn,mn,dn,path,bin_path,ww): 

	print ('convert nwp grib2 to binary data')
	file_list = glob.glob(path+'gmf.639.'+yn+'*grb2')
	file_list = sorted(file_list)
	y = int(yn)
	m = int(mn)
	d = int(dn)
	#TRIM(TRIM(code_root_path)//'/wgrib/NCEP_grib2_to_bin.sh')
	#gfs0p50_20180408_03_00
	for filename in file_list:
		temp_fn = os.path.split(filename)
		fname = temp_fn[1]
		path_nwp0 = temp_fn[0]
		len1 = len(path_nwp0)
		path_nwp = path_nwp0[0:len1-4]
		#print ('path_nwp = ',path_nwp)
		#gmf.639.2018093012024.grb2
		tnn = fname[8:16]
		ps = int(fname[16:18])
		tl = int(fname[18:21])
		#print ('tnn ps tl = ', tnn, ps, tl)
		h1 = ps + tl - 24
		if (h1 >=0) :
			hn = int_to_str(h1)
			bin_name = bin_path+'T639_'+yn+mn+dn+'_'+hn+'_00'
			#print ('filename = ',filename)
			print ('bin_name = ',bin_name)
			file_bin = glob.glob(bin_name)
			print ('file id = ',len(file_bin))
			if (len(file_bin) < 1):
				print (ww+' '+filename+' '+path_nwp)
				os.system(ww+' '+filename+' '+path_nwp)
				
		#call system(wgrib2_exe//' '//TRIM(grib_name)//' '//TRIM(bin_name))
	print ('convert t639 nwp data over!!!')


def convert_nwp3(yn,mn,dn,path,bin_path,ww): 

	print ('convert nwp grib2 to binary data')
	file_list = glob.glob(path+'gmf.gra.'+yn+'*grb2')
	file_list = sorted(file_list)
	y = int(yn)
	m = int(mn)
	d = int(dn)
	#TRIM(TRIM(code_root_path)//'/wgrib/NCEP_grib2_to_bin.sh')
	#gfs0p50_20180408_03_00
	for filename in file_list:
		temp_fn = os.path.split(filename)
		fname = temp_fn[1]
		path_nwp0 = temp_fn[0]
		len1 = len(path_nwp0)
		path_nwp = path_nwp0[0:len1-4]
		#print ('path_nwp = ',path_nwp)
		#gmf.639.2018093012024.grb2
		tnn = fname[8:16]
		ps = int(fname[16:18])
		tl = int(fname[18:21])
		print ('tnn ps tl = ', tnn, ps, tl)
		h1 = ps + tl - 24
		if (h1 >=0) :
			hn = int_to_str(h1)
			bin_name = bin_path+'GRAPES_GFS_'+yn+mn+dn+'_'+hn+'_00'
			#print ('filename = ',filename)
			print ('bin_name = ',bin_name)
			file_bin = glob.glob(bin_name)
			print ('file id = ',len(file_bin))
			if (len(file_bin) < 1):
				print (ww+' '+filename+' '+bin_name)
				os.system(ww+' '+filename+' '+bin_name)
				
		#call system(wgrib2_exe//' '//TRIM(grib_name)//' '//TRIM(bin_name))
	print ('convert grapes nwp data over!!!')

def convert_nwp4(yn,mn,dn,path,bin_path,ww): 

	print ('convert nwp grib2 to binary data')
	file_list = glob.glob(path+'gfs.t06z.pgrb2.0p25*')
	file_list = sorted(file_list)
	y = int(yn)
	m = int(mn)
	d = int(dn)
	#TRIM(TRIM(code_root_path)//'/wgrib/NCEP_grib2_to_bin_uv_0p25.sh')
	#gfs0p25_20180408_03_00
	for filename in file_list:
		temp_fn = os.path.split(filename)
		fname = temp_fn[1]
		#gfs.t06z.pgrb2.0p25.f018
		ps = int(fname[5:7])
		tl = int(fname[21:24])
		h1 = ps + tl - 24
		if (h1 >=0) :
			hn = int_to_str(h1)
			bin_name = bin_path+'gfs0p25_'+yn+mn+dn+'_'+hn+'_00'
			print (filename)
			print (bin_name)
			file_bin = glob.glob(bin_name)
			print ('file id = ',len(file_bin))
			if (len(file_bin) < 1):
				os.system(ww+' '+filename+' '+bin_name)
		#call system(wgrib2_exe//' '//TRIM(grib_name)//' '//TRIM(bin_name))
	print ('convert nwp data over!!!')

def convert_nwp5(yn,mn,dn,path,bin_path,ww): 

	print ('convert nwp grib2 to binary data')
	file_list = glob.glob(path+'gfs.t06z.pgrb2.0p50*')
	file_list = sorted(file_list)
	y = int(yn)
	m = int(mn)
	d = int(dn)
	#TRIM(TRIM(code_root_path)//'/wgrib/NCEP_grib2_to_bin.sh')
	#gfs0p50_20180408_03_00
	for filename in file_list:
		temp_fn = os.path.split(filename)
		fname = temp_fn[1]
		#gfs.t06z.pgrb2.0p50.f018
		ps = int(fname[5:7])
		tl = int(fname[21:24])
		h1 = ps + tl - 24
		if (h1 >=0) :
			hn = int_to_str(h1)
			bin_name = bin_path+'gfs0p50_41L_'+yn+mn+dn+'_'+hn+'_00'
			print (filename)
			print (bin_name)
			file_bin = glob.glob(bin_name)
			print ('file id = ',len(file_bin))
			if (len(file_bin) < 1):
				os.system(ww+' '+filename+' '+bin_name)
		#call system(wgrib2_exe//' '//TRIM(grib_name)//' '//TRIM(bin_name))
	print ('convert nwp data over!!!')

def convert_nwp6(yn,mn,dn,path,bin_path,ww): 

	print ('convert nwp grib2 to binary data')
	file_list = glob.glob(path+'gfs.t06z.pgrb2.0p25*')
	file_list = sorted(file_list)
	y = int(yn)
	m = int(mn)
	d = int(dn)
	#TRIM(TRIM(code_root_path)//'/wgrib/NCEP_grib2_to_bin.sh')
	#gfs0p50_20180408_03_00
	for filename in file_list:
		temp_fn = os.path.split(filename)
		fname = temp_fn[1]
		#gfs.t06z.pgrb2.0p50.f018
		ps = int(fname[5:7])
		tl = int(fname[21:24])
		h1 = ps + tl - 24
		if (h1 >=0) :
			hn = int_to_str(h1)
			bin_name = bin_path+'gfs0p25_41L_'+yn+mn+dn+'_'+hn+'_00'
			print (filename)
			print (bin_name)
			file_bin = glob.glob(bin_name)
			print ('file id = ',len(file_bin))
			if (len(file_bin) < 1):
				os.system(ww+' '+filename+' '+bin_name)
		#call system(wgrib2_exe//' '//TRIM(grib_name)//' '//TRIM(bin_name))
	print ('convert nwp data over!!!')


import sys
sys.dont_write_bytecode = True
import os 
import glob
import time
#import multiprocessing 
from multiprocessing import Pool
#import numpy as np 
#reload(sys)
#sys.setdefaultencoding('utf8')


#if __name__ == '__main__':
# ---------- MAIN ----------
# step 0: set path, time and sensor id
print ('   ')
print ('  --- step 0: input data saving path and result output path')
L1_path   = '/data/Data_minmin/fy3cloud/data/L1/1000M/'  #'/PUBLICDATAFY3PGS/FTPDATA/SAT/Data_fy3cloud/FY3D_MERSI_II_test/L1/'
GEO_path  = '/data/Data_minmin/fy3cloud/data/L1/GEO1K/' #'/PUBLICDATAFY3PGS/FTPDATA/SAT/Data_fy3cloud/FY3D_MERSI_II_test/GEO1K/'
nwp_path  = '/data/nwp/'
sst_path  = '/data/Data_minmin/oisst/'
L2_path   = '/data/Data_minmin/fy3cloud/data/L2/'
code_path = '/home/minm/Research/fy3/fy3cloud/retrieval_system_V3.0/'
print ('  L1_path   = '+L1_path)
print ('  GEO_path  = '+GEO_path)
print ('  nwp_path  = '+nwp_path)
print ('  sst_path  = '+sst_path)
print ('  L2_path   = '+L2_path)
print ('  code_path = '+code_path)
sensor_id = 21  # 21=fy3d; 22=fy3e
nwp_id    = 5  # 41=h8; 31=fy4a 
'''
  !*************************************************
  !   1 = ncep reanalysis 1*1 (grib1)
  !   2 = gfs1p00 1*1 (grib2)
  !   3 = T639 0.125*0.125 (grib2)
  !   4 = ncep reanalysis 1*1 (grib2) 
  !   5 = gfs0p50 0.5*0.5 (grib2) 
  !   6 = grapes 0.25*0.25 (grib2)
  !   7 = gdas1 0.25*0.25 (grib2)
  !   8 = gfs 0.25*0.25 (grib2)
  !   9 = gfs 0.5*0.5 @41-layers (grib2)
  !   10 = gfs 0.25*0.25 @41-layers (grib2)
  !*************************************************
'''
#        start  end
year  = [2020, 2020]
month = [   3,   3]
day   = [   30, 30]
hour_lmt  = [  0,  24]
nthread = 1 


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
	jd_p1 = jd + 1
	#convert date
	[y1,m1,d1] = julian_to_date(jd)
	[y2,m2,d2] = julian_to_date(jd_p1)
	print ('Date = ', y1, m1, d1)
	yn = int_to_str(y1)
	mn = int_to_str(m1)
	dn = int_to_str(d1)
	yn2 = int_to_str(y2)
	mn2 = int_to_str(m2)
	dn2 = int_to_str(d2)
		
	date_fn = yn+mn+dn
	date_fn_p1 = yn2+mn2+dn2
	
	year_fn = yn
	#directory_name = L1_path+date_fn+'/'
	directory_name = L1_path+year_fn+'/'+date_fn+'/'
	print ('directory_name = ',directory_name)
	
	directory_nwp = nwp_path+date_fn+'/ORG/'
	bin_nwp = nwp_path+date_fn+'/'
	if (nwp_id == 2 or nwp_id == 4 or nwp_id == 5):
		wgrib2_commend = code_path+'wgrib/NCEP_grib2_to_bin_uv.sh'
		convert_nwp(yn,mn,dn,directory_nwp,bin_nwp,wgrib2_commend)
	if (nwp_id == 3):
		wgrib2_commend = code_path+'wgrib/T639_to_bin.sh'
		convert_nwp2(yn,mn,dn,directory_nwp,bin_nwp,wgrib2_commend)
	if (nwp_id == 6):
		wgrib2_commend = code_path+'wgrib/GRAPES_GFS_grib2_to_bin_uv.sh'
		convert_nwp3(yn,mn,dn,directory_nwp,bin_nwp,wgrib2_commend)
	if (nwp_id == 7 or nwp_id == 8):
		wgrib2_commend = code_path+'wgrib/NCEP_grib2_to_bin_uv_0p25.sh'
		convert_nwp4(yn,mn,dn,directory_nwp,bin_nwp,wgrib2_commend)
	if (nwp_id == 9):
		wgrib2_commend = code_path+'wgrib/NCEP_grib2_to_bin_uv_41Layers.sh'
		convert_nwp5(yn,mn,dn,directory_nwp,bin_nwp,wgrib2_commend)
	if (nwp_id == 10):
		wgrib2_commend = code_path+'wgrib/NCEP_grib2_to_bin_uv_0p25_41Layers.sh'
		convert_nwp6(yn,mn,dn,directory_nwp,bin_nwp,wgrib2_commend)
	
	if os.path.isdir(directory_name):
	
		if sensor_id == 21:  # fy3d mersi_ii  ;FY3D_MERSI_GBAL_L1_20180220_2355_1000M_MS.HDF
			file_list = glob.glob(directory_name+'FY3D_MERSI_GBAL_L1_'+date_fn+'_*_1000M_MS.HDF')
			file_list = sorted(file_list)
		if sensor_id == 11:  # fy3d mersi_ii using modis
			file_list = glob.glob(directory_name+'FY3D_MERSI_GBAL_L1_'+date_fn+'_*_1000M_MS.HDF')
			file_list = sorted(file_list)
		
		#print(file_list)
			
		i=1
		p = Pool(nthread)
		for filename in file_list:
			#driver_fylat(filename, date_fn, jd, sensor_id, nwp_id, L1_path, GEO_path, L2_path, code_path, sst_path, nwp_path, hour_lmt ):
			#p.apply_async(write_new_data, args=(filename, out_path,))
			p.apply_async(driver_fylat, args=(filename, year_fn, date_fn, date_fn_p1, jd, sensor_id, nwp_id, L1_path, GEO_path, L2_path, code_path, sst_path, nwp_path, hour_lmt,))
		
		p.close()
		p.join()
		
	else:
		print ('find next directory!!!')

else:
    print ('day cycle over!!!!')
	
