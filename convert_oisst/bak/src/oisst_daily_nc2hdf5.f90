program oisst_daily_nc2hdf5

!---------------------------------------------------
! The propose of this program is to convert 
! oisst daily netcdf file to hdf5 format
!
! Author: Min Min
! Uint  : National Satellite Meteorological Center
!
!---------------------------------------------------


! use modules
!use data_arrays_module
use io_module

implicit none



!character(len=1000) :: temp, dummy     ! temporay name   
integer :: narg, i, j, k,ix,iy, hdf5_qa
character(len=300)  :: arg1, arg2, nc_path, out_path, outfile
character(len=  8)  :: arg3, tname1
! input time
character(len=  4)  :: year0
character(len=  2)  :: month0, day0
integer(kind=4)     :: year00, month00, day00
real(kind=8)        :: jday0
! oisst obs time
character(len=  4)  :: year1
character(len=  2)  :: month1, day1
real(kind=8)        :: jday, jday1
integer(kind=4)     :: year11, month11, day11, hour11, mint11
integer(kind=4)     :: numday

character(len=300)  :: nc_file, out_file, out_file2
!integer :: arg5, nameb
!character(len=  2)  :: bname0
integer :: bnum, bandnum

integer :: LENGTH,IERR

!--------------------------------------
!               start  
!--------------------------------------
print*,'--------------------------------------'
print*,'                start   '
print*,'--------------------------------------'
print*,' '
print*,' Convert oisst daily netcdf file (0.25*0.25 degree) to hdf5 format !!! '

! Step [1] : read input filename using command line
print*,' ' 
print*,' Step [1] : read input information using command line '
narg=IARGC()

if (narg < 3) then 
   print*, 'ERROR: args input number is wrong!'
   stop
else
   call getarg(1,arg1)
   call getarg(2,arg2)
   call getarg(3,arg3)
endif
nc_path   = arg1
out_path  = arg2
tname1    = arg3
year0  = tname1(1:4)
month0 = tname1(5:6)
day0   = tname1(7:8)
read(year0,'(I4)')year00
read(month0,'(I2)')month00
read(day0,'(I2)')day00
print*,'nc_path  = ',trim(nc_path)
print*,'out_path = ',trim(out_path)
!jday0 = julday(month0, day0, year0)
call julian (year00, month00, day00, 0, 0, jday0)
print*,'satellite obs time = ', year00, month00, day00,' julian day =',jday0
jday = jday0 - 3
call julian_to_date ( jday, year11, month11, day11, hour11, mint11)
print*,'oisst obs time = ', year11, month11, day11,' julian day =',jday
!jday1 = julday(1, 1, year)
call julian (year11, 1, 1, 0, 0, jday1)
numday = int(jday - jday1 + 1)
print*,'number of day in one year for oisst obs time =', numday

! Step [2] : find nc file name and define output hdf5 file name 
print*,' ' 
print*,' Step [2] : find nc file name and define output hdf5 file name '
!read(year11,'(A4)')year1
call ICNVRT(0,year11,year1,LENGTH,IERR)
nc_file = trim(nc_path)//'sst.day.mean.'//trim(year1)//'.v2.nc'
print*,'nc_file = ',trim(nc_file)
call ICNVRT(0,month11,month1,LENGTH,IERR)
if (month11 < 10) then 
   month1 = trim('0'//month1)
else
   month1 = month1
endif
call ICNVRT(0,day11,day1,LENGTH,IERR)
if (day11 < 10) then 
   day1 = trim('0'//day1)
else
   day1 = trim(day1)
endif
out_file  = trim(out_path)//'temp/sst.day.mean.'//year1//month1//day1//'.hdf5'
out_file2 = trim(out_path)//'sst.day.mean.'//year1//month1//day1//'.hdf5'
print*,'out_file = ',trim(out_file)


! Step [3] : read nc file 
print*,' ' 
print*,' Step [3] : read nc file  '
!sst0 = fltarr(1440,720,366)
call read_nc_file(nc_file,'sst',numday)


! Step [4] : write daily hdf5 oisst file 
print*,' ' 
print*,' Step [4] : write daily hdf5 oisst file   '
call write_hdf5_file(out_file,'sst')

!----------------
! Compress HDF5 File
hdf5_qa = 1
if (hdf5_qa == 1) then ! hdf5 file is ok
   print*,'Compress HDF5 File Start'
   CALL system('h5repack -f GZIP=5 -v '//trim(out_file)//' '//trim(out_file2)//'>/dev/null')
   CALL system('rm -rf '//trim(out_file))
   print*,'Compress HDF5 File End'
endif
!--------------------------------------
!               end  
!--------------------------------------
print*,' '
print*,'--------------------------------------'
print*,'                end   '
print*,'--------------------------------------'
    
end
