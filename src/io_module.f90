module io_module


!C-----------------------------------------------------------------------
!C !F90                                                                  
!C
!C !Description: 
!C    This module is to FY3/MERSI-II data.
!C
!C !Input parameters
!C    none
!C 
!C !Output parameters
!C    none
!C
!C !Author's information
!C    Author: Min Min
!C    E-mail: minmin@cma.gov.cn
!C    Tel   : 86-010-68406763
!C    National Satellite Meteorological Center 
!C  
!C !END
!C----------------------------------------------------------------------

! USE modules
use data_arrays_module
use names_module
use planck_module
use constant
use platform_module
use HDF5


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


CONTAINS
!+++++++++++++++++++ step 2: subroutines +++++++++++++++++++++++++++++++
!~~~~~~~~~~~~~~~~ Subroutine  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine fylat_read_fy3_mersi_geo_data(filename)

character(len=1000) :: filename

!--------------- define variables ----------------------
integer ::  error                         ! error 

! Science Data Set id 
integer (HID_T) :: file_id              ! file id 
integer (HID_T) :: sds_id_var           ! sds id all variables
integer (HID_T) :: attr_id              ! attr id all variables
integer (HID_T) :: asp_id               !   [unit:degree]
integer (HSIZE_T), dimension(1)         :: dims_attr
integer (HSIZE_T), dimension(2)         :: dims_geo2  


! buffer variables
integer(kind=1), allocatable, dimension(:,:)      :: var_int1
integer(kind=4), allocatable, dimension(:,:)      :: var_int4
real(kind=4),    allocatable, dimension(:,:)      :: var_real4
character,       allocatable, dimension(:,:)      :: var_char
real(kind=4)                                      :: interp, slope
real(kind=4)                                      :: vzar,  szar, razr, cossna, MDPHI
integer(kind=4)                                   :: i,j,k,nLine0


!--------------- read fy3 mersi-II variables ----------------------
print*,'  ... read GEO HDF5 data'
! allocate  variables
!call allocate_fylat_fy3mersi_geo_data()

!allocate (var_int1(sat%nElem, sat%nLine),  &
!          var_int4(sat%nElem, sat%nLine),  &
!          var_real4(sat%nElem, sat%nLine), &
!          var_char(sat%nElem, sat%nLine))

! open geo hdf5 file
call h5open_f(error)
print*,'geo=',trim(filename)
call h5fopen_f (trim(filename), H5F_ACC_RDONLY_F, file_id, error)  
!call h5fopen_f ('/PUBLICDATAFY3PGS/FTPDATA/SAT/Data_fy3cloud/MERSI_II2/2014_06_30/Geo/FY3D_MERSI_GBAL_L1_20140630_0350_GEO1K_MS.HDF', H5F_ACC_RDONLY_F, file_id, error)

call h5aopen_f (file_id, 'Number Of Scans', attr_id, error)
call h5aread_f (attr_id , H5T_NATIVE_INTEGER, nLine0, dims_attr, error) !
if (error /= 0) then 
    print*,'ERROR: A dataset Number Of Scans read failed !'
    stop
endif
sat%nLine = nLine0*10
call h5aclose_f(attr_id , error)    ! close
print *, 'New number of scans  = ', nLine0,'*10 = ', sat%nLine

! allocate  variables
call allocate_fylat_fy3mersi_geo_data()

allocate (var_int1(sat%nElem, sat%nLine),  &
          var_int4(sat%nElem, sat%nLine),  &
          var_real4(sat%nElem, sat%nLine), &
          var_char(sat%nElem, sat%nLine))

      
! (1) read lon/lat data
call h5dopen_f (file_id, '/Geolocation/Longitude', sds_id_var, error)
call h5dread_f (sds_id_var, H5T_NATIVE_REAL, var_real4, dims_geo2, error) !
if (error /= 0) then 
    print*,'ERROR: dataset Longitude read failed !'
    stop
endif
geo%lon(:,:) = var_real4
call h5dclose_f(sds_id_var, error)    ! close 
var_real4 = 0.
!print*,ref1(1:5,1:5,1)

call h5dopen_f (file_id, '/Geolocation/Latitude', sds_id_var, error)
call h5dread_f (sds_id_var, H5T_NATIVE_REAL, var_real4, dims_geo2, error) ! 
if (error /= 0) then 
    print*,'ERROR: dataset Latitude read failed !'
    stop
endif
geo%lat(:,:) = var_real4
call h5dclose_f(sds_id_var, error)    ! close 
var_real4 = 0.

! (2) read geometry data
!SZA
call h5dopen_f (file_id, '/Geolocation/SolarZenith', sds_id_var, error)
call h5dread_f (sds_id_var, H5T_NATIVE_INTEGER, var_int4, dims_geo2, error) ! sza
if (error /= 0) then 
    print*,'ERROR: dataset SolarZenith read failed !'
    stop
endif
  
! read attribute start
call h5aopen_f (sds_id_var, 'Slope', attr_id, error)
call h5aread_f (attr_id, H5T_NATIVE_REAL, slope, dims_attr, error)
if (error /= 0) then 
    print*,'ERROR: attribute SolarZenith/Slope read failed !'
    stop
endif
call h5aclose_f (attr_id, error)
call h5aopen_f (sds_id_var, 'Intercept', attr_id, error)
call h5aread_f (attr_id, H5T_NATIVE_REAL, interp, dims_attr, error)
if (error /= 0) then 
    print*,'ERROR: attribute SolarZenith/Intercept read failed !'
    stop
endif
call h5aclose_f (attr_id, error)
! read attribute end
geo%SolarZenith(:,:) = (var_int4+interp)*slope
call h5dclose_f(sds_id_var, error)    ! close 
var_int4 = 0

 
!SAA
call h5dopen_f (file_id, '/Geolocation/SolarAzimuth', sds_id_var, error)
call h5dread_f (sds_id_var, H5T_NATIVE_INTEGER, var_int4, dims_geo2, error) ! sza
if (error /= 0) then 
    print*,'ERROR: dataset SolarAzimuth read failed !'
    stop
endif
! read attribute start
call h5aopen_f (sds_id_var, 'Slope', attr_id, error)
call h5aread_f (attr_id, H5T_NATIVE_REAL, slope, dims_attr, error)
if (error /= 0) then 
    print*,'ERROR: attribute SolarAzimuth/Slope read failed !'
    stop
endif
call h5aclose_f (attr_id, error)
call h5aopen_f (sds_id_var, 'Intercept', attr_id, error)
call h5aread_f (attr_id, H5T_NATIVE_REAL, interp, dims_attr, error)
if (error /= 0) then 
    print*,'ERROR: attribute SolarAzimuth/Intercept read failed !'
    stop
endif
call h5aclose_f (attr_id, error)
! read attribute end
geo%SolarAzimuth(:,:) = (var_int4+interp)*slope
call h5dclose_f(sds_id_var, error)    ! close 
var_int4 = 0

! new fy3d mersi_II data
if (fylat_sensor_id >= 2 .and. fylat_sensor_id < 10) then 
   where(geo%SolarAzimuth >180.0) geo%SolarAzimuth = geo%SolarAzimuth - 360.0  
endif

!VZA
call h5dopen_f (file_id, '/Geolocation/SensorZenith', sds_id_var, error)
call h5dread_f (sds_id_var, H5T_NATIVE_INTEGER, var_int4, dims_geo2, error) ! sza
if (error /= 0) then 
    print*,'ERROR: dataset SensorZenith read failed !'
    stop
endif
! read attribute start
call h5aopen_f (sds_id_var, 'Slope', attr_id, error)
call h5aread_f (attr_id, H5T_NATIVE_REAL, slope, dims_attr, error)
if (error /= 0) then 
    print*,'ERROR: attribute SensorZenith/Slope read failed !'
    stop
endif
call h5aclose_f (attr_id, error)
call h5aopen_f (sds_id_var, 'Intercept', attr_id, error)
call h5aread_f (attr_id, H5T_NATIVE_REAL, interp, dims_attr, error)
if (error /= 0) then 
    print*,'ERROR: attribute SensorZenith/Intercept read failed !'
    stop
endif
call h5aclose_f (attr_id, error)
! read attribute end
geo%SensorZenith(:,:) = (var_int4+interp)*slope
call h5dclose_f(sds_id_var, error)    ! close 
var_int4 = 0

!VAA
call h5dopen_f (file_id, '/Geolocation/SensorAzimuth', sds_id_var, error)
call h5dread_f (sds_id_var, H5T_NATIVE_INTEGER, var_int4, dims_geo2, error) ! sza
if (error /= 0) then 
    print*,'ERROR: dataset SensorAzimuth read failed !'
    stop
endif
! read attribute start
call h5aopen_f (sds_id_var, 'Slope', attr_id, error)
call h5aread_f (attr_id, H5T_NATIVE_REAL, slope, dims_attr, error)
if (error /= 0) then 
    print*,'ERROR: attribute SensorAzimuth/Slope read failed !'
    stop
endif
call h5aclose_f (attr_id, error)
call h5aopen_f (sds_id_var, 'Intercept', attr_id, error)
call h5aread_f (attr_id, H5T_NATIVE_REAL, interp, dims_attr, error)
if (error /= 0) then 
    print*,'ERROR: attribute SensorAzimuth/Intercept read failed !'
    stop
endif
call h5aclose_f (attr_id, error)
! read attribute end
geo%SensorAzimuth(:,:) = (var_int4+interp)*slope
call h5dclose_f(sds_id_var, error)    ! close 
var_int4 = 0

! new fy3d mersi_II data
if (fylat_sensor_id >= 2 .and. fylat_sensor_id < 10) then 
   where(geo%SensorAzimuth >180.0) geo%SensorAzimuth = geo%SensorAzimuth - 360.0  
endif


! (3) read dem and lsm
call h5dopen_f (file_id, '/Geolocation/DEM', sds_id_var, error)
call h5dread_f (sds_id_var, H5T_NATIVE_INTEGER, var_int4, dims_geo2, error) !
if (error /= 0) then 
    print*,'ERROR: dataset DEM read failed !'
    stop
endif
geo%dem(:,:) = var_int4*1.0
call h5dclose_f(sds_id_var, error)    ! close 
var_int4 = 0
!print*,ref1(1:5,1:5,1)

call h5dopen_f (file_id, '/Geolocation/LandSeaMask', sds_id_var, error)
call h5dread_f (sds_id_var, H5T_STD_U8LE, var_char, dims_geo2, error) ! 
if (error /= 0) then 
    print*,'ERROR: dataset LandSeaMask read failed !'
    stop
endif
geo%lsm(:,:) = ichar(var_char)
call h5dclose_f(sds_id_var, error)    ! close 
!print*,geo%lsm(1:3,1:5) 

! close file
call h5fclose_f(file_id, error)       ! close 
call h5close_f(error)
 
deallocate (var_int1,  &
            var_int4,  &
            var_real4, &
            var_char)


! (4) calculate relative azimuth and glint angle
! ...    get relative azimuth angle
if (fylat_sensor_id == 1) then  ! -180 to 180 (MODIS)
   do i = 1 , sat%nElem
   do j = 1 , sat%nLine
      geo%RelAzimuth(i,j) = abs(180.0-abs(geo%SensorAzimuth(i,j)-geo%SolarAzimuth(i,j)))
   enddo
   enddo
endif

if (fylat_sensor_id >= 2 ) then   ! 0-360 (FY-3) proxy
   do i = 1 , sat%nElem
   do j = 1 , sat%nLine
      MDPHI = abs(geo%SolarAzimuth(i,j)- geo%SensorAzimuth(i,j))
      IF(MDPHI.GT.360.0) MDPHI=amod(MDPHI,360.0)          
      IF(MDPHI.GT.180.0) MDPHI=360.0-MDPHI
      !MDPHI = abs(180.0-abs(geo%SensorAzimuth(i,j)-geo%SolarAzimuth(i,j)))
      !geo%RelAzimuth(i,j) = MDPHI
      geo%RelAzimuth(i,j) = abs(180.0-MDPHI)     
   enddo
   enddo
endif

!if (fylat_sensor_id == 21 .or. fylat_sensor_id == 22) then  ! 0 to 180 (fy3d/mersi_ii)
!   do i = 1 , sat%nElem
!   do j = 1 , sat%nLine
!      !MDPHI = abs(geo%SolarAzimuth(i,j)- geo%SensorAzimuth(i,j))
!      !IF(MDPHI.GT.360.0) MDPHI=amod(MDPHI,360.0)          
!      !IF(MDPHI.GT.180.0) MDPHI=360.0-MDPHI
!      geo%RelAzimuth(i,j) = MDPHI
!   enddo
!   enddo
!endif

! ...    Now calculate the relative angle (value that sun glint is based upon.)
do i = 1 , sat%nElem
do j = 1 , sat%nLine
   if (nint(geo%SensorZenith(i,j))  .ne. nint(missing_value_real4) .and.  &
       nint(geo%SolarAzimuth(i,j))  .ne. nint(missing_value_real4) .and.  &
       nint(geo%SolarZenith(i,j))   .ne. nint(missing_value_real4) .and.  &
       nint(geo%SensorAzimuth(i,j)) .ne. nint(missing_value_real4)) then
       vzar = geo%SensorZenith(i,j) * DTOR
       szar = geo%SolarZenith(i,j) * DTOR
       razr = geo%RelAzimuth(i,j) * DTOR
       cossna = min( (sin(vzar) * sin(szar) * cos(razr) + cos(vzar) * cos(szar)) , 1.0)
       geo%GlintAngle(i,j) = acos(cossna) * RTOD
   endif
enddo
enddo
   
print*,'longitude test    = ',geo%lon(1000:1003,1000:1003)
print*,'latitude test     = ',geo%lat(1000:1003,1000:1003)
print*,'DEM data test     = ',geo%dem(1000:1003,1000:1003)
print*,'SolZen angle test = ',geo%SolarZenith(1000:1003,1000:1003)
print*,'SenZen angle test = ',geo%SensorZenith(1000:1003,1000:1003)
print*,'SolAzi angle test = ',geo%SolarAzimuth(1000:1003,1000:1003)
print*,'SenAzi angle test = ',geo%SensorAzimuth(1000:1003,1000:1003)
print*,'ReAzi angle test  = ',geo%RelAzimuth(1000:1003,1000:1003)
print*,'Glint angle test  = ',geo%GlintAngle(1000:1003,1000:1003)

end subroutine fylat_read_fy3_mersi_geo_data
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!~~~~~~~~~~~~~~~~ Subroutine  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine fylat_read_fy3_mersi_L1b_data(filename)

character(len=1000) :: filename
character(len=100)  :: dummy
integer             :: IC

!--------------- define variables ----------------------
integer   error                         ! error 

! Science Data Set id 
integer (HID_T) :: file_id              ! file id 
integer (HID_T) :: sds_id_var           ! sds id all variables
integer (HID_T) :: attr_id              ! attr id all variables
integer (HID_T) :: asp_id               !   [unit:degree]
integer (HSIZE_T), dimension(1)         :: dims_attr
integer (HSIZE_T), dimension(2)         :: dims_geo2  
integer (HSIZE_T), dimension(3)         :: dims_geo3

! buffer variables
integer(kind=4), allocatable, dimension(:,:,:)      :: var_int_tmp
integer(kind=4), allocatable, dimension(:,:,:)      :: var_int4
real(kind=4), allocatable, dimension(:,:)           :: var_real4_ti
real(kind=4), allocatable, dimension(:,:,:)         :: var_real4
real(kind=4),    dimension( 3, 19)                  :: var_real4_vis
real(kind=4),    dimension( 3,  6)                  :: var_real4_ir
! operational format
real(kind=4),    dimension( 200, 4,  6)             :: var_real4_ir2

real(kind=4)                                        :: interp, slope
real(kind=4), dimension(2)                          :: interp2, slope2
real(kind=4), dimension(4)                          :: interp4, slope4
integer(kind=4)                                     :: i,j,k

integer(kind=4)                                     :: vis_id ! added by minmin

!--------------- read fy3 mersi-II variables ----------------------
print*,'  ... read L1b HDF5 data'
! allocate  variables
call allocate_fylat_fy3mersi_L1b_data()

do i = 1 , sat%nElem
do j = 1 , sat%nLine
   sat%zsfc(i,j) = geo%dem(i,j) 
enddo
enddo
!sat%zsfc(:,:) = geo%dem(:,:)
print*,'L1b=',trim(filename)

! open geo hdf5 file
call h5open_f(error)
call h5fopen_f (trim(filename), H5F_ACC_RDONLY_F, file_id, error)  
sat%vis_cal_coef = 0.0
sat%ir_cal_coef  = 0.0

! (1) read calibration coefficient 
if (fylat_sensor_id == 1) then !'Sensor_id=1 / FY3D-MERSI_II data (convert modis to mersi_II)'

   call h5dopen_f (file_id, '/Calibration/VIS_Cal_Coeff', sds_id_var, error)
   call h5dread_f (sds_id_var, H5T_NATIVE_REAL, var_real4_vis, dims_geo2, error) !
   if (error /= 0) then 
       print*,'ERROR: dataset VIS_Cal_Coeff read failed !'
       stop
   endif
   sat%vis_cal_coef(1:3,:) = var_real4_vis(1:3,:)
   call h5dclose_f(sds_id_var, error)    ! close 

   call h5dopen_f (file_id, '/Calibration/IR_Cal_Coeff', sds_id_var, error)
   call h5dread_f (sds_id_var, H5T_NATIVE_REAL, var_real4_ir, dims_geo2, error) !
   if (error /= 0) then 
       print*,'ERROR: dataset IR_Cal_Coeff read failed !'
       stop
   endif
   sat%ir_cal_coef(1:3,:) = var_real4_ir(1:3,:)
   call h5dclose_f(sds_id_var, error)    ! close 
   
endif

if (fylat_sensor_id == 2 .or. fylat_sensor_id == 21) then  !'Sensor_id=2 / FY3D-MERSI_II data (convert modis to mersi_II in mersi_II format)'

   call h5dopen_f (file_id, '/Calibration/VIS_Cal_Coeff', sds_id_var, error)
   call h5dread_f (sds_id_var, H5T_NATIVE_REAL, var_real4_vis, dims_geo2, error) !
   if (error /= 0) then 
       print*,'ERROR: dataset VIS_Cal_Coeff read failed !'
       stop
   endif
   sat%vis_cal_coef(1:3,:) = var_real4_vis(1:3,:)
   call h5dclose_f(sds_id_var, error)    ! close 
   
   vis_id = 0 ! install new calibration coefficients
   if (vis_id == 1) then
   	  ! --- B02
      !sat%vis_cal_coef(1,2) = -3.4390
      !sat%vis_cal_coef(2,2) = 0.0255
      ! --- B03
      !sat%vis_cal_coef(1,3) = -6.7329
      !sat%vis_cal_coef(2,3) = 0.0267
      ! --- B05
      !sat%vis_cal_coef(1,5) = -3.3302
      !sat%vis_cal_coef(2,5) = 0.0232
      ! --- B07
      !sat%vis_cal_coef(1,7) = -3.2447
      !sat%vis_cal_coef(2,7) = 0.0207
      ! --- B08
      !sat%vis_cal_coef(1,8) = -1.9084
      !sat%vis_cal_coef(2,8) = 0.0101
      open(22,file='VIS_Cal_Coeff.xcfg',status='old')
      read(22,*) dummy 
      do i = 1, 19
         read(22,*) IC,sat%vis_cal_coef(1,i),sat%vis_cal_coef(2,i)
         print*,'test = ',IC,sat%vis_cal_coef(1,i),sat%vis_cal_coef(2,i)
      enddo
      close(22)
   endif
   
   call h5dopen_f (file_id, '/Calibration/IR_Cal_Coeff', sds_id_var, error)
   call h5dread_f (sds_id_var, H5T_NATIVE_REAL, var_real4_ir2, dims_geo2, error) !
   if (error /= 0) then 
       print*,'ERROR: dataset IR_Cal_Coeff read failed !'
       stop
   endif
   sat%ir_cal_coef(1:3,:) = var_real4_ir2(100,1:3,:)
   call h5dclose_f(sds_id_var, error)    ! close 

endif

!print*,'vis cal=',sat%vis_cal_coef(1:4,8)
!print*,'ir  cal=',sat%ir_cal_coef(1:4,2)
        
! (2) read ref/rad
!		level1b_buffer = (temp*1.0 - offsets(band_index+1)) *scale_factors(band_index+1)
! visible band 1-4
allocate (var_int4(sat%nElem, sat%nLine, 4))
call h5dopen_f (file_id, '/Data/EV_250_Aggr.1KM_RefSB', sds_id_var, error)
call h5dread_f (sds_id_var, H5T_NATIVE_INTEGER, var_int4, dims_geo3, error) !
if (error /= 0) then 
    print*,'ERROR: dataset EV_250_Aggr.1KM_RefSB read failed !'
    stop
endif

allocate (var_real4_ti(sat%nElem, sat%nLine))
do i = 1, 4 
   var_real4_ti = float(var_int4(:,:,i))
   !sat%ref(:,:,i) = (var_int4(:,:,i)*1.0-sat%vis_cal_coef(1,i)) * sat%vis_cal_coef(2,i)
   call dn_to_ref_rad(var_real4_ti(:,:),      &
                      sat%vis_cal_coef(1,i),  &
                      sat%vis_cal_coef(2,i),  &
                      sat%vis_cal_coef(3,i),  &
                      fylat_sensor_id,        &
                      sat%ref_vis(:,:,i),     &
                      0)
enddo
!deallocate(var_real4_ti)
call h5dclose_f(sds_id_var, error)    ! close 
deallocate (var_int4)
print*,'vis coeff = ',sat%vis_cal_coef(1,2),sat%vis_cal_coef(2,2), sat%vis_cal_coef(3,2)
print*,'ref1 test = ',sat%ref_vis(1000:1003,1000:1003,1)
print*,'ref2 test = ',sat%ref_vis(1000:1003,1000:1003,2)
print*,'ref3 test = ',sat%ref_vis(1000:1003,1000:1003,3)
print*,'ref4 test = ',sat%ref_vis(1000:1003,1000:1003,4)

! visible band 5-19
allocate (var_int4(sat%nElem, sat%nLine, 15))
call h5dopen_f (file_id, '/Data/EV_1KM_RefSB', sds_id_var, error)
call h5dread_f (sds_id_var, H5T_NATIVE_INTEGER, var_int4, dims_geo3, error) ! 
if (error /= 0) then 
    print*,'ERROR: dataset EV_1KM_RefSB read failed !'
    stop
endif
!allocate (var_real4_ti(sat%nElem, sat%nLine))
do i = 5, 19
   var_real4_ti = float(var_int4(:,:,i-4))
   call dn_to_ref_rad(var_real4_ti(:,:),      &
                      sat%vis_cal_coef(1,i),  &
                      sat%vis_cal_coef(2,i),  &
                      sat%vis_cal_coef(3,i),  &
                      fylat_sensor_id,        &
                      sat%ref_vis(:,:,i),     &
                      0)   ! CHECK by minmin
enddo
call h5dclose_f(sds_id_var, error)    ! close 
deallocate (var_int4)
deallocate(var_real4_ti)

if (fylat_sensor_id == 1 .or. fylat_sensor_id == 2) then

! ir band 20-23
allocate (var_real4(sat%nElem, sat%nLine, 6))
allocate (var_int_tmp(sat%nElem, sat%nLine, 4))
call h5dopen_f (file_id, '/Data/EV_1KM_Emissive', sds_id_var, error)
!call h5dread_f (sds_id_var, H5T_NATIVE_INTEGER, var_int4(1:sat%nElem,1:sat%nLine,1:4), dims_geo3, error) ! 
call h5dread_f (sds_id_var, H5T_NATIVE_INTEGER, var_int_tmp, dims_geo3, error) ! 

do i = 1,4 
   var_real4(1:sat%nElem,1:sat%nLine,i) = float(var_int_tmp(1:sat%nElem,1:sat%nLine,i))
enddo
!print*,'dn20',var_real4(1000:1003,1000:1003,1),var_int_tmp(1000:1003,1000:1003,1)
!print*,'dn21',var_real4(1000:1003,1000:1003,2),var_int_tmp(1000:1003,1000:1003,2)
deallocate(var_int_tmp)
if (error /= 0) then 
    print*,'ERROR: dataset EV_1KM_Emissive read failed !'
    stop
endif
call h5dclose_f(sds_id_var, error)    ! close 


! ir band 24-25
allocate (var_int_tmp(sat%nElem, sat%nLine, 2))
call h5dopen_f (file_id, '/Data/EV_250_Aggr.1KM_Emissive', sds_id_var, error)
!call h5dread_f (sds_id_var, H5T_NATIVE_INTEGER, var_int4(1:sat%nElem,1:sat%nLine,5:6), dims_geo3, error) ! 
call h5dread_f (sds_id_var, H5T_NATIVE_INTEGER, var_int_tmp, dims_geo3, error) ! 
do i = 1,2
   var_real4(1:sat%nElem,1:sat%nLine,4+i) = float(var_int_tmp(1:sat%nElem,1:sat%nLine,i))
enddo
deallocate(var_int_tmp)
if (error /= 0) then 
    print*,'ERROR: dataset EV_250_Aggr.1KM_Emissive read failed !'
    stop
endif
!call h5dclose_f(sds_id_var, error)    ! close 
endif


if (fylat_sensor_id == 21) then

! ir band 20-23
allocate (var_real4(sat%nElem, sat%nLine, 6))
allocate (var_int_tmp(sat%nElem, sat%nLine, 4))
call h5dopen_f (file_id, '/Data/EV_1KM_Emissive', sds_id_var, error)
!call h5dread_f (sds_id_var, H5T_NATIVE_INTEGER, var_int4(1:sat%nElem,1:sat%nLine,1:4), dims_geo3, error) ! 
call h5dread_f (sds_id_var, H5T_NATIVE_INTEGER, var_int_tmp, dims_geo3, error) ! 
! read attribute start
call h5aopen_f (sds_id_var, 'Slope', attr_id, error)
call h5aread_f (attr_id, H5T_NATIVE_REAL, slope4, dims_attr, error)
if (error /= 0) then 
    print*,'ERROR: attribute /Data/EV_1KM_Emissive/Slope read failed !'
    stop
endif
call h5aclose_f (attr_id, error)
call h5aopen_f (sds_id_var, 'Intercept', attr_id, error)
call h5aread_f (attr_id, H5T_NATIVE_REAL, interp4, dims_attr, error)
if (error /= 0) then 
    print*,'ERROR: attribute /Data/EV_1KM_Emissive/Intercept read failed !'
    stop
endif
call h5aclose_f (attr_id, error)
print*,'TBB slope = ', slope4(1:4),interp4(1:4)
do i = 1,4 
   var_real4(1:sat%nElem,1:sat%nLine,i) = (float(var_int_tmp(1:sat%nElem,1:sat%nLine,i))+interp4(i))*slope4(i)
enddo
!print*,'dn20',var_real4(1000:1003,1000:1003,1),var_int_tmp(1000:1003,1000:1003,1)
!print*,'dn21',var_real4(1000:1003,1000:1003,2),var_int_tmp(1000:1003,1000:1003,2)
deallocate(var_int_tmp)
if (error /= 0) then 
    print*,'ERROR: dataset EV_1KM_Emissive read failed !'
    stop
endif
call h5dclose_f(sds_id_var, error)    ! close 

! ir band 24-25
allocate (var_int_tmp(sat%nElem, sat%nLine, 2))
call h5dopen_f (file_id, '/Data/EV_250_Aggr.1KM_Emissive', sds_id_var, error)
!call h5dread_f (sds_id_var, H5T_NATIVE_INTEGER, var_int4(1:sat%nElem,1:sat%nLine,5:6), dims_geo3, error) ! 
call h5dread_f (sds_id_var, H5T_NATIVE_INTEGER, var_int_tmp, dims_geo3, error) ! 
! read attribute start
call h5aopen_f (sds_id_var, 'Slope', attr_id, error)
call h5aread_f (attr_id, H5T_NATIVE_REAL, slope2, dims_attr, error)
if (error /= 0) then 
    print*,'ERROR: attribute /Data/EV_1KM_Emissive/Slope read failed !'
    stop
endif
call h5aclose_f (attr_id, error)
call h5aopen_f (sds_id_var, 'Intercept', attr_id, error)
call h5aread_f (attr_id, H5T_NATIVE_REAL, interp2, dims_attr, error)
if (error /= 0) then 
    print*,'ERROR: attribute /Data/EV_1KM_Emissive/Intercept read failed !'
    stop
endif
call h5aclose_f (attr_id, error)
do i = 1,2
   var_real4(1:sat%nElem,1:sat%nLine,4+i) = (float(var_int_tmp(1:sat%nElem,1:sat%nLine,i))+interp2(i))*slope2(i)
enddo
deallocate(var_int_tmp)
if (error /= 0) then 
    print*,'ERROR: dataset EV_250_Aggr.1KM_Emissive read failed !'
    stop
endif
!call h5dclose_f(sds_id_var, error)    ! close 
endif


print*,'  ... convert radiance to tbb'
do i = 1, 6
   call dn_to_ref_rad(var_real4(:,:,i),       &
                      sat%ir_cal_coef(1,i),   &
                      sat%ir_cal_coef(2,i),   &
                      sat%ir_cal_coef(3,i),   &
                      fylat_sensor_id,        &
                      sat%rad_ir(:,:,i),      &
                      1)
   call rad_to_tbb(sat%rad_ir(:,:,i), sat%tbb_ir(:,:,i), i)
enddo
call h5dclose_f(sds_id_var, error)    ! close 
deallocate (var_real4)


if (fylat_sensor_id == 2 .or. fylat_sensor_id == 21) then
   print*,'  ... for fy3mersi_ii convert radiance (mW/m^2*sr*cm) to radiance (W/m^2*um*sr) after retrieving tbb '
   do i = 1, 6
      call ir_radcm_to_radum(sat%rad_ir(:,:,i), i)
   enddo
endif

print*,'rad1 test = ',sat%rad_ir(1000:1003,1000:1003,1)
print*,'rad2 test = ',sat%rad_ir(1000:1003,1000:1003,2)
print*,'rad3 test = ',sat%rad_ir(1000:1003,1000:1003,3)
print*,'rad4 test = ',sat%rad_ir(1000:1003,1000:1003,4)
print*,'rad5 test = ',sat%rad_ir(1000:1003,1000:1003,5)
print*,'rad6 test = ',sat%rad_ir(1000:1003,1000:1003,6)
print*,'tbb1 test = ',sat%tbb_ir(1000:1003,1000:1003,1)
print*,'tbb2 test = ',sat%tbb_ir(1000:1003,1000:1003,2)
print*,'tbb3 test = ',sat%tbb_ir(1000:1003,1000:1003,3)
print*,'tbb4 test = ',sat%tbb_ir(1000:1003,1000:1003,4)
print*,'tbb5 test = ',sat%tbb_ir(1000:1003,1000:1003,5)
print*,'tbb6 test = ',sat%tbb_ir(1000:1003,1000:1003,6)

! close file
call h5fclose_f(file_id, error)       ! close 
call h5close_f(error)      

end subroutine fylat_read_fy3_mersi_L1b_data
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~ Subroutine  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine ir_radcm_to_radum(var_inout,ir_id)

!integer(kind=4), dimension(:,:), intent(in)  :: var_in
real(kind=4), dimension(:,:), intent(inout)     :: var_inout
!real(kind=4),    dimension(:,:), intent(out) :: var_out
integer(kind=4) :: opt,ir_id,iband
integer(kind=4) :: i,j
real(kind=4)    :: wl, wn
!real(kind=4)    :: coef0, coef1, coef2, cos_ang

iband = 19+ir_id  ! for fy3mersi_ii data

wl = sat%midwave(iband)
wn = sat%midwnum(iband)

var_inout(:,:) = 1e-3*(var_inout(:,:)*wn/wl)   ! W/m^2*um*sr

!wl   = midwave(ichan) 
!wn   = 1./(w*1e-4)
!rad1 = C1/((pi*w**5)*(exp(C2/(w*T))-1))        ! W/m^2*um*sr
!rad  = 1e3*rad1*wl/wn                           ! mW/m^2*sr*cm
!rad1 = 1e-3*(rad*wn/wl)

end subroutine ir_radcm_to_radum
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~ Subroutine  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine dn_to_ref_rad(var_in,coef0,coef1,coef2,opt,var_out,ir_id)

!integer(kind=4), dimension(:,:), intent(in)  :: var_in
real(kind=4), dimension(:,:), intent(in)     :: var_in
real(kind=4),    dimension(:,:), intent(out) :: var_out
integer(kind=1) :: opt,ir_id
integer(kind=4) :: i,j
real(kind=4)    :: coef0, coef1, coef2, cos_ang

! vis band
if (opt == 1 .and. ir_id == 0) then  ! convert modis to fy3/mersi_II
   var_out(:,:) = (var_in(:,:)*1.0-coef0) *coef1  
endif

if (opt == 2 .and. ir_id == 0) then   ! convert modis to fy3/mersi_II in mersi_II format
   var_out(:,:) = (var_in(:,:)*1.0-coef0) *coef1 
endif

if (opt == 3 .and. ir_id == 0) then  ! convert viirs to fy3/mersi_II
   var_out(:,:) = coef0 + (var_in(:,:)*1.0) *coef1  + (var_in(:,:)*1.0)**2 *coef2  
endif

if (opt == 21 .and. ir_id == 0) then  ! fy3d-mersi_II 
   var_out(:,:) = (coef0 + coef1*var_in(:,:)+ coef2*var_in(:,:)**2) * 0.01 
endif

! ir band 
if (opt == 1 .and. ir_id == 1) then  ! convert modis to fy3/mersi_II
   var_out(:,:) = (var_in(:,:)*1.0-coef0) *coef1  
endif

if (opt == 2 .and. ir_id == 1) then  ! convert modis to fy3/mersi_II in mersi_II format
   var_out(:,:) = var_in(:,:)*0.01   ! ir rad = mW/m2cm-1 sr
endif

if (opt == 3 .and. ir_id == 1) then  ! convert viirs to fy3/mersi_II
   var_out(:,:) = coef0 + (var_in(:,:)*1.0) *coef1  + (var_in(:,:)*1.0)**2 *coef2  
endif

if (opt == 21 .and. ir_id == 1) then ! fy3d-mersi_II 
   var_out(:,:) = var_in(:,:)        ! ir rad = mW/m2cm-1 sr
endif


!if (opt > 10) then  ! real fy3/mersi_II
!   var_out(:,:) = coef0 + (var_in(:,:)*1.0) *coef1  + (var_in(:,:)*1.0)**2 *coef2  
!endif

if (ir_id == 0) then  !1 means is ir band

   do j = 1, sat%nLine
   do i = 1, sat%nElem
      cos_ang = cos( geo%SolarZenith(i,j) * DTOR )
      var_out(i,j) = var_out(i,j) / cos_ang
   enddo
   enddo
   
endif

end subroutine dn_to_ref_rad
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~ subroutine~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine rad_to_tbb(var_in, var_out, band)

real(kind=4), dimension(:,:), intent(in)  :: var_in
real(kind=4), dimension(:,:), intent(out) :: var_out  
real(kind=4)                              :: rad 
integer(kind=4)                           :: band, band_fylat, units
integer(kind=4)                           :: i,j,k

if (fylat_sensor_id == 1 ) then 
    units = 1 ! ...   Watts per square meter per steradian per micron
endif
if (fylat_sensor_id == 2 .or. fylat_sensor_id == 21) then 
    units = 2 ! ...   milliWatts per square meter per steradian per wavenumber
endif
band_fylat = band + 19

do i = 1, sat%nElem
do j = 1, sat%nLine

   rad = var_in(i,j)
   var_out(i,j) = fylat_planck_rad2tbb(rad, band_fylat, units)
   
enddo
enddo

end subroutine rad_to_tbb
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~ Subroutine  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine allocate_fylat_fy3mersi_geo_data()

integer :: error

allocate (geo%lon(sat%nElem, sat%nLine),           &
          geo%lat(sat%nElem, sat%nLine),           &        
          geo%SolarZenith(sat%nElem, sat%nLine),   & 
          geo%SolarAzimuth(sat%nElem, sat%nLine),  & 
          geo%SensorZenith(sat%nElem, sat%nLine),  & 
          geo%SensorAzimuth(sat%nElem, sat%nLine), &   
          geo%RelAzimuth(sat%nElem, sat%nLine),    &
          geo%GlintAngle(sat%nElem, sat%nLine),    &                    	
          geo%dem(sat%nElem, sat%nLine),           &     
          geo%lsm(sat%nElem, sat%nLine),           &   
          geo%flag(sat%nElem, sat%nLine,5),        & 
          geo%Cos_Satzen(sat%nElem, sat%nLine),    &
          geo%Cos_Solzen(sat%nElem, sat%nLine),    &        
          geo%Scatzen(sat%nElem, sat%nLine),       &       	         
          stat=error)
            
if (error /= 0) then
    print *,"(a,'Not enough memory to allocate fylat geo data arrays.')"
    stop
endif


end subroutine allocate_fylat_fy3mersi_geo_data
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~ Subroutine  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine deallocate_fylat_fy3mersi_geo_data()

integer :: error

deallocate (geo%lon,           &
            geo%lat,           &        
            geo%SolarZenith,   & 
            geo%SolarAzimuth,  & 
            geo%SensorZenith,  & 
            geo%SensorAzimuth, &  
            geo%RelAzimuth,    &
            geo%GlintAngle,    &                   	
            geo%dem,           &     
            geo%lsm,           & 
            geo%flag,          &
            geo%Cos_Satzen,    &
            geo%Cos_Solzen,    &        
            geo%Scatzen,       &       	         
            stat=error)
            
if (error /= 0) then
    print *,"(ERROR: 'Not enough memory to deallocate fylat geo data arrays.')"
    stop
endif

end subroutine deallocate_fylat_fy3mersi_geo_data
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~ Subroutine  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine allocate_fylat_fy3mersi_L1b_data()

integer :: error

allocate (sat%ref_vis(sat%nElem, sat%nLine, sat%nvis),           &
          sat%rad_ir(sat%nElem, sat%nLine,  sat%nir ),           &     
          sat%tbb_ir(sat%nElem, sat%nLine,  sat%nir ),           &    
          sat%x_nwp(sat%nElem, sat%nLine ),                      &      	   
          sat%y_nwp(sat%nElem, sat%nLine ),                      &  
          sat%zsfc(sat%nElem, sat%nLine ),                       &     
          stat=error)
    
if (error /= 0) then
    print *,"(ERROR: 'Not enough memory to allocate fylat L1b data arrays.')"
    stop
endif

if (fylat_rtm_opt > 0 ) then

allocate (sat%ivza(sat%nElem, sat%nLine ),                       &     
          sat%isfc(sat%nElem, sat%nLine ),                       &  
          sat%rad_clr38(sat%nElem, sat%nLine ),                  &
          sat%rad_clr40(sat%nElem, sat%nLine ),                  &
          sat%rad_clr73(sat%nElem, sat%nLine ),                  &
          sat%rad_clr86(sat%nElem, sat%nLine ),                  &
          sat%rad_clr11(sat%nElem, sat%nLine ),                  &
          sat%rad_clr12(sat%nElem, sat%nLine ),                  &
          sat%bt_clr38(sat%nElem, sat%nLine ),                   &
          sat%bt_clr40(sat%nElem, sat%nLine ),                   &
          sat%bt_clr73(sat%nElem, sat%nLine ),                   &
          sat%bt_clr86(sat%nElem, sat%nLine ),                   &
          sat%bt_clr11(sat%nElem, sat%nLine ),                   &
          sat%bt_clr12(sat%nElem, sat%nLine ),                   &        
          stat=error)
    
if (error /= 0) then
    print *,"(ERROR: 'Not enough memory to allocate fylat L1b rtm data arrays.')"
    stop
endif

endif

end subroutine allocate_fylat_fy3mersi_L1b_data
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~ Subroutine  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine deallocate_fylat_fy3mersi_L1b_data()

integer :: error

deallocate (sat%ref_vis,     &
            sat%rad_ir,      &   
            sat%tbb_ir,      &       
            sat%x_nwp,       &      	   
            sat%y_nwp,       &    
            sat%zsfc,        &          
            stat=error)
            
if (error /= 0) then
    print *,"(ERROR: 'Not enough memory to deallocate fylat L1b data arrays.')"
    stop
endif

if (fylat_rtm_opt > 0 ) then

deallocate (sat%ivza,                       &
            sat%isfc,                       &     	   
            sat%rad_clr38,                  &
            sat%rad_clr40,                  &
            sat%rad_clr73,                  &
            sat%rad_clr86,                  &
            sat%rad_clr11,                  &
            sat%rad_clr12,                  &
            sat%bt_clr38,                   &
            sat%bt_clr40,                   &
            sat%bt_clr73,                   &
            sat%bt_clr86,                   &
            sat%bt_clr11,                   &
            sat%bt_clr12,                   &       
            stat=error)
            
if (error /= 0) then
    print *,"(ERROR: 'Not enough memory to deallocate fylat L1b rtm data arrays.')"
    stop
endif

endif

end subroutine deallocate_fylat_fy3mersi_L1b_data
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~ Subroutine  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine fylat_write_out_cloud_mask()

!-----------------------------------------------------------------------
! !F90 write_out_arrays
!
! !Description:
!    The main program for write cloud mask array out.
!
! !Input  parameters:
!    filename        = 
!
! !Output parameters:
!    none
!-----------------------------------------------------------------------


IMPLICIT NONE

! 1. define variables
! cloud mask
integer :: error
integer (HID_T) :: file_id                 ! File identifier for cloud mask

integer (HID_T) :: sds_id_cloudmask        ! sds id for cloud mask
integer (HID_T) :: sds_id_cloudmask_qa     ! sds id for cloud mask qa

integer (HID_T) :: dsp_id_cloudmask        ! dsp id for cloud optical depth
integer (HID_T) :: dsp_id_cloudmask_qa     ! dsp id for cloud particle effective radius [unit:um]

character (LEN=10), parameter           :: dset_cloudmask_name        = "Cloud_Mask"        
character (LEN=17), parameter           :: dset_cloudmask_qa_name     = "Quality_Assurance" 

integer (HSIZE_T), dimension(3)         :: dims_cloudmask
integer (HSIZE_T), dimension(3)         :: dims_cloudmask_qa


!integer(kind=4) :: error     ! Error flag for hdf5 open
integer(kind=4), parameter :: RANK2 = 2 ! Dataset rank
integer(kind=4), parameter :: RANK3 = 3 ! Dataset rank
character,dimension(:,:,:),allocatable :: output_cloudmask
character,dimension(:,:,:),allocatable :: cloudmask_qa
character(len=20) :: satellite_name

!---------------------------------------------------
! attribute
integer(HSIZE_T),dimension(1) :: a_dims
integer(kind=4)               :: a_rank
integer(HID_T)                :: output_a_id,output_aspace_id,output_atype_id
! 2. begin program

!+++++++++++++++++++++++++++ Step 2: Write product start +++++++++++++++++++++++++++
!== initialized

!--- variables
allocate(output_cloudmask(sat%nElem,sat%nLine,cm_byte_dim),  &
         cloudmask_qa(sat%nElem,sat%nLine,cm_qa_dim))

!------         
dims_cloudmask(1) = sat%nElem
dims_cloudmask(2) = sat%nLine
dims_cloudmask(3) = cm_byte_dim
dims_cloudmask_qa(1) = sat%nElem
dims_cloudmask_qa(2) = sat%nLine
dims_cloudmask_qa(3) = cm_qa_dim

!== 2.1. write start
print*,'    ... fylat write out fy3/MERSI_II Cloud Mask HDF5 product at first !!! '
	 
call h5open_f(error)

call h5fcreate_f(trim(fy3_mersi_CLM_data), H5F_ACC_TRUNC_F, file_id, error)

!------------------------------------
! --- Write Cloud Mask
output_cloudmask = char(cm_bitarray)!cm_bitarray, cm_qa_bitarray
call h5screate_simple_f(RANK3, dims_cloudmask, dsp_id_cloudmask, error)
call h5dcreate_f(file_id, dset_cloudmask_name, H5T_STD_U8LE, dsp_id_cloudmask, sds_id_cloudmask, error)
call h5dwrite_f (sds_id_cloudmask, H5T_STD_U8LE, output_cloudmask, dims_cloudmask, error)

! ... write the attributes
      ! .. units attribute
      a_rank = 1
      a_dims = 1
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cloudmask,'units',H5T_NATIVE_CHARACTER, output_aspace_id, output_a_id, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. valid_range attribute
      a_dims = 2
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cloudmask,'valid_range',H5T_NATIVE_INTEGER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_INTEGER, (/0,255/), a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. _FillValue attibute
      a_dims = 1 
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cloudmask,'_FillValue',H5T_NATIVE_INTEGER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_INTEGER, 0 , a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. Intercept attribute
      a_dims = 1 
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cloudmask,'Intercept',H5T_NATIVE_REAL, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_REAL, 0.0 , a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)     
      ! .. Slope attribute
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cloudmask,'Slope',H5T_NATIVE_REAL, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_REAL, 1.0 , a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. long_name
      a_dims = 1
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cloudmask,'long_name',H5T_NATIVE_CHARACTER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_CHARACTER, 'fylat MERSI_II Cloud Mask', a_dims, error)
      !call h5awrite_f(output_a_id, H5T_NATIVE_SHORT, 'fylat MERSI_II Cloud Mask', a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. band_name
      a_dims = 1 
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cloudmask,'band_name',H5T_NATIVE_CHARACTER, output_aspace_id, output_a_id, error)
!      call h5awrite_f(output_a_id, H5T_NATIVE_CHARACTER, 'fylat MERSI Cloud Mask', a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      
      ! Close output dataset ,dataspace, file

call h5dclose_f (sds_id_cloudmask, error)      
call h5sclose_f (dsp_id_cloudmask, error)

!------------------------------------
! --- Write Quality Assurance
cloudmask_qa = char(cm_qa_bitarray) !cm_bitarray, cm_qa_bitarray
call h5screate_simple_f(RANK3, dims_cloudmask_qa, dsp_id_cloudmask_qa, error)
call h5dcreate_f(file_id, dset_cloudmask_qa_name, H5T_STD_U8LE, dsp_id_cloudmask_qa, sds_id_cloudmask_qa, error)
call h5dwrite_f (sds_id_cloudmask_qa, H5T_STD_U8LE, cloudmask_qa, dims_cloudmask_qa, error)

! ... write the attributes
      ! .. units attribute
      a_rank = 1
      a_dims = 1
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cloudmask_qa,'units',H5T_NATIVE_CHARACTER, output_aspace_id, output_a_id, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. valid_range attribute
      a_dims = 2
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cloudmask_qa,'valid_range',H5T_NATIVE_INTEGER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_INTEGER, (/0,255/), a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. _FillValue attibute
      a_dims = 1 
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cloudmask_qa,'_FillValue',H5T_NATIVE_INTEGER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_INTEGER, 0 , a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. Intercept attribute
      a_dims = 1 
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cloudmask_qa,'Intercept',H5T_NATIVE_REAL, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_REAL, 0.0 , a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)     
      ! .. Slope attribute
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cloudmask_qa,'Slope',H5T_NATIVE_REAL, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_REAL, 1.0 , a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. long_name
      a_dims = 25
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cloudmask_qa,'long_name',H5T_NATIVE_CHARACTER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_CHARACTER, 'fylat MERSI_II Cloud Mask Quality Assurance', a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. band_name
      a_dims = 1
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cloudmask_qa,'band_name',H5T_NATIVE_CHARACTER, output_aspace_id, output_a_id, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)

call h5dclose_f (sds_id_cloudmask_qa, error)      
call h5sclose_f (dsp_id_cloudmask_qa, error)
      
! close file
call h5fclose_f(file_id, error)

!----------------------------------------------------

call h5close_f(error)

deallocate(output_cloudmask, cloudmask_qa)

end subroutine fylat_write_out_cloud_mask
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!~~~~~~~~~~~~~~~~ Subroutine  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine fylat_write_out_cloud_amount()

!-----------------------------------------------------------------------
! !F90 write_out_arrays
!
! !Description:
!    The main program for write cloud mask array out.
!
! !Input  parameters:
!    filename        = 
!
! !Output parameters:
!    none
!-----------------------------------------------------------------------


IMPLICIT NONE

! 1. define variables
! cloud mask
integer :: error
integer (HID_T) :: file_id                 ! File identifier for cloud mask

integer (HID_T) :: sds_id_cloudamount        ! sds id for cloud mask
integer (HID_T) :: sds_id_cloudamount_qa     ! sds id for cloud mask qa
integer (HID_T) :: sds_id_lon5km             ! sds id for cloud mask
integer (HID_T) :: sds_id_lat5km             ! sds id for cloud mask qa

integer (HID_T) :: dsp_id_cloudamount        ! dsp id for cloud optical depth
integer (HID_T) :: dsp_id_cloudamount_qa     ! dsp id for cloud particle effective radius [unit:um]
integer (HID_T) :: dsp_id_lon5km             ! dsp id for cloud optical depth
integer (HID_T) :: dsp_id_lat5km             ! dsp id for cloud particle effective radius [unit:um]

character (LEN=28), parameter           :: dset_cloudamount_name        = "5-min granule Cloud Fraction"        
character (LEN=37), parameter           :: dset_cloudamount_qa_name     = "5-min granule Cloud Fraction QA_flags" 
character (LEN= 8), parameter           :: dset_lon5km_name             = "Latitude"        
character (LEN= 9), parameter           :: dset_lat5km_name             = "Longitude" 

integer (HSIZE_T), dimension(2)         :: dims_cloudamount
integer (HSIZE_T), dimension(2)         :: dims_cloudamount_qa
integer (HSIZE_T), dimension(2)         :: dims_lon5km
integer (HSIZE_T), dimension(2)         :: dims_lat5km

!integer(kind=4) :: error     ! Error flag for hdf5 open
integer(kind=4), parameter :: RANK2 = 2 ! Dataset rank
integer(kind=4), parameter :: RANK3 = 3 ! Dataset rank
integer(kind=4),dimension(:,:),allocatable :: output_cloudamount
character,dimension(:,:),allocatable :: cloudamount_qa
character(len=20) :: satellite_name

!---------------------------------------------------
! attribute
integer(HSIZE_T),dimension(1) :: a_dims
integer(kind=4)               :: a_rank
integer(HID_T)                :: output_a_id,output_aspace_id,output_atype_id
! 2. begin program

!+++++++++++++++++++++++++++ Step 2: Write product start +++++++++++++++++++++++++++
!== initialized

!--- variables
allocate(output_cloudamount(ix_5km,iy_5km),  &
         cloudamount_qa(ix_5km,iy_5km))

!------         
dims_cloudamount(1) = ix_5km
dims_cloudamount(2) = iy_5km
dims_cloudamount_qa(1) = ix_5km
dims_cloudamount_qa(2) = iy_5km
dims_lon5km(1) = ix_5km
dims_lon5km(2) = iy_5km
dims_lat5km(1) = ix_5km
dims_lat5km(2) = iy_5km

!== 2.1. write start
print*,'    ... fylat write out fy3/MERSI_II Cloud Amount HDF5 product !!! '
	 
call h5open_f(error)

call h5fcreate_f(trim(fy3_mersi_CLA_data), H5F_ACC_TRUNC_F, file_id, error)

!------------------------------------
! --- Write Cloud Mask
!call h5dcreate_f(group_id, "EV_1KM_RefSB", H5T_STD_U16LE, dsp_id, sds_id, error)
!call h5dwrite_f (sds_id, H5T_NATIVE_INTEGER, rad(:,:,5:19), dims_sp32, error)
output_cloudamount = cloud_amount !cm_bitarray, cm_qa_bitarray
call h5screate_simple_f(RANK2, dims_cloudamount, dsp_id_cloudamount, error)
call h5dcreate_f(file_id, dset_cloudamount_name, H5T_STD_I16LE, dsp_id_cloudamount, sds_id_cloudamount, error)
call h5dwrite_f (sds_id_cloudamount, H5T_NATIVE_INTEGER, output_cloudamount, dims_cloudamount, error)

! ... write the attributes
      ! .. units attribute
      a_rank = 1
      a_dims = 1
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cloudamount,'units',H5T_NATIVE_CHARACTER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_CHARACTER, 'none', a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. valid_range attribute
      a_dims = 2
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cloudamount,'valid_range',H5T_NATIVE_INTEGER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_INTEGER, (/0,100/), a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. _FillValue attibute
      a_dims = 1 
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cloudamount,'_FillValue',H5T_NATIVE_INTEGER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_INTEGER, -999 , a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. Intercept attribute
      a_dims = 1 
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cloudamount,'Intercept',H5T_NATIVE_REAL, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_REAL, 0.0 , a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)     
      ! .. Slope attribute
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cloudamount,'Slope',H5T_NATIVE_REAL, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_REAL, 1.0 , a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. long_name
      a_dims = 1
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cloudamount,'long_name',H5T_NATIVE_CHARACTER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_CHARACTER, '5-min granule Cloud Fraction', a_dims, error)
      !call h5awrite_f(output_a_id, H5T_NATIVE_SHORT, 'fylat MERSI_II Cloud Mask', a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. band_name
      a_dims = 1 
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cloudamount,'band_name',H5T_NATIVE_CHARACTER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_CHARACTER, 'NANA', a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      
      ! Close output dataset ,dataspace, file

call h5dclose_f (sds_id_cloudamount, error)      
call h5sclose_f (dsp_id_cloudamount, error)

!------------------------------------
! --- Write Quality Assurance
cloudamount_qa = char(cloud_amount_qa) !cm_bitarray, cm_qa_bitarray
call h5screate_simple_f(RANK2, dims_cloudamount_qa, dsp_id_cloudamount_qa, error)
call h5dcreate_f(file_id, dset_cloudamount_qa_name, H5T_STD_U8LE, dsp_id_cloudamount_qa, sds_id_cloudamount_qa, error)
call h5dwrite_f (sds_id_cloudamount_qa, H5T_STD_U8LE, cloudamount_qa, dims_cloudamount_qa, error)

! ... write the attributes
      ! .. units attribute
      a_rank = 1
      a_dims = 1
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cloudamount_qa,'units',H5T_NATIVE_CHARACTER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_CHARACTER, 'none', a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. valid_range attribute
      a_dims = 2
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cloudamount_qa,'valid_range',H5T_NATIVE_INTEGER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_INTEGER, (/0,2/), a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. _FillValue attibute
      a_dims = 1 
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cloudamount_qa,'_FillValue',H5T_NATIVE_INTEGER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_INTEGER, 0 , a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. Intercept attribute
      a_dims = 1 
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cloudamount_qa,'Intercept',H5T_NATIVE_REAL, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_REAL, 0.0 , a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)     
      ! .. Slope attribute
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cloudamount_qa,'Slope',H5T_NATIVE_REAL, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_REAL, 1.0 , a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. long_name
      a_dims = 37
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cloudamount_qa,'long_name',H5T_NATIVE_CHARACTER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_CHARACTER, '5-min granule Cloud Fraction QA_flags', a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. band_name
      a_dims = 1
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cloudamount_qa,'band_name',H5T_NATIVE_CHARACTER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_CHARACTER, 'NANA', a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)

call h5dclose_f (sds_id_cloudamount_qa, error)      
call h5sclose_f (dsp_id_cloudamount_qa, error)

!------------------------------------
! --- Write lon and lat
call h5screate_simple_f(RANK2, dims_lon5km, dsp_id_lon5km, error)
call h5dcreate_f(file_id, dset_lon5km_name, H5T_STD_I16LE, dsp_id_lon5km, sds_id_lon5km, error)
call h5dwrite_f (sds_id_lon5km, H5T_NATIVE_INTEGER, lon_5km, dims_lon5km, error)
! ... write the attributes
      ! .. units attribute
      a_rank = 1
      a_dims = 1
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_lon5km,'units',H5T_NATIVE_CHARACTER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_CHARACTER, 'degree', a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. valid_range attribute
      a_dims = 2
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_lon5km,'valid_range',H5T_NATIVE_INTEGER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_INTEGER, (/-18000,18000/), a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. _FillValue attibute
      a_dims = 1 
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_lon5km,'_FillValue',H5T_NATIVE_INTEGER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_INTEGER, 65535 , a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. Intercept attribute
      a_dims = 1 
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_lon5km,'Intercept',H5T_NATIVE_REAL, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_REAL, 0.0 , a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)     
      ! .. Slope attribute
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_lon5km,'Slope',H5T_NATIVE_REAL, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_REAL, 0.01, a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. long_name
      a_dims = 1
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cloudamount,'long_name',H5T_NATIVE_CHARACTER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_CHARACTER, 'Longitude', a_dims, error)
      !call h5awrite_f(output_a_id, H5T_NATIVE_SHORT, 'fylat MERSI_II Cloud Mask', a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. band_name
      a_dims = 1 
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_lon5km,'band_name',H5T_NATIVE_CHARACTER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_CHARACTER, ' ', a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      
call h5dclose_f (sds_id_lon5km, error)      
call h5sclose_f (dsp_id_lon5km, error)

call h5screate_simple_f(RANK2, dims_lat5km, dsp_id_lat5km, error)
call h5dcreate_f(file_id, dset_lat5km_name, H5T_STD_I16LE, dsp_id_lat5km, sds_id_lat5km, error)
call h5dwrite_f (sds_id_lat5km, H5T_NATIVE_INTEGER, lat_5km, dims_lat5km, error)
! ... write the attributes
      ! .. units attribute
      a_rank = 1
      a_dims = 1
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_lat5km,'units',H5T_NATIVE_CHARACTER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_CHARACTER, 'degree', a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. valid_range attribute
      a_dims = 2
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cloudamount,'valid_range',H5T_NATIVE_INTEGER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_INTEGER, (/-9000,9000/), a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. _FillValue attibute
      a_dims = 1 
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_lat5km,'_FillValue',H5T_NATIVE_INTEGER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_INTEGER, 65535 , a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. Intercept attribute
      a_dims = 1 
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_lat5km,'Intercept',H5T_NATIVE_REAL, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_REAL, 0.0 , a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)     
      ! .. Slope attribute
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_lat5km,'Slope',H5T_NATIVE_REAL, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_REAL, 0.01, a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. long_name
      a_dims = 1
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_lat5km,'long_name',H5T_NATIVE_CHARACTER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_CHARACTER, 'Latitude', a_dims, error)
      !call h5awrite_f(output_a_id, H5T_NATIVE_SHORT, 'fylat MERSI_II Cloud Mask', a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. band_name
      a_dims = 1 
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_lat5km,'band_name',H5T_NATIVE_CHARACTER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_CHARACTER, ' ', a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      
call h5dclose_f (sds_id_lat5km, error)      
call h5sclose_f (dsp_id_lat5km, error)
      
! close file
call h5fclose_f(file_id, error)

!----------------------------------------------------

call h5close_f(error)

deallocate(output_cloudamount, cloudamount_qa)

end subroutine fylat_write_out_cloud_amount
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!~~~~~~~~~~~~~~~~ Subroutine  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine fylat_write_out_cloud_phase()

!-----------------------------------------------------------------------
! !F90 write_out_arrays
!
! !Description:
!    The main program for write cloud mask array out.
!
! !Input  parameters:
!    filename        = 
!
! !Output parameters:
!    none
!-----------------------------------------------------------------------


IMPLICIT NONE

! 1. define variables
! cloud mask
integer :: error
integer (HID_T) :: file_id                 ! File identifier for cloud mask
  
integer (HID_T) :: sds_id_cloudtype        ! sds id for cloud mask
integer (HID_T) :: sds_id_cloudphase       ! sds id for cloud mask qa
integer (HID_T) :: sds_id_Cldphase_Qf      ! sds id for cloud mask
integer (HID_T) :: sds_id_Cldphase_Qpi     ! sds id for cloud mask qa
!--others
integer (HID_T) :: sds_id_Cldtype_Tmpy     ! sds id for cloud mask qa
integer (HID_T) :: sds_id_Xgrad_Emiss14     ! sds id for cloud mask qa
integer (HID_T) :: sds_id_Ygrad_Emiss14     ! sds id for cloud mask qa
integer (HID_T) :: sds_id_Num_Steps_Gradient        ! sds id for cloud mask qa
integer (HID_T) :: sds_id_Cldphase_Lrc_Mask         ! sds id for cloud mask qa
integer (HID_T) :: sds_id_Emiss_Chn10_Tot           ! sds id for cloud mask qa
integer (HID_T) :: sds_id_Emiss_Chn10_Tot_Lrc       ! sds id for cloud mask qa
integer (HID_T) :: sds_id_Emiss_Chn14_Tot           ! sds id for cloud mask qa
integer (HID_T) :: sds_id_Emiss_Chn14_Tot_Multi     ! sds id for cloud mask qa
integer (HID_T) :: sds_id_Cldbeta7311_Tot     ! sds id for cloud mask qa
integer (HID_T) :: sds_id_Cldbeta7311_Tot_Multi           ! sds id for cloud mask qa
integer (HID_T) :: sds_id_Cldbeta1112_Tot       ! sds id for cloud mask qa
integer (HID_T) :: sds_id_Cldbeta1112_Tot_Lrc            ! sds id for cloud mask qa
integer (HID_T) :: sds_id_Cldbeta1112_Tot_Multi    ! sds id for cloud mask qa
integer (HID_T) :: sds_id_Cldbeta1112_Opaque     ! sds id for cloud mask qa
integer (HID_T) :: sds_id_Cldbeta1112_Opaque_Multi            ! sds id for cloud mask qa
integer (HID_T) :: sds_id_Cldbeta8511_Tot       ! sds id for cloud mask qa
integer (HID_T) :: sds_id_Cldbeta8511_Tot_Multi            ! sds id for cloud mask qa
integer (HID_T) :: sds_id_Cldbeta8511_Opaque     ! sds id for cloud mask qa
integer (HID_T) :: sds_id_Cldbeta8511_Opaque_Lrc       ! sds id for cloud mask qa
integer (HID_T) :: sds_id_Cldbeta8511_Opaque_Multi            ! sds id for cloud mask qa
integer (HID_T) :: sds_id_Opaque_Cld_Temp_Chn10      ! sds id for cloud mask qa
integer (HID_T) :: sds_id_Opaque_Cld_Temp_Chn14      ! sds id for cloud mask qa  


!-----
integer (HID_T) :: dsp_id_cloudtype        ! dsp id for cloud optical depth
integer (HID_T) :: dsp_id_cloudphase       ! dsp id for cloud particle effective radius [unit:um]
integer (HID_T) :: dsp_id_Cldphase_Qf      ! dsp id for cloud optical depth
integer (HID_T) :: dsp_id_Cldphase_Qpi     ! dsp id for cloud particle effective radius [unit:um]

integer (HID_T) :: dsp_id_Cldtype_Tmpy          ! dsp id for cloud particle effective radius [unit:um]
integer (HID_T) :: dsp_id_Xgrad_Emiss14         ! dsp id for cloud particle effective radius [unit:um]
integer (HID_T) :: dsp_id_Ygrad_Emiss14         ! dsp id for cloud particle effective radius [unit:um]
integer (HID_T) :: dsp_id_Num_Steps_Gradient    ! dsp id for cloud particle effective radius [unit:um]
integer (HID_T) :: dsp_id_Cldphase_Lrc_Mask     ! dsp id for cloud particle effective radius [unit:um]
integer (HID_T) :: dsp_id_Emiss_Chn10_Tot           ! sds id for cloud mask qa
integer (HID_T) :: dsp_id_Emiss_Chn10_Tot_Lrc       ! sds id for cloud mask qa
integer (HID_T) :: dsp_id_Emiss_Chn14_Tot           ! sds id for cloud mask qa
integer (HID_T) :: dsp_id_Emiss_Chn14_Tot_Multi     ! sds id for cloud mask qa
integer (HID_T) :: dsp_id_Cldbeta7311_Tot     ! sds id for cloud mask qa
integer (HID_T) :: dsp_id_Cldbeta7311_Tot_Multi           ! sds id for cloud mask qa
integer (HID_T) :: dsp_id_Cldbeta1112_Tot       ! sds id for cloud mask qa
integer (HID_T) :: dsp_id_Cldbeta1112_Tot_Lrc            ! sds id for cloud mask qa
integer (HID_T) :: dsp_id_Cldbeta1112_Tot_Multi    ! sds id for cloud mask qa
integer (HID_T) :: dsp_id_Cldbeta1112_Opaque     ! sds id for cloud mask qa
integer (HID_T) :: dsp_id_Cldbeta1112_Opaque_Multi            ! sds id for cloud mask qa
integer (HID_T) :: dsp_id_Cldbeta8511_Tot       ! sds id for cloud mask qa
integer (HID_T) :: dsp_id_Cldbeta8511_Tot_Multi            ! sds id for cloud mask qa
integer (HID_T) :: dsp_id_Cldbeta8511_Opaque     ! sds id for cloud mask qa
integer (HID_T) :: dsp_id_Cldbeta8511_Opaque_Lrc       ! sds id for cloud mask qa
integer (HID_T) :: dsp_id_Cldbeta8511_Opaque_Multi            ! sds id for cloud mask qa
integer (HID_T) :: dsp_id_Opaque_Cld_Temp_Chn10      ! sds id for cloud mask qa
integer (HID_T) :: dsp_id_Opaque_Cld_Temp_Chn14      ! sds id for cloud mask qa  

!-----
character (LEN=10), parameter           :: dset_cloudtype_name      = "Cloud_Type"        
character (LEN=11), parameter           :: dset_cloudphase_name     = "Cloud_Phase" 
character (LEN=11), parameter           :: dset_Cldphase_Qf_name    = "Cldphase_Qf"        
character (LEN=12), parameter           :: dset_Cldphase_Qpi_name   = "Cldphase_Qpi" 
! others
character (LEN=11), parameter           :: dset_Cldtype_Tmpy_name        = "Cldtype_Tmpy" 
character (LEN=13), parameter           :: dset_Xgrad_Emiss14_name       = "Xgrad_Emiss14" 
character (LEN=13), parameter           :: dset_Ygrad_Emiss14_name       = "Ygrad_Emiss14" 
character (LEN=18), parameter           :: dset_Num_Steps_Gradient_name  = "Num_Steps_Gradient" 
character (LEN=17), parameter           :: dset_Cldphase_Lrc_Mask_name   = "Cldphase_Lrc_Mask" 

character (LEN=15), parameter           :: dset_Emiss_Chn10_Tot_name          = "Emiss_Chn10_Tot" 
character (LEN=19), parameter           :: dset_Emiss_Chn10_Tot_Lrc_name      = "Emiss_Chn10_Tot_Lrc" 
character (LEN=15), parameter           :: dset_Emiss_Chn14_Tot_name          = "Emiss_Chn14_Tot" 
character (LEN=19), parameter           :: dset_Emiss_Chn14_Tot_Multi_name    = "Emiss_Chn14_Tot_Multi" 
character (LEN=15), parameter           :: dset_Cldbeta7311_Tot_name          = "Cldbeta7311_Tot" 
character (LEN=21), parameter           :: dset_Cldbeta7311_Tot_Multi_name    = "Cldbeta7311_Tot_Multi" 
character (LEN=15), parameter           :: dset_Cldbeta1112_Tot_name          = "Cldbeta1112_Tot" 
character (LEN=19), parameter           :: dset_Cldbeta1112_Tot_Lrc_name      = "Cldbeta1112_Tot_Lrc" 
character (LEN=32), parameter           :: dset_Cldbeta1112_Tot_Multi_name    = "Cldbeta1112_Tot_Multi" 
character (LEN=18), parameter           :: dset_Cldbeta1112_Opaque_name       = "Cldbeta1112_Opaque" 
character (LEN=24), parameter           :: dset_Cldbeta1112_Opaque_Multi_name = "Cldbeta1112_Opaque_Multi" 
character (LEN=15), parameter           :: dset_Cldbeta8511_Tot_name          = "Cldbeta8511_Tot" 
character (LEN=19), parameter           :: dset_Cldbeta8511_Tot_Multi_name    = "Cldbeta8511_Tot_Multi" 
character (LEN=18), parameter           :: dset_Cldbeta8511_Opaque_name       = "Cldbeta8511_Opaque" 
character (LEN=22), parameter           :: dset_Cldbeta8511_Opaque_Lrc_name   = "Cldbeta8511_Opaque_Lrc" 
character (LEN=24), parameter           :: dset_Cldbeta8511_Opaque_Multi_name = "Cldbeta8511_Opaque_Multi"
character (LEN=21), parameter           :: dset_Opaque_Cld_Temp_Chn10_name    = "Opaque_Cld_Temp_Chn10" 
character (LEN=21), parameter           :: dset_Opaque_Cld_Temp_Chn14_name    = "Opaque_Cld_Temp_Chn14"

integer (HSIZE_T), dimension(2)         :: dims_cloudtype
integer (HSIZE_T), dimension(2)         :: dims_cloudphase
integer (HSIZE_T), dimension(3)         :: dims_Cldphase_Qf
integer (HSIZE_T), dimension(3)         :: dims_Cldphase_Qpi

!integer(kind=4) :: error     ! Error flag for hdf5 open
integer(kind=4), parameter :: RANK2 = 2 ! Dataset rank
integer(kind=4), parameter :: RANK3 = 3 ! Dataset rank
character,dimension(:,:),allocatable   :: output_cloudp
character,dimension(:,:,:),allocatable :: Cldphase_Qf
character,dimension(:,:,:),allocatable :: Cldphase_Qpi
character(len=20) :: satellite_name

!---------------------------------------------------
! attribute
integer(HSIZE_T),dimension(1) :: a_dims
integer(kind=4)               :: a_rank
integer(HID_T)                :: output_a_id,output_aspace_id,output_atype_id
! 2. begin program

!+++++++++++++++++++++++++++ Step 2: Write product start +++++++++++++++++++++++++++
!== initialized

!--- variables
allocate(output_cloudp(sat%nElem,sat%nLine),    &
         Cldphase_Qf(6,sat%nElem,sat%nLine),    &
         Cldphase_Qpi(20,sat%nElem,sat%nLine) )

!------         
dims_cloudtype(1) = sat%nElem
dims_cloudtype(2) = sat%nLine
dims_cloudphase(1) = sat%nElem
dims_cloudphase(2) = sat%nLine
dims_Cldphase_Qf(1) = 6
dims_Cldphase_Qf(2) = sat%nElem
dims_Cldphase_Qf(3) = sat%nLine
dims_Cldphase_Qpi(1) = 20
dims_Cldphase_Qpi(2) = sat%nElem
dims_Cldphase_Qpi(3) = sat%nLine

!== 2.1. write start
print*,'    ... fylat write out fy3/MERSI_II Cloud Phase HDF5 product !!! '
	 
call h5open_f(error)

call h5fcreate_f(trim(fy3_mersi_CLP_data), H5F_ACC_TRUNC_F, file_id, error)

!------------------------------------
! --- Write Cloud Phase and Type
!call h5dcreate_f(group_id, "EV_1KM_RefSB", H5T_STD_U16LE, dsp_id, sds_id, error)
!call h5dwrite_f (sds_id, H5T_NATIVE_INTEGER, rad(:,:,5:19), dims_sp32, error)
output_cloudp = char(clp%Cldphase)!cm_bitarray, cm_qa_bitarray
call h5screate_simple_f(RANK2, dims_cloudphase, dsp_id_cloudphase, error)
call h5dcreate_f(file_id, dset_cloudphase_name, H5T_STD_U8LE, dsp_id_cloudphase, sds_id_cloudphase, error)
call h5dwrite_f (sds_id_cloudphase, H5T_STD_U8LE, output_cloudp, dims_cloudphase, error)

! ... write the attributes
      ! .. units attribute
      a_rank = 1
      a_dims = 1
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cloudphase,'units',H5T_NATIVE_CHARACTER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_CHARACTER, 'none', a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. valid_range attribute
      a_dims = 2
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cloudphase,'valid_range',H5T_NATIVE_INTEGER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_INTEGER, (/0,4/), a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. _FillValue attibute
      a_dims = 1 
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cloudphase,'_FillValue',H5T_NATIVE_INTEGER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_INTEGER, 5 , a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. Intercept attribute
      a_dims = 1 
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cloudphase,'Intercept',H5T_NATIVE_REAL, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_REAL, 0.0 , a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)     
      ! .. Slope attribute
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cloudphase,'Slope',H5T_NATIVE_REAL, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_REAL, 1.0 , a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. long_name
      a_dims = 1
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cloudphase,'long_name',H5T_NATIVE_CHARACTER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_CHARACTER, '5-min granule Cloud Phase', a_dims, error)
      !call h5awrite_f(output_a_id, H5T_NATIVE_SHORT, 'fylat MERSI_II Cloud Mask', a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. band_name
      a_dims = 1 
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cloudphase,'band_name',H5T_NATIVE_CHARACTER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_CHARACTER, 'NANA', a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      
      ! Close output dataset ,dataspace, file

call h5dclose_f (sds_id_cloudphase, error)      
call h5sclose_f (dsp_id_cloudphase, error)

! cloud type
output_cloudp = char(clp%Cldtype)!cm_bitarray, cm_qa_bitarray
call h5screate_simple_f(RANK2, dims_cloudtype, dsp_id_cloudtype, error)
call h5dcreate_f(file_id, dset_cloudtype_name, H5T_STD_U8LE, dsp_id_cloudtype, sds_id_cloudtype, error)
call h5dwrite_f (sds_id_cloudtype, H5T_STD_U8LE, output_cloudp, dims_cloudtype, error)

! ... write the attributes
      ! .. units attribute
      a_rank = 1
      a_dims = 1
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cloudtype,'units',H5T_NATIVE_CHARACTER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_CHARACTER, 'none', a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. valid_range attribute
      a_dims = 2
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cloudtype,'valid_range',H5T_NATIVE_INTEGER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_INTEGER, (/0,8/), a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. _FillValue attibute
      a_dims = 1 
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cloudtype,'_FillValue',H5T_NATIVE_INTEGER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_INTEGER, 9 , a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. Intercept attribute
      a_dims = 1 
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cloudtype,'Intercept',H5T_NATIVE_REAL, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_REAL, 0.0 , a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)     
      ! .. Slope attribute
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cloudtype,'Slope',H5T_NATIVE_REAL, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_REAL, 1.0 , a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. long_name
      a_dims = 1
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cloudphase,'long_name',H5T_NATIVE_CHARACTER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_CHARACTER, '5-min granule Cloud Phase', a_dims, error)
      !call h5awrite_f(output_a_id, H5T_NATIVE_SHORT, 'fylat MERSI_II Cloud Mask', a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. band_name
      a_dims = 1 
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cloudtype,'band_name',H5T_NATIVE_CHARACTER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_CHARACTER, 'NANA', a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      
      ! Close output dataset ,dataspace, file

call h5dclose_f (sds_id_cloudtype, error)      
call h5sclose_f (dsp_id_cloudtype, error)

!------------------------------------
! --- Write Quality Assurance
Cldphase_Qf = char(clp%Cldphase_Qf) !cm_bitarray, cm_qa_bitarray
call h5screate_simple_f(RANK3, dims_Cldphase_Qf, dsp_id_Cldphase_Qf, error)
call h5dcreate_f(file_id, dset_Cldphase_Qf_name, H5T_STD_U8LE, dsp_id_Cldphase_Qf, sds_id_Cldphase_Qf, error)
call h5dwrite_f (sds_id_Cldphase_Qf, H5T_STD_U8LE, Cldphase_Qf, dims_Cldphase_Qf, error)

! ... write the attributes
      ! .. units attribute
      a_rank = 1
      a_dims = 1
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_Cldphase_Qf,'units',H5T_NATIVE_CHARACTER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_CHARACTER, 'none', a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. valid_range attribute
      a_dims = 2
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_Cldphase_Qf,'valid_range',H5T_NATIVE_INTEGER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_INTEGER, (/0,2/), a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. _FillValue attibute
      a_dims = 1 
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_Cldphase_Qf,'_FillValue',H5T_NATIVE_INTEGER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_INTEGER, 0 , a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. Intercept attribute
      a_dims = 1 
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_Cldphase_Qf,'Intercept',H5T_NATIVE_REAL, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_REAL, 0.0 , a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)     
      ! .. Slope attribute
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_Cldphase_Qf,'Slope',H5T_NATIVE_REAL, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_REAL, 1.0 , a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. long_name
      a_dims = 37
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_Cldphase_Qf,'long_name',H5T_NATIVE_CHARACTER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_CHARACTER, '5-min granule Cloud Fraction QA_flags', a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. band_name
      a_dims = 1
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_Cldphase_Qf,'band_name',H5T_NATIVE_CHARACTER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_CHARACTER, 'NANA', a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)

call h5dclose_f (sds_id_Cldphase_Qf, error)      
call h5sclose_f (dsp_id_Cldphase_Qf, error)

!
Cldphase_Qpi = char(clp%Cldphase_Qpi) !cm_bitarray, cm_qa_bitarray
call h5screate_simple_f(RANK3, dims_Cldphase_Qpi, dsp_id_Cldphase_Qpi, error)
call h5dcreate_f(file_id, dset_Cldphase_Qpi_name, H5T_STD_U8LE, dsp_id_Cldphase_Qpi, sds_id_Cldphase_Qpi, error)
call h5dwrite_f (sds_id_Cldphase_Qpi, H5T_STD_U8LE, Cldphase_Qpi, dims_Cldphase_Qpi, error)

! ... write the attributes
      ! .. units attribute
      a_rank = 1
      a_dims = 1
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_Cldphase_Qpi,'units',H5T_NATIVE_CHARACTER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_CHARACTER, 'none', a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. valid_range attribute
      a_dims = 2
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_Cldphase_Qpi,'valid_range',H5T_NATIVE_INTEGER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_INTEGER, (/0,2/), a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. _FillValue attibute
      a_dims = 1 
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_Cldphase_Qpi,'_FillValue',H5T_NATIVE_INTEGER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_INTEGER, 0 , a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. Intercept attribute
      a_dims = 1 
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_Cldphase_Qpi,'Intercept',H5T_NATIVE_REAL, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_REAL, 0.0 , a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)     
      ! .. Slope attribute
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_Cldphase_Qpi,'Slope',H5T_NATIVE_REAL, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_REAL, 1.0 , a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. long_name
      a_dims = 37
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_Cldphase_Qpi,'long_name',H5T_NATIVE_CHARACTER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_CHARACTER, '5-min granule Cloud Fraction QA_flags', a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. band_name
      a_dims = 1
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_Cldphase_Qpi,'band_name',H5T_NATIVE_CHARACTER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_CHARACTER, 'NANA', a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)

call h5dclose_f (sds_id_Cldphase_Qpi, error)      
call h5sclose_f (dsp_id_Cldphase_Qpi, error)

!------------------------------------
! --- Write others out

if (fylat_alg_opt%cloudphase_index > 1) then

output_cloudp = char(clp%Cldtype_Tmpy) !cm_bitarray, cm_qa_bitarray
call h5screate_simple_f(RANK2, dims_cloudphase, dsp_id_Cldtype_Tmpy, error)
call h5dcreate_f(file_id, dset_Cldtype_Tmpy_name, H5T_STD_U8LE, dsp_id_Cldtype_Tmpy, sds_id_Cldtype_Tmpy, error)
call h5dwrite_f (sds_id_Cldtype_Tmpy, H5T_STD_U8LE, output_cloudp, dims_cloudphase, error)
call h5dclose_f (sds_id_Cldtype_Tmpy, error)      
call h5sclose_f (dsp_id_Cldtype_Tmpy, error)

output_cloudp = char(clp%Cldphase_Lrc_Mask) !cm_bitarray, cm_qa_bitarray
call h5screate_simple_f(RANK2, dims_cloudphase, dsp_id_Cldphase_Lrc_Mask, error)
call h5dcreate_f(file_id, dset_Cldphase_Lrc_Mask_name, H5T_STD_U8LE, dsp_id_Cldphase_Lrc_Mask, sds_id_Cldphase_Lrc_Mask, error)
call h5dwrite_f (sds_id_Cldphase_Lrc_Mask, H5T_STD_U8LE, output_cloudp, dims_cloudphase, error)
call h5dclose_f (sds_id_Cldphase_Lrc_Mask, error)      
call h5sclose_f (dsp_id_Cldphase_Lrc_Mask, error)

call h5screate_simple_f(RANK2, dims_cloudphase, dsp_id_Xgrad_Emiss14, error)
call h5dcreate_f(file_id, dset_Xgrad_Emiss14_name, H5T_NATIVE_INTEGER, dsp_id_Xgrad_Emiss14, sds_id_Xgrad_Emiss14, error)
call h5dwrite_f (sds_id_Xgrad_Emiss14, H5T_NATIVE_INTEGER, clp%Xgrad_Emiss14, dims_cloudphase, error)
call h5dclose_f (sds_id_Xgrad_Emiss14, error)      
call h5sclose_f (dsp_id_Xgrad_Emiss14, error)

call h5screate_simple_f(RANK2, dims_cloudphase, dsp_id_Ygrad_Emiss14, error)
call h5dcreate_f(file_id, dset_Ygrad_Emiss14_name, H5T_NATIVE_INTEGER, dsp_id_Ygrad_Emiss14, sds_id_Ygrad_Emiss14, error)
call h5dwrite_f (sds_id_Ygrad_Emiss14, H5T_NATIVE_INTEGER, clp%Ygrad_Emiss14, dims_cloudphase, error)
call h5dclose_f (sds_id_Ygrad_Emiss14, error)      
call h5sclose_f (dsp_id_Ygrad_Emiss14, error)

call h5screate_simple_f(RANK2, dims_cloudphase, dsp_id_Num_Steps_Gradient, error)
call h5dcreate_f(file_id, dset_Num_Steps_Gradient_name, H5T_NATIVE_INTEGER, dsp_id_Num_Steps_Gradient, sds_id_Num_Steps_Gradient, error)
call h5dwrite_f (sds_id_Num_Steps_Gradient, H5T_NATIVE_INTEGER, clp%Num_Steps_Gradient, dims_cloudphase, error)
call h5dclose_f (sds_id_Num_Steps_Gradient, error)      
call h5sclose_f (dsp_id_Num_Steps_Gradient, error)

call h5screate_simple_f(RANK2, dims_cloudphase, dsp_id_Emiss_Chn10_Tot, error)
call h5dcreate_f(file_id, dset_Emiss_Chn10_Tot_name, H5T_NATIVE_REAL, dsp_id_Emiss_Chn10_Tot, sds_id_Emiss_Chn10_Tot, error)
call h5dwrite_f (sds_id_Emiss_Chn10_Tot, H5T_NATIVE_REAL, clp%Emiss_Chn10_Tot, dims_cloudphase, error)
call h5dclose_f (sds_id_Emiss_Chn10_Tot, error)      
call h5sclose_f (dsp_id_Emiss_Chn10_Tot, error)

call h5screate_simple_f(RANK2, dims_cloudphase, dsp_id_Emiss_Chn10_Tot_Lrc, error)
call h5dcreate_f(file_id, dset_Emiss_Chn10_Tot_Lrc_name, H5T_NATIVE_REAL, dsp_id_Emiss_Chn10_Tot_Lrc, sds_id_Emiss_Chn10_Tot_Lrc, error)
call h5dwrite_f (sds_id_Emiss_Chn10_Tot_Lrc, H5T_NATIVE_REAL, clp%Emiss_Chn10_Tot_Lrc, dims_cloudphase, error)
call h5dclose_f (sds_id_Emiss_Chn10_Tot_Lrc, error)      
call h5sclose_f (dsp_id_Emiss_Chn10_Tot_Lrc, error)

call h5screate_simple_f(RANK2, dims_cloudphase, dsp_id_Emiss_Chn14_Tot, error)
call h5dcreate_f(file_id, dset_Emiss_Chn14_Tot_name, H5T_NATIVE_REAL, dsp_id_Emiss_Chn14_Tot, sds_id_Emiss_Chn14_Tot, error)
call h5dwrite_f (sds_id_Emiss_Chn14_Tot, H5T_NATIVE_REAL, clp%Emiss_Chn14_Tot, dims_cloudphase, error)
call h5dclose_f (sds_id_Emiss_Chn14_Tot, error)      
call h5sclose_f (dsp_id_Emiss_Chn14_Tot, error)

call h5screate_simple_f(RANK2, dims_cloudphase, dsp_id_Emiss_Chn14_Tot_Multi, error)
call h5dcreate_f(file_id, dset_Emiss_Chn14_Tot_Multi_name, H5T_NATIVE_REAL, dsp_id_Emiss_Chn14_Tot_Multi, sds_id_Emiss_Chn14_Tot_Multi, error)
call h5dwrite_f (sds_id_Emiss_Chn14_Tot_Multi, H5T_NATIVE_REAL, clp%Emiss_Chn14_Tot_Multi, dims_cloudphase, error)
call h5dclose_f (sds_id_Emiss_Chn14_Tot_Multi, error)      
call h5sclose_f (dsp_id_Emiss_Chn14_Tot_Multi, error)

call h5screate_simple_f(RANK2, dims_cloudphase, dsp_id_Cldbeta7311_Tot, error)
call h5dcreate_f(file_id, dset_Cldbeta7311_Tot_name, H5T_NATIVE_REAL, dsp_id_Cldbeta7311_Tot, sds_id_Cldbeta7311_Tot, error)
call h5dwrite_f (sds_id_Cldbeta7311_Tot, H5T_NATIVE_REAL, clp%Cldbeta7311_Tot, dims_cloudphase, error)
call h5dclose_f (sds_id_Cldbeta7311_Tot, error)      
call h5sclose_f (dsp_id_Cldbeta7311_Tot, error)

call h5screate_simple_f(RANK2, dims_cloudphase, dsp_id_Cldbeta7311_Tot_Multi, error)
call h5dcreate_f(file_id, dset_Cldbeta7311_Tot_Multi_name, H5T_NATIVE_REAL, dsp_id_Cldbeta7311_Tot_Multi, sds_id_Cldbeta7311_Tot_Multi, error)
call h5dwrite_f (sds_id_Cldbeta7311_Tot_Multi, H5T_NATIVE_REAL, clp%Cldbeta7311_Tot_Multi, dims_cloudphase, error)
call h5dclose_f (sds_id_Cldbeta7311_Tot_Multi, error)      
call h5sclose_f (dsp_id_Cldbeta7311_Tot_Multi, error)

call h5screate_simple_f(RANK2, dims_cloudphase, dsp_id_Cldbeta1112_Tot, error)
call h5dcreate_f(file_id, dset_Cldbeta1112_Tot_name, H5T_NATIVE_REAL, dsp_id_Cldbeta1112_Tot, sds_id_Cldbeta1112_Tot, error)
call h5dwrite_f (sds_id_Cldbeta1112_Tot, H5T_NATIVE_REAL, clp%Cldbeta1112_Tot, dims_cloudphase, error)
call h5dclose_f (sds_id_Cldbeta1112_Tot, error)      
call h5sclose_f (dsp_id_Cldbeta1112_Tot, error)

call h5screate_simple_f(RANK2, dims_cloudphase, dsp_id_Cldbeta1112_Tot_Lrc, error)
call h5dcreate_f(file_id, dset_Cldbeta1112_Tot_Lrc_name, H5T_NATIVE_REAL, dsp_id_Cldbeta1112_Tot_Lrc, sds_id_Cldbeta1112_Tot_Lrc, error)
call h5dwrite_f (sds_id_Cldbeta1112_Tot_Lrc, H5T_NATIVE_REAL, clp%Cldbeta1112_Tot_Lrc, dims_cloudphase, error)
call h5dclose_f (sds_id_Cldbeta1112_Tot_Lrc, error)      
call h5sclose_f (dsp_id_Cldbeta1112_Tot_Lrc, error)

call h5screate_simple_f(RANK2, dims_cloudphase, dsp_id_Cldbeta1112_Tot_Multi, error)
call h5dcreate_f(file_id, dset_Cldbeta1112_Tot_Multi_name, H5T_NATIVE_REAL, dsp_id_Cldbeta1112_Tot_Multi, sds_id_Cldbeta1112_Tot_Multi, error)
call h5dwrite_f (sds_id_Cldbeta1112_Tot_Multi, H5T_NATIVE_REAL, clp%Cldbeta1112_Tot_Multi, dims_cloudphase, error)
call h5dclose_f (sds_id_Cldbeta1112_Tot_Multi, error)      
call h5sclose_f (dsp_id_Cldbeta1112_Tot_Multi, error)

call h5screate_simple_f(RANK2, dims_cloudphase, dsp_id_Cldbeta1112_Opaque, error)
call h5dcreate_f(file_id, dset_Cldbeta1112_Opaque_name, H5T_NATIVE_REAL, dsp_id_Cldbeta1112_Opaque, sds_id_Cldbeta1112_Opaque, error)
call h5dwrite_f (sds_id_Cldbeta1112_Opaque, H5T_NATIVE_REAL, clp%Cldbeta1112_Opaque, dims_cloudphase, error)
call h5dclose_f (sds_id_Cldbeta1112_Opaque, error)      
call h5sclose_f (dsp_id_Cldbeta1112_Opaque, error)

call h5screate_simple_f(RANK2, dims_cloudphase, dsp_id_Cldbeta1112_Opaque_Multi, error)
call h5dcreate_f(file_id, dset_Cldbeta1112_Opaque_Multi_name, H5T_NATIVE_REAL, dsp_id_Cldbeta1112_Opaque_Multi, sds_id_Cldbeta1112_Opaque_Multi, error)
call h5dwrite_f (sds_id_Cldbeta1112_Opaque_Multi, H5T_NATIVE_REAL, clp%Cldbeta1112_Opaque_Multi, dims_cloudphase, error)
call h5dclose_f (sds_id_Cldbeta1112_Opaque_Multi, error)      
call h5sclose_f (dsp_id_Cldbeta1112_Opaque_Multi, error)

call h5screate_simple_f(RANK2, dims_cloudphase, dsp_id_Cldbeta8511_Tot, error)
call h5dcreate_f(file_id, dset_Cldbeta8511_Tot_name, H5T_NATIVE_REAL, dsp_id_Cldbeta8511_Tot, sds_id_Cldbeta8511_Tot, error)
call h5dwrite_f (sds_id_Cldbeta8511_Tot, H5T_NATIVE_REAL, clp%Cldbeta8511_Tot, dims_cloudphase, error)
call h5dclose_f (sds_id_Cldbeta8511_Tot, error)      
call h5sclose_f (dsp_id_Cldbeta8511_Tot, error)

call h5screate_simple_f(RANK2, dims_cloudphase, dsp_id_Cldbeta8511_Tot_Multi, error)
call h5dcreate_f(file_id, dset_Cldbeta8511_Tot_Multi_name, H5T_NATIVE_REAL, dsp_id_Cldbeta8511_Tot_Multi, sds_id_Cldbeta8511_Tot_Multi, error)
call h5dwrite_f (sds_id_Cldbeta8511_Tot_Multi, H5T_NATIVE_REAL, clp%Cldbeta8511_Tot_Multi, dims_cloudphase, error)
call h5dclose_f (sds_id_Cldbeta8511_Tot_Multi, error)      
call h5sclose_f (dsp_id_Cldbeta8511_Tot_Multi, error)

call h5screate_simple_f(RANK2, dims_cloudphase, dsp_id_Cldbeta8511_Opaque, error)
call h5dcreate_f(file_id, dset_Cldbeta8511_Opaque_name, H5T_NATIVE_REAL, dsp_id_Cldbeta8511_Opaque, sds_id_Cldbeta8511_Opaque, error)
call h5dwrite_f (sds_id_Cldbeta8511_Opaque, H5T_NATIVE_REAL, clp%Cldbeta8511_Opaque, dims_cloudphase, error)
call h5dclose_f (sds_id_Cldbeta8511_Opaque, error)      
call h5sclose_f (dsp_id_Cldbeta8511_Opaque, error)

call h5screate_simple_f(RANK2, dims_cloudphase, dsp_id_Cldbeta8511_Opaque_Lrc, error)
call h5dcreate_f(file_id, dset_Cldbeta8511_Opaque_Lrc_name, H5T_NATIVE_REAL, dsp_id_Cldbeta8511_Opaque_Lrc, sds_id_Cldbeta8511_Opaque_Lrc, error)
call h5dwrite_f (sds_id_Cldbeta8511_Opaque_Lrc, H5T_NATIVE_REAL, clp%Cldbeta8511_Opaque_Lrc, dims_cloudphase, error)
call h5dclose_f (sds_id_Cldbeta8511_Opaque_Lrc, error)      
call h5sclose_f (dsp_id_Cldbeta8511_Opaque_Lrc, error)

call h5screate_simple_f(RANK2, dims_cloudphase, dsp_id_Cldbeta8511_Opaque_Multi, error)
call h5dcreate_f(file_id, dset_Cldbeta8511_Opaque_Multi_name, H5T_NATIVE_REAL, dsp_id_Cldbeta8511_Opaque_Multi, sds_id_Cldbeta8511_Opaque_Multi, error)
call h5dwrite_f (sds_id_Cldbeta8511_Opaque_Multi, H5T_NATIVE_REAL, clp%Cldbeta8511_Opaque_Multi, dims_cloudphase, error)
call h5dclose_f (sds_id_Cldbeta8511_Opaque_Multi, error)      
call h5sclose_f (dsp_id_Cldbeta8511_Opaque_Multi, error)

call h5screate_simple_f(RANK2, dims_cloudphase, dsp_id_Opaque_Cld_Temp_Chn10, error)
call h5dcreate_f(file_id, dset_Opaque_Cld_Temp_Chn10_name, H5T_NATIVE_REAL, dsp_id_Opaque_Cld_Temp_Chn10, sds_id_Opaque_Cld_Temp_Chn10, error)
call h5dwrite_f (sds_id_Opaque_Cld_Temp_Chn10, H5T_NATIVE_REAL, clp%Opaque_Cld_Temp_Chn10, dims_cloudphase, error)
call h5dclose_f (sds_id_Opaque_Cld_Temp_Chn10, error)      
call h5sclose_f (dsp_id_Opaque_Cld_Temp_Chn10, error)

call h5screate_simple_f(RANK2, dims_cloudphase, dsp_id_Opaque_Cld_Temp_Chn14, error)
call h5dcreate_f(file_id, dset_Opaque_Cld_Temp_Chn14_name, H5T_NATIVE_REAL, dsp_id_Opaque_Cld_Temp_Chn14, sds_id_Opaque_Cld_Temp_Chn14, error)
call h5dwrite_f (sds_id_Opaque_Cld_Temp_Chn14, H5T_NATIVE_REAL, clp%Opaque_Cld_Temp_Chn14, dims_cloudphase, error)
call h5dclose_f (sds_id_Opaque_Cld_Temp_Chn14, error)      
call h5sclose_f (dsp_id_Opaque_Cld_Temp_Chn14, error)


endif 
! close file
call h5fclose_f(file_id, error)

!----------------------------------------------------

call h5close_f(error)

deallocate(output_cloudp,   &
           Cldphase_Qf,     &
           Cldphase_Qpi)

end subroutine fylat_write_out_cloud_phase
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~ Subroutine  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine fylat_write_out_cloud_height()

!-----------------------------------------------------------------------
! !F90 write_out_arrays
!
! !Description:
!    The main program for write cloud mask array out.
!
! !Input  parameters:
!    filename        = 
!
! !Output parameters:
!    none
!-----------------------------------------------------------------------


IMPLICIT NONE

! 1. define variables
! cloud mask
integer :: error
integer (HID_T) :: file_id                 ! File identifier for cloud mask
  
integer (HID_T) :: sds_id_cldp        ! sds id for cloud mask
integer (HID_T) :: sds_id_cldt        ! sds id for cloud mask qa
integer (HID_T) :: sds_id_cldz        ! sds id for cloud mask
integer (HID_T) :: sds_id_cldemiss

!--others
integer (HID_T) :: sds_id_cod_vis      ! sds id for cloud mask qa
integer (HID_T) :: sds_id_cldbeta1112  ! sds id for cloud mask qa

!-----
integer (HID_T) :: dsp_id_cldp         ! dsp id for cloud optical depth
integer (HID_T) :: dsp_id_cldt         ! dsp id for cloud particle effective radius [unit:um]
integer (HID_T) :: dsp_id_cldz         ! dsp id for cloud optical depth
integer (HID_T) :: dsp_id_cldemiss     ! dsp id for cloud particle effective radius [unit:um]

integer (HID_T) :: dsp_id_cod_vis          ! dsp id for cloud particle effective radius [unit:um]
integer (HID_T) :: dsp_id_cldbeta1112         ! dsp id for cloud particle effective radius [unit:um]


!-----
character (LEN=18), parameter           :: dset_cldp_name       = "Cloud_Top_Pressure"        
character (LEN=21), parameter           :: dset_cldt_name       = "Cloud_Top_Temperature" 
character (LEN=16), parameter           :: dset_cldz_name       = "Cloud_Top_Height"        
character (LEN=24), parameter           :: dset_cldemiss_name   = "Cloud_Emissivity_at_11um" 
! others
character (LEN=13), parameter           :: dset_cod_vis_name        = "cod_vis_11_12" 
character (LEN=11), parameter           :: dset_cldbeta1112_name    = "cldbeta1112" 


integer (HSIZE_T), dimension(2)         :: dims_ctp


!integer(kind=4) :: error     ! Error flag for hdf5 open
integer(kind=4), parameter :: RANK2 = 2 ! Dataset rank
integer(kind=4), parameter :: RANK3 = 3 ! Dataset rank
integer(kind=4),dimension(:,:),allocatable   :: output_int
character,dimension(:,:),allocatable   :: output_cloudp
character,dimension(:,:,:),allocatable :: Cldphase_Qf
character,dimension(:,:,:),allocatable :: Cldphase_Qpi
character(len=20) :: satellite_name

!---------------------------------------------------
! attribute
integer(HSIZE_T),dimension(1) :: a_dims
integer(kind=4)               :: a_rank
integer(HID_T)                :: output_a_id,output_aspace_id,output_atype_id
! 2. begin program

!+++++++++++++++++++++++++++ Step 2: Write product start +++++++++++++++++++++++++++
!== initialized

!--- variables
!allocate(output_cloudp(sat%nElem,sat%nLine),    &
!         Cldphase_Qf(6,sat%nElem,sat%nLine),    &
!         Cldphase_Qpi(20,sat%nElem,sat%nLine) )
allocate(output_int(sat%nElem,sat%nLine))

!------         
dims_ctp(1) = sat%nElem
dims_ctp(2) = sat%nLine


!== 2.1. write start
print*,'    ... fylat write out fy3/MERSI_II Cloud Height HDF5 product !!! '
	 
call h5open_f(error)

call h5fcreate_f(trim(fy3_mersi_CTP_data), H5F_ACC_TRUNC_F, file_id, error)

!------------------------------------
! --- Write Cloud Phase and Type
!call h5dcreate_f(group_id, "EV_1KM_RefSB", H5T_STD_U16LE, dsp_id, sds_id, error)
!call h5dwrite_f (sds_id, H5T_NATIVE_INTEGER, rad(:,:,5:19), dims_sp32, error)
!output_cloudp = char(clp%Cldphase)!cm_bitarray, cm_qa_bitarray
output_int = int(ctp%cldp*100.0)
call h5screate_simple_f(RANK2, dims_ctp, dsp_id_cldp, error)
call h5dcreate_f(file_id, dset_cldp_name, H5T_NATIVE_INTEGER, dsp_id_cldp, sds_id_cldp, error)
!call h5dwrite_f (sds_id_cldp, H5T_NATIVE_REAL, ctp%cldp, dims_ctp, error)
call h5dwrite_f (sds_id_cldp, H5T_NATIVE_INTEGER, output_int, dims_ctp, error)

! ... write the attributes
      ! .. units attribute
      a_rank = 1
      a_dims = 1
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cldp,'units',H5T_NATIVE_CHARACTER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_CHARACTER, 'hPa', a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. valid_range attribute
      a_dims = 2
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cldp,'valid_range',H5T_NATIVE_INTEGER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_INTEGER, (/0,2000/), a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. _FillValue attibute
      a_dims = 1 
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cldp,'_FillValue',H5T_NATIVE_INTEGER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_INTEGER, -999.9 , a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. Intercept attribute
      a_dims = 1 
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cldp,'Intercept',H5T_NATIVE_REAL, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_REAL, 0.0 , a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)     
      ! .. Slope attribute
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cldp,'Slope',H5T_NATIVE_REAL, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_REAL, 0.01 , a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. long_name
      a_dims = 1
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cldp,'long_name',H5T_NATIVE_CHARACTER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_CHARACTER, '5-min granule Cloud Top Pressure', a_dims, error)
      !call h5awrite_f(output_a_id, H5T_NATIVE_SHORT, 'fylat MERSI_II Cloud Mask', a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. band_name
      a_dims = 1 
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cldp,'band_name',H5T_NATIVE_CHARACTER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_CHARACTER, 'NANA', a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      
      ! Close output dataset ,dataspace, file

call h5dclose_f (sds_id_cldp, error)      
call h5sclose_f (dsp_id_cldp, error)

! cloud type
output_int = int(ctp%cldt*100.0)
call h5screate_simple_f(RANK2, dims_ctp, dsp_id_cldt, error)
call h5dcreate_f(file_id, dset_cldt_name, H5T_NATIVE_INTEGER, dsp_id_cldt, sds_id_cldt, error)
!call h5dwrite_f (sds_id_cldt, H5T_NATIVE_REAL, ctp%cldt, dims_ctp, error)
call h5dwrite_f (sds_id_cldt, H5T_NATIVE_INTEGER, output_int, dims_ctp, error)

! ... write the attributes
      ! .. units attribute
      a_rank = 1
      a_dims = 1
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cldt,'units',H5T_NATIVE_CHARACTER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_CHARACTER, 'K', a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. valid_range attribute
      a_dims = 2
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cldt,'valid_range',H5T_NATIVE_INTEGER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_INTEGER, (/0,40000/), a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. _FillValue attibute
      a_dims = 1 
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cldt,'_FillValue',H5T_NATIVE_INTEGER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_INTEGER, -999.9 , a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. Intercept attribute
      a_dims = 1 
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cldt,'Intercept',H5T_NATIVE_REAL, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_REAL, 0.0 , a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)     
      ! .. Slope attribute
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cldt,'Slope',H5T_NATIVE_REAL, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_REAL, 0.01 , a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. long_name
      a_dims = 1
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cldt,'long_name',H5T_NATIVE_CHARACTER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_CHARACTER, '5-min granule Cloud Top Temperature', a_dims, error)
      !call h5awrite_f(output_a_id, H5T_NATIVE_SHORT, 'fylat MERSI_II Cloud Mask', a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. band_name
      a_dims = 1 
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cldt,'band_name',H5T_NATIVE_CHARACTER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_CHARACTER, 'NANA', a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      
      ! Close output dataset ,dataspace, file

call h5dclose_f (sds_id_cldt, error)      
call h5sclose_f (dsp_id_cldt, error)


output_int = int(ctp%cldz*10.0)
call h5screate_simple_f(RANK2, dims_ctp, dsp_id_cldz, error)
call h5dcreate_f(file_id, dset_cldz_name, H5T_NATIVE_INTEGER, dsp_id_cldz, sds_id_cldz, error)
call h5dwrite_f (sds_id_cldz, H5T_NATIVE_INTEGER, output_int, dims_ctp, error)
! ... write the attributes
      ! .. units attribute
      a_rank = 1
      a_dims = 1
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cldz,'units',H5T_NATIVE_CHARACTER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_CHARACTER, 'm', a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. valid_range attribute
      a_dims = 2
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cldz,'valid_range',H5T_NATIVE_INTEGER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_INTEGER, (/0,30000/), a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. _FillValue attibute
      a_dims = 1 
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cldz,'_FillValue',H5T_NATIVE_INTEGER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_INTEGER, -999.9 , a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. Intercept attribute
      a_dims = 1 
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cldz,'Intercept',H5T_NATIVE_REAL, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_REAL, 0.0 , a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)     
      ! .. Slope attribute
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cldz,'Slope',H5T_NATIVE_REAL, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_REAL, 0.1 , a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. long_name
      a_dims = 1
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cldz,'long_name',H5T_NATIVE_CHARACTER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_CHARACTER, '5-min granule Cloud Top Height', a_dims, error)
      !call h5awrite_f(output_a_id, H5T_NATIVE_SHORT, 'fylat MERSI_II Cloud Mask', a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. band_name
      a_dims = 1 
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cldz,'band_name',H5T_NATIVE_CHARACTER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_CHARACTER, 'NANA', a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      
      ! Close output dataset ,dataspace, file

call h5dclose_f (sds_id_cldz, error)      
call h5sclose_f (dsp_id_cldz, error)


output_int = int(ctp%cldemiss*10000.0)
call h5screate_simple_f(RANK2, dims_ctp, dsp_id_cldemiss, error)
call h5dcreate_f(file_id, dset_cldemiss_name, H5T_NATIVE_INTEGER, dsp_id_cldemiss, sds_id_cldemiss, error)
call h5dwrite_f (sds_id_cldemiss, H5T_NATIVE_INTEGER, output_int, dims_ctp, error)
! ... write the attributes
      ! .. units attribute
      a_rank = 1
      a_dims = 1
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cldemiss,'units',H5T_NATIVE_CHARACTER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_CHARACTER, 'm', a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. valid_range attribute
      a_dims = 2
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cldemiss,'valid_range',H5T_NATIVE_INTEGER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_INTEGER, (/0,10000/), a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. _FillValue attibute
      a_dims = 1 
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cldemiss,'_FillValue',H5T_NATIVE_INTEGER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_INTEGER, -999.9 , a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. Intercept attribute
      a_dims = 1 
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cldemiss,'Intercept',H5T_NATIVE_REAL, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_REAL, 1.0 , a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)     
      ! .. Slope attribute
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cldemiss,'Slope',H5T_NATIVE_REAL, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_REAL, 0.0001 , a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. long_name
      a_dims = 1
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cldemiss,'long_name',H5T_NATIVE_CHARACTER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_CHARACTER, '5-min granule Cloud Emissivity', a_dims, error)
      !call h5awrite_f(output_a_id, H5T_NATIVE_SHORT, 'fylat MERSI_II Cloud Mask', a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. band_name
      a_dims = 1 
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cldemiss,'band_name',H5T_NATIVE_CHARACTER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_CHARACTER, 'NANA', a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      
      ! Close output dataset ,dataspace, file

call h5dclose_f (sds_id_cldemiss, error)      
call h5sclose_f (dsp_id_cldemiss, error)


!------------------------------------
! --- Write others out

if (fylat_alg_opt%cloudtopz_index > 1) then

call h5screate_simple_f(RANK2, dims_ctp, dsp_id_cod_vis, error)
call h5dcreate_f(file_id, dset_cod_vis_name, H5T_NATIVE_REAL, dsp_id_cod_vis, sds_id_cod_vis, error)
call h5dwrite_f (sds_id_cod_vis, H5T_NATIVE_REAL, ctp%cod_vis, dims_ctp, error)
call h5dclose_f (sds_id_cod_vis, error)      
call h5sclose_f (dsp_id_cod_vis, error)

call h5screate_simple_f(RANK2, dims_ctp, dsp_id_cldbeta1112, error)
call h5dcreate_f(file_id, dset_cldbeta1112_name, H5T_NATIVE_REAL, dsp_id_cldbeta1112, sds_id_cldbeta1112, error)
call h5dwrite_f (sds_id_cldbeta1112, H5T_NATIVE_REAL, ctp%cldbeta1112, dims_ctp, error)
call h5dclose_f (sds_id_cldbeta1112, error)      
call h5sclose_f (dsp_id_cldbeta1112, error)

endif 
! close file
call h5fclose_f(file_id, error)

!----------------------------------------------------

call h5close_f(error)

deallocate(output_int)


end subroutine fylat_write_out_cloud_height
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~ Subroutine  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine fylat_write_out_cloud_micro_day()

!-----------------------------------------------------------------------
! !F90 write_out_arrays
!
! !Description:
!    The main program for write cloud mask array out.
!
! !Input  parameters:
!    filename        = 
!
! !Output parameters:
!    none
!-----------------------------------------------------------------------


IMPLICIT NONE

! 1. define variables
! cloud mask
integer :: error
integer (HID_T) :: file_id                 ! File identifier for cloud mask

integer (HID_T) :: sds_id_cod_vis        ! sds id for cloud mask
integer (HID_T) :: sds_id_cldreff        ! sds id for cloud mask qa
integer (HID_T) :: sds_id_cldlwp         ! sds id for cloud mask
integer (HID_T) :: sds_id_cldiwp
integer (HID_T) :: sds_id_qcflg_cotd

!-----
integer (HID_T) :: dsp_id_cod_vis         ! dsp id for cloud optical depth
integer (HID_T) :: dsp_id_cldreff         ! dsp id for cloud particle effective radius [unit:um]
integer (HID_T) :: dsp_id_cldlwp          ! dsp id for cloud optical depth
integer (HID_T) :: dsp_id_cldiwp          ! dsp id for cloud particle effective radius [unit:um]
integer (HID_T) :: dsp_id_qcflg_cotd          ! dsp id for cloud particle effective radius [unit:um]


!-----
character (LEN=19), parameter           :: dset_cod_vis_name       = "Cloud_Optical_Depth"        
character (LEN=22), parameter           :: dset_cldreff_name       = "Cloud_Effective_Radius" 
character (LEN=23), parameter           :: dset_cldlwp_name        = "Cloud_Liquid_Water_Path"        
character (LEN=20), parameter           :: dset_cldiwp_name        = "Cloud_Ice_Water_Path" 
character (LEN=12), parameter           :: dset_qcflg_cotd_name    = "Cldot_qcflag"  


integer (HSIZE_T), dimension(2)         :: dims_cot
integer (HSIZE_T), dimension(3)         :: dims_qf

!integer(kind=4) :: error     ! Error flag for hdf5 open
integer(kind=4),dimension(:,:),allocatable   :: output_int
integer(kind=4), parameter :: RANK2 = 2 ! Dataset rank
integer(kind=4), parameter :: RANK3 = 3 ! Dataset rank
character,dimension(:,:),allocatable   :: output_cloudp
character,dimension(:,:,:),allocatable :: Cld_Qf
character,dimension(:,:,:),allocatable :: Cldphase_Qpi
character(len=20) :: satellite_name

!---------------------------------------------------
! attribute
integer(HSIZE_T),dimension(1) :: a_dims
integer(kind=4)               :: a_rank
integer(HID_T)                :: output_a_id,output_aspace_id,output_atype_id
! 2. begin program

!+++++++++++++++++++++++++++ Step 2: Write product start +++++++++++++++++++++++++++
!== initialized

!--- variables
allocate( Cld_Qf(2,sat%nElem,sat%nLine), &
          output_int(sat%nElem,sat%nLine))

!------         
dims_cot(1) = sat%nElem
dims_cot(2) = sat%nLine
dims_qf(1)  = 2
dims_qf(2)  = sat%nElem
dims_qf(3)  = sat%nLine

!== 2.1. write start
print*,'    ... fylat write out fy3/MERSI_II Cloud Microphysical and Optical Properties at daytime HDF5 product !!! '
	 
call h5open_f(error)

call h5fcreate_f(trim(fy3_mersi_COT_data), H5F_ACC_TRUNC_F, file_id, error)

!------------------------------------
! --- Write Cloud Phase and Type
!call h5dcreate_f(group_id, "EV_1KM_RefSB", H5T_STD_U16LE, dsp_id, sds_id, error)
!call h5dwrite_f (sds_id, H5T_NATIVE_INTEGER, rad(:,:,5:19), dims_sp32, error)
!output_cloudp = char(clp%Cldphase)!cm_bitarray, cm_qa_bitarray
output_int = int(cot%cod_vis*100.0)
call h5screate_simple_f(RANK2, dims_cot, dsp_id_cod_vis, error)
call h5dcreate_f(file_id, dset_cod_vis_name, H5T_NATIVE_INTEGER, dsp_id_cod_vis, sds_id_cod_vis, error)
call h5dwrite_f (sds_id_cod_vis, H5T_NATIVE_INTEGER, output_int, dims_cot, error)

! ... write the attributes
      ! .. units attribute
      a_rank = 1
      a_dims = 1
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cod_vis,'units',H5T_NATIVE_CHARACTER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_CHARACTER, ' ', a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. valid_range attribute
      a_dims = 2
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cod_vis,'valid_range',H5T_NATIVE_INTEGER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_INTEGER, (/0,15000/), a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. _FillValue attibute
      a_dims = 1 
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cod_vis,'_FillValue',H5T_NATIVE_INTEGER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_INTEGER, -999.9 , a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. Intercept attribute
      a_dims = 1 
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cod_vis,'Intercept',H5T_NATIVE_REAL, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_REAL, 0.0 , a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)     
      ! .. Slope attribute
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cod_vis,'Slope',H5T_NATIVE_REAL, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_REAL, 0.01 , a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. long_name
      a_dims = 1
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cod_vis,'long_name',H5T_NATIVE_CHARACTER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_CHARACTER, '5-min granule Cloud Optical Depth', a_dims, error)
      !call h5awrite_f(output_a_id, H5T_NATIVE_SHORT, 'fylat MERSI_II Cloud Mask', a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. band_name
      a_dims = 1 
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cod_vis,'band_name',H5T_NATIVE_CHARACTER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_CHARACTER, 'NANA', a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      
      ! Close output dataset ,dataspace, file

call h5dclose_f (sds_id_cod_vis, error)      
call h5sclose_f (dsp_id_cod_vis, error)

! cloud reff
output_int = int(cot%cldreff*100.0)
call h5screate_simple_f(RANK2, dims_cot, dsp_id_cldreff, error)
call h5dcreate_f(file_id, dset_cldreff_name, H5T_NATIVE_INTEGER, dsp_id_cldreff, sds_id_cldreff, error)
call h5dwrite_f (sds_id_cldreff, H5T_NATIVE_INTEGER, output_int, dims_cot, error)

! ... write the attributes
      ! .. units attribute
      a_rank = 1
      a_dims = 1
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cldreff,'units',H5T_NATIVE_CHARACTER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_CHARACTER, 'um', a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. valid_range attribute
      a_dims = 2
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cldreff,'valid_range',H5T_NATIVE_INTEGER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_INTEGER, (/0,10000/), a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. _FillValue attibute
      a_dims = 1 
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cldreff,'_FillValue',H5T_NATIVE_INTEGER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_INTEGER, -999.9 , a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. Intercept attribute
      a_dims = 1 
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cldreff,'Intercept',H5T_NATIVE_REAL, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_REAL, 0.0 , a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)     
      ! .. Slope attribute
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cldreff,'Slope',H5T_NATIVE_REAL, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_REAL, 0.01 , a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. long_name
      a_dims = 1
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cldreff,'long_name',H5T_NATIVE_CHARACTER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_CHARACTER, '5-min granule Cloud Top Temperature', a_dims, error)
      !call h5awrite_f(output_a_id, H5T_NATIVE_SHORT, 'fylat MERSI_II Cloud Mask', a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. band_name
      a_dims = 1 
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cldreff,'band_name',H5T_NATIVE_CHARACTER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_CHARACTER, 'NANA', a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      
      ! Close output dataset ,dataspace, file

call h5dclose_f (sds_id_cldreff, error)      
call h5sclose_f (dsp_id_cldreff, error)


output_int = int(cot%cldlwp*10.0)
call h5screate_simple_f(RANK2, dims_cot, dsp_id_cldlwp, error)
call h5dcreate_f(file_id, dset_cldlwp_name, H5T_NATIVE_INTEGER, dsp_id_cldlwp, sds_id_cldlwp, error)
call h5dwrite_f (sds_id_cldlwp, H5T_NATIVE_INTEGER, output_int, dims_cot, error)
! ... write the attributes
      ! .. units attribute
      a_rank = 1
      a_dims = 1
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cldlwp,'units',H5T_NATIVE_CHARACTER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_CHARACTER, 'g/m2', a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. valid_range attribute
      a_dims = 2
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cldlwp,'valid_range',H5T_NATIVE_INTEGER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_INTEGER, (/0,10000/), a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. _FillValue attibute
      a_dims = 1 
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cldlwp,'_FillValue',H5T_NATIVE_INTEGER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_INTEGER, -999.9 , a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. Intercept attribute
      a_dims = 1 
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cldlwp,'Intercept',H5T_NATIVE_REAL, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_REAL, 0.0 , a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)     
      ! .. Slope attribute
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cldlwp,'Slope',H5T_NATIVE_REAL, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_REAL, 0.1 , a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. long_name
      a_dims = 1
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cldlwp,'long_name',H5T_NATIVE_CHARACTER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_CHARACTER, '5-min granule Cloud Liquid Water Path', a_dims, error)
      !call h5awrite_f(output_a_id, H5T_NATIVE_SHORT, 'fylat MERSI_II Cloud Mask', a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. band_name
      a_dims = 1 
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cldlwp,'band_name',H5T_NATIVE_CHARACTER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_CHARACTER, 'NANA', a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      
      ! Close output dataset ,dataspace, file

call h5dclose_f (sds_id_cldlwp, error)      
call h5sclose_f (dsp_id_cldlwp, error)


output_int = int(cot%cldiwp*10.0)
call h5screate_simple_f(RANK2, dims_cot, dsp_id_cldiwp, error)
call h5dcreate_f(file_id, dset_cldiwp_name, H5T_NATIVE_INTEGER, dsp_id_cldiwp, sds_id_cldiwp, error)
call h5dwrite_f (sds_id_cldiwp, H5T_NATIVE_INTEGER, output_int, dims_cot, error)
! ... write the attributes
      ! .. units attribute
      a_rank = 1
      a_dims = 1
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cldiwp,'units',H5T_NATIVE_CHARACTER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_CHARACTER, 'g/m2', a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. valid_range attribute
      a_dims = 2
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cldiwp,'valid_range',H5T_NATIVE_INTEGER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_INTEGER, (/0,10000/), a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. _FillValue attibute
      a_dims = 1 
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cldiwp,'_FillValue',H5T_NATIVE_INTEGER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_INTEGER, -999.9 , a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. Intercept attribute
      a_dims = 1 
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cldiwp,'Intercept',H5T_NATIVE_REAL, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_REAL, 0.0 , a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)     
      ! .. Slope attribute
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cldiwp,'Slope',H5T_NATIVE_REAL, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_REAL, 0.1 , a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. long_name
      a_dims = 1
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cldiwp,'long_name',H5T_NATIVE_CHARACTER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_CHARACTER, '5-min granule Cloud Ice Water Path', a_dims, error)
      !call h5awrite_f(output_a_id, H5T_NATIVE_SHORT, 'fylat MERSI_II Cloud Mask', a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      ! .. band_name
      a_dims = 1 
      call h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      call h5acreate_f(sds_id_cldiwp,'band_name',H5T_NATIVE_CHARACTER, output_aspace_id, output_a_id, error)
      call h5awrite_f(output_a_id, H5T_NATIVE_CHARACTER, 'NANA', a_dims, error)
      call h5sclose_f(output_aspace_id, error)
      call h5aclose_f(output_a_id, error)
      
      ! Close output dataset ,dataspace, file

call h5dclose_f (sds_id_cldiwp, error)      
call h5sclose_f (dsp_id_cldiwp, error)

Cld_Qf = char(cot%qcflg_cotd)
call h5screate_simple_f(RANK3, dims_qf, dsp_id_qcflg_cotd, error)
call h5dcreate_f(file_id, dset_qcflg_cotd_name, H5T_STD_U8LE, dsp_id_qcflg_cotd, sds_id_qcflg_cotd, error)
call h5dwrite_f (sds_id_qcflg_cotd, H5T_STD_U8LE, Cld_Qf, dims_qf, error)
call h5dclose_f (sds_id_qcflg_cotd, error)      
call h5sclose_f (dsp_id_qcflg_cotd, error)


deallocate(Cld_Qf, output_int)
           
! close file
call h5fclose_f(file_id, error)

!----------------------------------------------------

call h5close_f(error)


end subroutine fylat_write_out_cloud_micro_day
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!~~~~~~~~~~~~~~~~ Subroutine  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine allocate_cloudphase_arrays()

integer :: error

allocate (clp%Cldtype(sat%nElem, sat%nLine),                  &
          clp%Cldtype_Tmpy(sat%nElem, sat%nLine),             &        
          clp%Cldphase(sat%nElem, sat%nLine),                 & 
          clp%Cldphase_Qf(6, sat%nElem, sat%nLine),           & 
          clp%Cldphase_Qpi(20,sat%nElem, sat%nLine),          & 
          clp%Xgrad_Emiss14(sat%nElem, sat%nLine),            &   
          clp%Ygrad_Emiss14(sat%nElem, sat%nLine),            &
          clp%Cldphase_Lrc_Mask(sat%nElem, sat%nLine),        &                    	
          clp%Num_Steps_Gradient(sat%nElem, sat%nLine),       &     
          clp%Emiss_Chn10_Tot(sat%nElem, sat%nLine),          &   
          clp%Emiss_Chn10_Tot_Lrc(sat%nElem, sat%nLine),      & 
          clp%Emiss_Chn14_Tot(sat%nElem, sat%nLine),          &
          clp%Emiss_Chn14_Tot_Multi(sat%nElem, sat%nLine),    &        
          clp%Cldbeta7311_Tot(sat%nElem, sat%nLine),          &  
          clp%Cldbeta7311_Tot_Multi(sat%nElem, sat%nLine),    &      	   
          clp%Cldbeta1112_Tot(sat%nElem, sat%nLine),          &  
          clp%Cldbeta1112_Tot_Lrc(sat%nElem, sat%nLine),      & 
          clp%Cldbeta1112_Tot_Multi(sat%nElem, sat%nLine),    &      
          clp%Cldbeta1112_Opaque(sat%nElem, sat%nLine),       &  
          clp%Cldbeta1112_Opaque_Multi(sat%nElem, sat%nLine), &          
          clp%Cldbeta8511_Tot(sat%nElem, sat%nLine),          &  
          clp%Cldbeta8511_Opaque_Lrc(sat%nElem, sat%nLine),   & 
          clp%Cldbeta8511_Tot_Multi(sat%nElem, sat%nLine),    &      
          clp%Cldbeta8511_Opaque(sat%nElem, sat%nLine),       &  
          clp%Cldbeta8511_Opaque_Multi(sat%nElem, sat%nLine), &     
          clp%Opaque_Cld_Temp_Chn10(sat%nElem, sat%nLine),    &  
          clp%Opaque_Cld_Temp_Chn14(sat%nElem, sat%nLine),    &     
          stat=error)
            
if (error /= 0) then
    print *,"(a,'Not enough memory to allocate fylat cloud phase data arrays.')"
    stop
endif

end subroutine allocate_cloudphase_arrays
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~ Subroutine  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine deallocate_cloudphase_arrays()

integer :: error

deallocate (clp%Cldtype,                  &
            clp%Cldtype_Tmpy,             &        
            clp%Cldphase,                 & 
            clp%Cldphase_Qf,              & 
            clp%Cldphase_Qpi,             & 
            clp%Xgrad_Emiss14,            &   
            clp%Ygrad_Emiss14,            &
            clp%Cldphase_Lrc_Mask,        &                    	
            clp%Num_Steps_Gradient,       &     
            clp%Emiss_Chn10_Tot,          &   
            clp%Emiss_Chn10_Tot_Lrc,      & 
            clp%Emiss_Chn14_Tot,          &
            clp%Emiss_Chn14_Tot_Multi,    &        
            clp%Cldbeta7311_Tot,          &  
            clp%Cldbeta7311_Tot_Multi,    &      	   
            clp%Cldbeta1112_Tot,          &  
            clp%Cldbeta1112_Tot_Lrc,      & 
            clp%Cldbeta1112_Tot_Multi,    &      
            clp%Cldbeta1112_Opaque,       &  
            clp%Cldbeta1112_Opaque_Multi, &          
            clp%Cldbeta8511_Tot,          &  
            clp%Cldbeta8511_Opaque_Lrc,   & 
            clp%Cldbeta8511_Tot_Multi,    &      
            clp%Cldbeta8511_Opaque,       &  
            clp%Cldbeta8511_Opaque_Multi, &     
            clp%Opaque_Cld_Temp_Chn10,    &  
            clp%Opaque_Cld_Temp_Chn14,    &     
            stat=error)
            
if (error /= 0) then
    print *,"(ERROR: 'Not enough memory to deallocate fylat cloud phase data arrays.')"
    stop
endif

end subroutine deallocate_cloudphase_arrays
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~ Subroutine  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine allocate_cloudheight_arrays()

integer :: error

allocate (ctp%cldt(sat%nElem, sat%nLine),                 &
          ctp%cldp(sat%nElem, sat%nLine),                 &        
          ctp%cldz(sat%nElem, sat%nLine),                 & 
          ctp%cldemiss(sat%nElem, sat%nLine),             & 
          ctp%cod_vis(sat%nElem, sat%nLine),              & 
          ctp%cldbeta1112(sat%nElem, sat%nLine),          &      
          stat=error)
            
if (error /= 0) then
    print *,"(a,'Not enough memory to allocate fylat cloud height data arrays.')"
    stop
endif

end subroutine allocate_cloudheight_arrays
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~ Subroutine  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine deallocate_cloudheight_arrays()

integer :: error

deallocate (ctp%cldt,                 &
            ctp%cldp,                 &        
            ctp%cldz,                 & 
            ctp%cldemiss,             & 
            ctp%cod_vis,              & 
            ctp%cldbeta1112,          &      
            stat=error)    
            
if (error /= 0) then
    print *,"(ERROR: 'Not enough memory to deallocate fylat cloud height data arrays.')"
    stop
endif

end subroutine deallocate_cloudheight_arrays
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~ Subroutine  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine allocate_cloud_micro_day_arrays()

integer :: error

  real(kind=4), dimension(:,:), pointer :: cod_vis               !cloud temperature (K)
  real(kind=4), dimension(:,:), pointer :: cldreff               !cloud emissivity at 11 microns
  real(kind=4), dimension(:,:), pointer :: cldlwp                !cloud pressure (hPa)
  real(kind=4), dimension(:,:), pointer :: cldiwp                !cloud height (km)
  integer(kind=1), dimension(:,:,:), pointer :: qcflg_cotd         !retrieval quality flags
  
allocate (cot%cod_vis(sat%nElem, sat%nLine),            &
          cot%cldreff(sat%nElem, sat%nLine),            &        
          cot%cldlwp(sat%nElem, sat%nLine),             & 
          cot%cldiwp(sat%nElem, sat%nLine),             & 
          cot%qcflg_cotd(2,sat%nElem, sat%nLine),       &      
          stat=error)
            
if (error /= 0) then
    print *,"(a,'Not enough memory to allocate fylat cloud micro and optical daytime data arrays.')"
    stop
endif

end subroutine allocate_cloud_micro_day_arrays
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~ Subroutine  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine deallocate_cloud_micro_day_arrays()

integer :: error

deallocate (cot%cod_vis,            &
            cot%cldreff,            &        
            cot%cldlwp,             & 
            cot%cldiwp,             & 
            cot%qcflg_cotd,         &      
            stat=error)  
            
if (error /= 0) then
    print *,"(ERROR: 'Not enough memory to deallocate fylat cloud micro and optical daytime data arrays.')"
    stop
endif

end subroutine deallocate_cloud_micro_day_arrays
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!~~~~~~~~~~~~~~~~ Subroutine  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine allocate_sst_arrays()

integer :: error

allocate (sfc%sst(sat%nElem, sat%nLine),                 &
          sfc%qcflg_sst(sat%nElem, sat%nLine),          &      
          stat=error)
            
if (error /= 0) then
    print *,"(a,'Not enough memory to allocate fylat Sea Surface Temperature data arrays.')"
    stop
endif

end subroutine allocate_sst_arrays
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~ Subroutine  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine deallocate_sst_arrays()

integer :: error

deallocate (sfc%sst,                 &
            sfc%qcflg_sst,           &      
            stat=error)    
            
if (error /= 0) then
    print *,"(ERROR: 'Not enough memory to deallocate fylat Sea Surface Temperature data arrays.')"
    stop
endif

end subroutine deallocate_sst_arrays
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~ Subroutine  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine fylat_write_out_surface_sst()

!-----------------------------------------------------------------------
! !F90 write_out_arrays
!
! !Description:
!    The main program for write aerosol detection array out.
!
! !Input  parameters:
!    filename        = 
!
! !Output parameters:
!    none
!-----------------------------------------------------------------------


IMPLICIT NONE

! 1. define variables
! cloud mask
integer :: error, iout
integer (HID_T) :: file_id                 ! File identifier for cloud mask

integer (HID_T) :: sds_id_sst       ! sds id for cloud mask
integer (HID_T) :: sds_id_qcflg_sst

 
!-----
integer (HID_T) :: dsp_id_sst               ! dsp id for cloud optical depth=
integer (HID_T) :: dsp_id_qcflg_sst         ! dsp id for cloud particle effective radius [unit:um]

!-----
!character (LEN=18 ), parameter          :: dset_dust_score_name    = "Aerosol_Dust_Score"         
!character (LEN=17), parameter           :: dset_qcflg_dsc_name     = "Dust_Score_qcflag"  
character (LEN=23), parameter            :: dset_sst_name           = "Sea_Surface_Temperature"
character (LEN=10), parameter            :: dset_qcflg_sst_name     = "SST_qcflag"

integer (HSIZE_T), dimension(2)          :: dims_sst
integer (HSIZE_T), dimension(2)          :: dims_fireRe
integer (HSIZE_T), dimension(2)          :: dims_qf

!integer(kind=4) :: error     ! Error flag for hdf5 open
!integer(kind=4),dimension(:,:),allocatable   :: output_int
integer(kind=4), parameter :: RANK2 = 2 ! Dataset rank
integer(kind=4), parameter :: RANK3 = 3 ! Dataset rank
character, dimension(:,:), allocatable   :: output_int
integer(kind=4), dimension(:,:), allocatable   :: sst_out
!character,dimension(:,:,:),allocatable :: Cldphase_Qpi
character(len=20) :: satellite_name

!---------------------------------------------------
! attribute
integer(HSIZE_T),dimension(1) :: a_dims
integer(kind=4)               :: a_rank
integer(HID_T)                :: output_a_id,output_aspace_id,output_atype_id
! 2. begin program

!+++++++++++++++++++++++++++ Step 2: Write product start +++++++++++++++++++++++++++
!== initialized

!--- variables
allocate( output_int(sat%nElem,sat%nLine),  sst_out(sat%nElem,sat%nLine))

!------         
dims_sst(1) = sat%nElem
dims_sst(2) = sat%nLine


!== 2.1. write start
print*,'    ... fylat write out Sea Surface Temperature (SST) HDF5 product !!! '
	 
call h5open_f(error)

call h5fcreate_f(trim(fy3_mersi_SST_data), H5F_ACC_TRUNC_F, file_id, error)

!------------------------------------
! --- Write Cloud Phase and Type
!call h5dcreate_f(group_id, "EV_1KM_RefSB", H5T_STD_U16LE, dsp_id, sds_id, error)
!call h5dwrite_f (sds_id, H5T_NATIVE_INTEGER, rad(:,:,5:19), dims_sp32, error)
!output_cloudp = char(clp%Cldphase)!cm_bitarray, cm_qa_bitarray
sst_out = int(sfc%sst*100.0)
CALL h5screate_simple_f(RANK2, dims_sst, dsp_id_sst, error)
CALL h5dcreate_f(file_id, dset_sst_name, H5T_STD_U16LE, dsp_id_sst, sds_id_sst, error)
CALL h5dwrite_f (sds_id_sst, H5T_NATIVE_INTEGER, sst_out, dims_sst, error)
! ... write the attributes
      ! .. units attribute
      a_rank = 1
      a_dims = 1
      CALL h5_write_attr_str(sds_id_sst, 'units', 'K')
      !CALL h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      !CALL h5acreate_f(sds_id_sst,'units',H5T_NATIVE_CHARACTER, output_aspace_id, output_a_id, error)
      !call h5awrite_f(output_a_id, H5T_NATIVE_CHARACTER, 'K', a_dims, error)
      !CALL h5sclose_f(output_aspace_id, error)
      !CALL h5aclose_f(output_a_id, error)
      ! .. valid_range attribute
      a_dims = 2
      CALL h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      CALL h5acreate_f(sds_id_sst,'valid_range',H5T_NATIVE_INTEGER, output_aspace_id, output_a_id, error)
      CALL h5awrite_f(output_a_id, H5T_NATIVE_INTEGER, (/0,40000/), a_dims, error)
      CALL h5sclose_f(output_aspace_id, error)
      CALL h5aclose_f(output_a_id, error)
      ! .. _FillValue attibute
      a_dims = 1 
      CALL h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      CALL h5acreate_f(sds_id_sst,'_FillValue',H5T_NATIVE_INTEGER, output_aspace_id, output_a_id, error)
      CALL h5awrite_f(output_a_id, H5T_NATIVE_INTEGER, 0, a_dims, error)
      CALL h5sclose_f(output_aspace_id, error)
      CALL h5aclose_f(output_a_id, error)
      ! .. Intercept attribute
      a_dims = 1 
      CALL h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      CALL h5acreate_f(sds_id_sst,'Intercept',H5T_NATIVE_REAL, output_aspace_id, output_a_id, error)
      CALL h5awrite_f(output_a_id, H5T_NATIVE_REAL, 0.0 , a_dims, error)
      CALL h5sclose_f(output_aspace_id, error)
      CALL h5aclose_f(output_a_id, error)     
      ! .. Slope attribute
      CALL h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      CALL h5acreate_f(sds_id_sst,'Slope',H5T_NATIVE_REAL, output_aspace_id, output_a_id, error)
      CALL h5awrite_f(output_a_id, H5T_NATIVE_REAL, 0.01 , a_dims, error)
      CALL h5sclose_f(output_aspace_id, error)
      CALL h5aclose_f(output_a_id, error)
      ! .. SDS_name
      a_dims = 41
      !CALL h5screate_simple_f(a_rank, a_dims, output_aspace_id, error)
      !CALL h5acreate_f(sds_id_assyocean,'SDS_name',H5T_NATIVE_CHARACTER, output_aspace_id, output_a_id, error)
      !CALL h5awrite_f(output_a_id, H5T_NATIVE_CHARACTER, 'Ocean Aerosol Assymetry Factor at 659nm', a_dims, error)
      !CALL h5sclose_f(output_aspace_id, error)
      !CALL h5aclose_f(output_a_id, error)
      CALL h5_write_attr_str(sds_id_sst, 'SDS_name', 'Sea Surface Temperature')
CALL h5sclose_f (dsp_id_sst, error)
CALL h5dclose_f (sds_id_sst, error)


output_int = char(sfc%qcflg_sst)
call h5screate_simple_f(RANK2, dims_sst, dsp_id_qcflg_sst, error)
call h5dcreate_f(file_id, dset_qcflg_sst_name, H5T_STD_U8LE, dsp_id_qcflg_sst, sds_id_qcflg_sst, error)
call h5dwrite_f (sds_id_qcflg_sst, H5T_STD_U8LE, output_int, dims_sst, error)
call h5dclose_f (sds_id_qcflg_sst, error)      
call h5sclose_f (dsp_id_qcflg_sst, error)



CALL h5screate_simple_f(RANK2, dims_sst, dsp_id_sst, error)
CALL h5dcreate_f(file_id, 'Latitude', H5T_NATIVE_REAL, dsp_id_sst, sds_id_sst, error)
CALL h5dwrite_f (sds_id_sst, H5T_NATIVE_REAL, geo%lat, dims_sst, error)
CALL h5sclose_f (dsp_id_sst, error)
CALL h5dclose_f (sds_id_sst, error)

CALL h5screate_simple_f(RANK2, dims_sst, dsp_id_sst, error)
CALL h5dcreate_f(file_id, 'Longitude', H5T_NATIVE_REAL, dsp_id_sst, sds_id_sst, error)
CALL h5dwrite_f (sds_id_sst, H5T_NATIVE_REAL, geo%lon, dims_sst, error)
CALL h5sclose_f (dsp_id_sst, error)
CALL h5dclose_f (sds_id_sst, error)

deallocate(output_int, sst_out)
           
! close file
call h5fclose_f(file_id, error)

!----------------------------------------------------

call h5close_f(error)


end subroutine fylat_write_out_surface_sst
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE h5_write_attr_str(dset_id, attr_name, attr_data)
! 写入字符型属性
!  USE HDF5
  INTEGER(HID_T) :: dset_id       ! Dataset identifier
  CHARACTER(*)   :: attr_name         ! Attribute name
  CHARACTER(*)   :: attr_data     ! Attribute data
  
  INTEGER(HID_T) :: attr_id       ! Attribute identifier
  INTEGER(HID_T) :: aspace_id     ! Attribute Dataspace identifier
  INTEGER(HID_T) :: atype_id      ! Attribute Dataspace identifier
  INTEGER(HSIZE_T), DIMENSION(1) :: adims = (/1/) ! Attribute dimension
  INTEGER        :: arank = 1                     ! Attribure rank
  INTEGER        :: error         ! Error flag
  INTEGER(SIZE_T) :: attrlen     ! Length of the attribute string
  INTEGER(HSIZE_T), DIMENSION(1) :: data_dims
  data_dims(1) = 1
  
  CALL h5screate_simple_f(arank, adims, aspace_id, error)
  CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, error)
  attrlen = LEN_TRIM(attr_data)
  if (attrlen > 0) then
    CALL h5tset_size_f(atype_id, attrlen, error)
  else
    attrlen = 1
    CALL h5tset_size_f(atype_id, attrlen, error)
  end if
  CALL h5acreate_f(dset_id, attr_name, atype_id, aspace_id, attr_id, error)
  CALL h5awrite_f(attr_id, atype_id, TRIM(attr_data), data_dims, error)
  CALL h5aclose_f(attr_id, error)
  CALL h5sclose_f(aspace_id, error) 
  
END SUBROUTINE h5_write_attr_str
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!-------------------------- END MODULE ---------------------------------
end module io_module


