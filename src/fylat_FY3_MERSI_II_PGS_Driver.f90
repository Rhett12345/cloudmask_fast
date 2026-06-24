program fylat_FY3_MERSI_II_Driver

!-----------------------------------------------------------------------
! f90 driver_program for fylat-fy3_mersi_ii
!
! Description:
!
!    Driver program for FY3/MERSI-II product.
!
! input parameters:
!
! output parameters:
!
!
!-----------------------------------------------------------------------

use names_module
use constant
use data_arrays_module
use get_ancil_data_module
use io_module 
use io_module_intermediate               !lyj 
use read_nwp_data_module
use nwp_utils_module
use frontend_module
use platform_module
use rtm_utils_module
! algorithm
use fylat_fy3mersi_cloud_mask
use fylat_fy3mersi_cloud_amount
!use fylat_fy3mersi_cloud_phase
!use fylat_fy3mersi_cloud_height
!use fylat_fy3mersi_cloud_micro_day
!use fylat_fy3mersi_sst

implicit none

!+++++++++++++++++++ step 1: define variables ++++++++++++++++++++++++++
character(len=1000) :: temp, dummy     ! temporay name   
character(len=1000) :: inCfgFile     ! input namelist File   
!integer(kind=4)     :: i, j, k
integer(kind=4)     :: clock_count_start, clock_count_end
integer(kind=4)     :: clock_rate  
real(kind=4)        :: total_processing_time    
integer(kind=1)     :: write_id


! file name         

namelist /config/ fy3_mersi_L1b_data,   &
                  fy3_mersi_GEO_data,   &
                  fy3_mersi_CLM_data,   &
                  fy3_mersi_CLA_data,   &
                  fy3_mersi_CLP_data,   &
                  fy3_mersi_CTP_data,   &
                  fy3_mersi_COT_data,   &
                  fy3_mersi_CON_data,   &
                  fy3_mersi_SST_data,   &
                  fy3_intermediate,     &       !lyj
                  code_root_path,       &
                  L1b_data_path,        &
                  nwp_data_path,        &
                  oisst_data_path,      &
                  nwp_grib_data1,       &
                  nwp_grib_data2,       &
                  oisst_data,           &
                  fylat_sensor_id,      &
                  fylat_nwp_opt,        &
                  fylat_rtm_opt,        &
                  cloudmask_id,         & ! cloud mask id
                  cloudamount_id,       & ! cloud amount id
                  cloudphase_id,        & ! cloud phase and type id
                  cloudtopz_id,         & ! cloud top height
                  cloudtau_day_id,      & ! cloud optical properties at daytime
                  cloudtau_night_id,    & ! cloud optical properties at nighttime
                  cloudtypeII_id,       & ! cloud type II
                  surface_sst_id,       & ! cloud type II
                  write_inter_id          ! write intermediate result
                  
!---- compute total start time -------
CALL system_clock(count=clock_count_start)

print*, '++++++++++++++++++++++++++++'
print*, '        FYLAT Start         '
print*, '++++++++++++++++++++++++++++'
print*, ' '

CALL getarg(1,inCfgFile)

!%%%%%%%%%%%%%%%
! STEP 0. read configure file
print*, ' '
print*, '-------------------------------------------------------------'
print*,' STEP 0. fylat read configure file and get sensor information '

!.. ... Get The Command_Line as input-Parameter
!open(1,file='fylat_config.nml',status='old')
open(1,file=inCfgFile,status='old')
read(1,NML=config)
close(1)

print*, '-------------------------------------------------------------'
print*, 'code_root_path = ',trim(code_root_path )
print*, ' '
print*, 'Input File:'
print*, 'fy3_mersi_GEO_data = ',trim(fy3_mersi_GEO_data)
print*, 'fy3_mersi_L1b_data = ',trim(fy3_mersi_L1b_data)
print*, 'nwp_grib_data1 = ',trim(nwp_grib_data1)
print*, 'nwp_grib_data2 = ',trim(nwp_grib_data2)
print*, 'oisst_data = ',trim(oisst_data)
print*, '  '
print*, 'Output File:'
print*, 'fy3_mersi_CLM_data = ',trim(fy3_mersi_CLM_data)
print*, 'fy3_mersi_CLA_data = ',trim(fy3_mersi_CLA_data)
print*, 'fy3_mersi_CLP_data = ',trim(fy3_mersi_CLP_data)
print*, 'fy3_mersi_CTP_data = ',trim(fy3_mersi_CTP_data)
print*, 'fy3_mersi_COT_data = ',trim(fy3_mersi_COT_data)
print*, 'fy3_mersi_CON_data = ',trim(fy3_mersi_CON_data)
print*, 'fy3_mersi_SST_data = ',trim(fy3_mersi_SST_data)
print*, 'fy3_intermediate = ',trim(fy3_intermediate)
print*, '  '
print*, '  '
print*, '*****************'
print*, ' fylat Algorithm '
print*, '*****************'

! assign algorithm id 
fylat_alg_opt%cloudmask_index      = cloudmask_id
fylat_alg_opt%cloudamount_index    = cloudamount_id
fylat_alg_opt%cloudphase_index     = cloudphase_id
fylat_alg_opt%cloudtopz_index      = cloudtopz_id
fylat_alg_opt%cloudtau_day_index   = cloudtau_day_id
fylat_alg_opt%cloudtau_night_index = cloudtau_night_id
fylat_alg_opt%cloudtypeII_index    = cloudtypeII_id
fylat_alg_opt%surface_sst_index    = surface_sst_id
write_id                           = write_inter_id

if (fylat_alg_opt%cloudmask_index == 1) then
   print*,'    Cloud Mask Algorithm is Selected!'
endif

if (fylat_alg_opt%cloudamount_index == 1) then
   print*,'    Cloud Amount Algorithm is Selected!'
endif

if (fylat_alg_opt%cloudphase_index == 1) then
   print*,'    Cloud Phase Algorithm is Selected!'
endif

if (fylat_alg_opt%cloudtopz_index == 1) then
   print*,'    Cloud Top Properties Algorithm is Selected!'
endif

if (fylat_alg_opt%cloudtau_day_index == 1) then
   print*,'    Cloud Optical Properties at Daytime Algorithm is Selected!'
endif

if (fylat_alg_opt%cloudtau_night_index == 1) then
   print*,'    Cloud Optical Properties at Nighttime Algorithm is Selected!'
endif

if (fylat_alg_opt%cloudtypeII_index == 1) then
   print*,'    Cloud Type-II Algorithm is Selected!'
endif

if (fylat_alg_opt%surface_sst_index == 1) then
   print*,'    Sea Surface Temperature Algorithm is Selected!'
endif

print*, '  &&& STEP 0. finished '
print*, '   '
!print*, '-----------------------------------------------'

! get element and line information from platform id
call fylat_platform_info(fylat_sensor_id)

! extract satellite time
call extract_sattime(trim(L1b_data_path),trim(fy3_mersi_L1b_data),  &
                     sat%year,sat%month,sat%day,sat%hour,sat%mint,sat%nday)
!print*, sat%year,sat%month,sat%day,sat%hour,sat%mint,sat%nday


!%%%%%%%%%%%%%%%
! STEP 1. read fy3/mersi_II GEO and L1b HDF5 data
print*, '  '
print*, '  '
print*, '------------------------------------------------------'
print*,' STEP 1. fylat read fy3/mersi GEO and L1b HDF5 data '
print*, '------------------------------------------------------'

call fylat_read_fy3_mersi_geo_data(fy3_mersi_GEO_data)

call fylat_read_fy3_mersi_L1b_data(fy3_mersi_L1b_data)

PRINT*,'  ... fylat calculate some angles'
! sun to earth distance 
geo%sun_earth_distance = compute_earth2sun(sat%nday)
 
! get cos satzen and solzen
call compute_cos_zenith_angles(geo%SensorZenith, geo%SolarZenith, &
                               geo%Cos_Satzen, geo%Cos_Solzen)

! get glintzen and scatzen
call compute_scattering_angles(geo%Cos_Satzen, geo%Cos_Solzen,                    &
                               geo%SensorZenith, geo%SolarZenith, geo%RelAzimuth, &
                               geo%Scatzen)
print*,'cos satzen test = ',geo%Cos_Satzen(1000:1003,1000:1003)
print*,'scat angle test = ',geo%Scatzen(1000:1003,1000:1003)
                             
print*, '  &&& STEP 1. finished '

!%%%%%%%%%%%%%%%
! STEP 2. read nwp data
print*, '  '
print*, '  '
print*, '----------------------------------------'
print*,' STEP 2. fylat read and process nwp data '
print*, '----------------------------------------'
  !  fylat_nwp_opt	          
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
  
  
call fylat_read_nwp_data()

call fylat_time_interp_nwp()

call fylat_interp101_nwp()

print*, '  &&& STEP 2. finished '

!%%%%%%%%%%%%%%%
! STEP 3. read fy3/mersi_II ancillary data
print*, '  '
print*, '  '
print*, '-----------------------------------------------'
print*,' STEP 3. fylat read fy3/mersi ancillary data '
print*, '-----------------------------------------------'

call allocate_fylat_ancil_data()

call read_sfc_snow_ice_mask(sat%month)

call read_ecosystem_file(ecosystem_name)

call read_oisst_file()

call read_emissivity_data(sat%month)

call read_albedo_data(sat%nday)

print*, '  &&& STEP 3. finished '
print*, ' '
print*, 'fylat_rtm_opt = ',fylat_rtm_opt

if (fylat_rtm_opt >0) then

   call planck_main(sat%chan_flag)
   
   call ir_rtm_driver()
   
endif

print*, '  &&& STEP 3.1 rtm finished '

!%%%%%%%%%%%%%%%
! STEP 4. run fylat algorithms 
print*, '  '
print*, '  '
print*, '-------------------------------'
print*,' STEP 4. run fylat algorithms  '
print*, '-------------------------------'

!+++ cloud mask +++
if (fylat_alg_opt%cloudmask_index == 1) then 

    ! allocate cloud mask arrays
    allocate(cm_bitarray(sat%nElem,sat%nLine,cm_byte_dim),  &
             cm_qa_bitarray(sat%nElem,sat%nLine,cm_qa_dim))
            
    print*, '  [1] cloud mask algorithm'
    call fy3mersi_cloud_mask()
   
    ! write out cloud mask result
    call fylat_write_out_cloud_mask()
   
    ! deallocate cloud mask output arrays
    !deallocate(cm_bitarray,     &
    !           cm_qa_bitarray)
              
    print*, '  &&& [1] cloud mask algorithm finished '
    print*, ' '
   
endif

if (fylat_alg_opt%cloudmask_index == 1    .and.   &
    (fylat_alg_opt%cloudamount_index == 1 .or.    &
     fylat_alg_opt%cloudphase_index == 1  .or.    &
     fylat_alg_opt%surface_sst_index == 1)) then 
     
    print*, '  [1.1] convert cloud mask product '  
    allocate(cm_tmp(sat%nElem,sat%nLine,2))
    call convert_cloud_mask()
    print*, '  ' 
    
endif


!+++ cloud amount +++
if (fylat_alg_opt%cloudamount_index == 1 .and.  &
    fylat_alg_opt%cloudmask_index == 1 ) then 
	
	ix_5km = int(sat%nElem/5)
	iy_5km = int(sat%nLine/5)
	print*,'  ix_5km, iy_5km = ',ix_5km, iy_5km
	
    ! allocate cloud amount arrays
    allocate(cloud_amount(ix_5km,iy_5km),    &
             cloud_amount_qa(ix_5km,iy_5km), &
             lon_5km(ix_5km,iy_5km),         &
             lat_5km(ix_5km,iy_5km))	
             
    print*,'  [2] cloud amount algorithm'
    call fy3mersi_cloud_amount()

    ! write out cloud amount result
    call fylat_write_out_cloud_amount()
    
    ! deallocate cloud amount arrays
    deallocate(cloud_amount,    &
               cloud_amount_qa, &
               lon_5km,         &
               lat_5km)	

    print*, '  &&& [2] cloud amount algorithm finished '
    print*, ' '
                
endif


!+++ cloud phase and type +++
if (fylat_alg_opt%cloudphase_index > 0 .and.  &   ! set for output 0 = no ; 1 = clp out; 2 = clp + others out
    fylat_alg_opt%cloudmask_index == 1 ) then 

    ! allocate cloud phase arrays
    !call allocate_cloudphase_arrays()
    
    print*,'  [3] cloud phase and type algorithm'
    !call fy3mersi_cloud_phase()

    ! write out cloud phase result
    !call fylat_write_out_cloud_phase()
    
    print*, '  &&& [3] cloud phase algorithm finished '
    print*, ' '
   
endif

!+++ cloud height +++
if (fylat_alg_opt%cloudtopz_index > 0 .and.   &
    fylat_alg_opt%cloudphase_index > 0 .and.  &   ! set for output 0 = no ; 1 = clp out; 2 = clp + others out
    fylat_alg_opt%cloudmask_index == 1 ) then 

    ! allocate cloud phase arrays
    !call allocate_cloudheight_arrays()
    
    print*,'  [4] cloud height algorithm using 11 and 12 um bt'
    !call fy3mersi_cloud_height_11_12()

    ! write out cloud phase result
    !call fylat_write_out_cloud_height()
    
    print*, '  &&& [4] cloud height algorithm using 11 and 12 um bt finished '
    print*, ' '
   
endif

!+++ cloud optical depth at daytime +++
if (fylat_alg_opt%cloudtau_day_index ==1 .and.   &
    fylat_alg_opt%cloudtopz_index > 0 .and.      &
    fylat_alg_opt%cloudphase_index > 0 .and.     &   ! set for output 0 = no ; 1 = clp out; 2 = clp + others out
    fylat_alg_opt%cloudmask_index == 1 ) then 

    ! allocate cloud phase arrays
    !call allocate_cloud_micro_day_arrays()
    
    print*,'  [5] cloud optical depth at daytime algorithm'
    !call fy3mersi_cloud_micro_day()

    ! write out cloud phase result
    !call fylat_write_out_cloud_micro_day()
    
    print*, '  &&& [5] cloud optical depth at daytime algorithm finished '
    print*, ' '
   
endif

!+++ surface SST +++
if (fylat_alg_opt%surface_sst_index > 0 .and.   &
    fylat_alg_opt%cloudmask_index == 1 ) then 

    ! allocate cloud phase arrays
    !call allocate_sst_arrays()
    
    print*,'  [6] SST algorithm '
    !call fy3mersi_sst()

    ! write out cloud phase result
    !call fylat_write_out_surface_sst()
    
    print*, '  &&& [6] SST algorithm finished '
    print*, ' '
   
endif

! --- deallocate output arrays!
if (fylat_alg_opt%cloudmask_index == 1) then 
    ! deallocate cloud mask output arrays
    deallocate(cm_bitarray,     &
               cm_qa_bitarray)   
endif

if (fylat_alg_opt%cloudmask_index == 1 .and.   &
    (fylat_alg_opt%cloudamount_index == 1 .or. &
     fylat_alg_opt%cloudphase_index > 0)) then 
    ! deallocate converted cloud mask output arrays
    deallocate(cm_tmp)
endif

if (fylat_alg_opt%cloudphase_index > 0 .and.  &
    fylat_alg_opt%cloudmask_index == 1 ) then 
    ! deallocate cloud phase output arrays 
    !call deallocate_cloudphase_arrays()
endif

if (fylat_alg_opt%cloudtopz_index > 0 .and.   &
    fylat_alg_opt%cloudphase_index > 0 .and.  &
    fylat_alg_opt%cloudmask_index == 1 ) then 
    ! deallocate cloud height output arrays 
    !call deallocate_cloudheight_arrays()
endif

if (fylat_alg_opt%cloudtau_day_index ==1 .and.   &
    fylat_alg_opt%cloudtopz_index > 0 .and.      &
    fylat_alg_opt%cloudphase_index > 0 .and.     &   ! set for output 0 = no ; 1 = clp out; 2 = clp + others out
    fylat_alg_opt%cloudmask_index == 1 ) then 
    ! deallocate cloud optical depth at daytime arrays
    !call deallocate_cloud_micro_day_arrays()
endif

if (fylat_alg_opt%surface_sst_index > 0 .and.   &
    fylat_alg_opt%cloudmask_index == 1 ) then 
    ! deallocate cloud phase arrays
    !call deallocate_sst_arrays()  
endif


!write_id = 0 
if (write_id == 1) then

    print*, "    ******PRINT INTERMEDIATE ARRAYS START******"        !lyj
    call fylat_write_out_intermediate()                              !lyj
    print*, "    ******PRINT INTERMEDIATE ARRAYS OVER*******"        !lyj
    
endif

!%%%%%%%%%%%%%%%
! STEP XX. deallocate data in memory
print*, '  '
print*, '  '
print*, '-----------------------------------------'
print*,' STEP 5. fylat deallocate data in memory '
print*, '-----------------------------------------'

if (fylat_nwp_opt == 1 .or. fylat_nwp_opt == 2 .or. fylat_nwp_opt == 4 .or. fylat_nwp_opt == 5 ) then 
   call deallocate_nwp26_arrays
endif

if (fylat_nwp_opt == 3 ) then  !T639
   call deallocate_nwp36_arrays
endif

if (fylat_nwp_opt == 6 ) then !grapes
   call deallocate_nwp40_arrays
endif

if (fylat_nwp_opt == 7 .or. fylat_nwp_opt == 8 ) then !gfs 0.25
   call deallocate_nwp31_arrays
endif

if (fylat_nwp_opt == 9 .or. fylat_nwp_opt == 10 ) then !gfs 0.25
   call deallocate_nwp41_arrays
endif



call deallocate_nwp101_arrays() 
call deallocate_fylat_ancil_data()
call deallocate_emiss_arrays()
call deallocate_alb_arrays()
call deallocate_fylat_fy3mersi_geo_data()
call deallocate_fylat_fy3mersi_L1b_data()

print*, '  &&& STEP 5. finished '

print*, '  '
print*,'++++++++++++++++++++++++++++'
print*,'         FYLAT Over         '
print*,'++++++++++++++++++++++++++++'
print*, '  '

!---- show total time cost-------
print*,'       Cost Time         '
call system_clock(count=clock_count_end,count_rate=clock_rate)
total_processing_time = time_elapsed(clock_count_start,clock_count_end,clock_rate,2)

!----------------
end