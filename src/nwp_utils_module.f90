module nwp_utils_module

!C-----------------------------------------------------------------------
!C !F90  nwp_utils_module                                                                
!C
!C !Description: 
!C    This moduleis to deal with Numerical Weather Prediction [nwp] data 
!C    arrays for fylat platform.
!C
!C !Input parameters
!C    none
!C 
!C !Output parameters
!C    none
!C    [all of arrays are stored as public variables 
!C     in data_arrays_module.f90.]
!C
!C !Author's information
!C    Author: Min Min
!C    E-mail: minmin@cma.gov.cn
!C    Tel   : 86-010-68406763
!C    National Satellite Meteorological Center 
!C  
!C !end
!C----------------------------------------------------------------------

! use modules
use data_arrays_module
use numerical
use read_nwp_data_module
use constant
use names_module

implicit none

!+++++++++++++++++++ step 1: define global variables +++++++++++++++++++
integer(kind=4), public, parameter :: NLEVELS_INTERP = 101
integer(kind=4), public, parameter :: PTOP = 25.0
integer(kind=4), public, parameter :: PBOT = 300.0
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


CONTAINS
!+++++++++++++++++++ step 2: subroutines +++++++++++++++++++++++++++++++

!~~~~~~~~~~~~~~~~~~~ subroutine 1: interp_NWP data ~~~~~~~~~~~~~~~~~~~~~
subroutine fylat_time_interp_nwp()

!-----------------------------------------------------------------------
! !F90 interp_NWP
!
! !Description:
!    This program is to interpolate NWP arrays to satellite observation
!    time.
!
! !Input  parameters:
!    none
!
! !Output parameters:
!    none
!
!-----------------------------------------------------------------------

! 1. define variables
real(kind=4)    :: w1, w2, w11, w22
real(kind=4)    :: buf_r4, wo,delta_1, delta_2, delta_p
integer(kind=4)  :: kfirst_rh, wlayer, nw, xn , yn


!*******
! 2. begin program
  print*,'  ... interpolate nwp data to satellite observation time'

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

!=== 2.1. allocate nwp26 arrays
if (fylat_nwp_opt == 1 .or. fylat_nwp_opt == 2 .or. fylat_nwp_opt == 4 .or. fylat_nwp_opt == 5) then   ! 1=ncep and 2=gfs
   call allocate_nwp26_arrays
endif

if (fylat_nwp_opt == 3) then   ! 3=T639
   call allocate_nwp36_arrays
endif

if (fylat_nwp_opt == 6) then   ! 6=grapes gfs 40 layers
   call allocate_nwp40_arrays
endif

IF (fylat_nwp_opt == 7 .or. fylat_nwp_opt == 8) THEN   ! 8 = gfs 0p25
   CALL allocate_nwp31_arrays
ENDIF

IF (fylat_nwp_opt == 9) THEN   ! 9 = gfs 0p50
   CALL allocate_nwp41_arrays
ENDIF

IF (fylat_nwp_opt == 10) THEN   ! 10 = gfs 0p25
   CALL allocate_nwp41_arrays
ENDIF

!=== 2.2. julian data calculation
call cal_julian

!=== 2.3. weight calcultion
call cal_weight(w1,w2)


!=== 2.4. interp start
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

if (fylat_nwp_opt == 1 .or. fylat_nwp_opt == 2 .or. fylat_nwp_opt == 4 .or. fylat_nwp_opt == 5) then   

nwp26%lon  = nwpo%lon(:,:)
nwp26%lat  = nwpo%lat(:,:)

nwp26%wlev = 0.0
kfirst_rh = (nwp26%nlevels-nwp26%nlevels_rh) + 1

xn = nwp26%nlon
yn = nwp26%nlat
if (fylat_nwp_opt == 5) then   
   xn = nwp26%nlon05
   yn = nwp26%nlat05 
endif

do i = 1, xn
do j = 1, yn

   nwp26%plev(i,j,:) = nwpo%plev_nointerp(:,1)
   
   !print*,nwpo%tsfc(i,j,1)
   
   if ((nwpo%psfc(i,j,1) < 0.0) .or. (nwpo%psfc(i,j,2) < 0.0)) then
      nwp26%psfc(i,j) = missing_value_real4
   else
      nwp26%psfc(i,j) = w1 * nwpo%psfc(i,j,1)   + w2 * nwpo%psfc(i,j,2)
   endif
   
   if ((nwpo%pmsl(i,j,1) < 0.0) .or. (nwpo%pmsl(i,j,2) < 0.0)) then
      nwp26%pmsl(i,j) = missing_value_real4
   else
       nwp26%pmsl(i,j) = w1 * nwpo%pmsl(i,j,1)   + w2 * nwpo%pmsl(i,j,2)
   endif

   if ((nwpo%tsfc(i,j,1) < 0.0) .or. (nwpo%tsfc(i,j,2) < 0.0)) then
      nwp26%tsfc(i,j) = missing_value_real4
   else
      nwp26%tsfc(i,j) = w1 * nwpo%tsfc(i,j,1)   + w2 * nwpo%tsfc(i,j,2)
   endif

   if ((nwpo%zsfc(i,j,1) < 0.0) .or. (nwpo%zsfc(i,j,2) < 0.0)) then
      nwp26%zsfc(i,j) = missing_value_real4
   else
      nwp26%zsfc(i,j) = w1 * nwpo%zsfc(i,j,1)   + w2 * nwpo%zsfc(i,j,2)
   endif
   
   if ((nwpo%albedo(i,j,1) < 0.0) .or. (nwpo%albedo(i,j,2) < 0.0)) then
      nwp26%albedo(i,j) = missing_value_real4
   else
      nwp26%albedo(i,j) = w1 * nwpo%albedo(i,j,1)   + w2 * nwpo%albedo(i,j,2)
   endif
   
   if ((nwpo%t_sigma(i,j,1) < 0.0) .or. (nwpo%t_sigma(i,j,2) < 0.0)) then
      nwp26%t_sigma(i,j) = missing_value_real4
   else
      nwp26%t_sigma(i,j) = w1 * nwpo%t_sigma(i,j,1)   + w2 * nwpo%t_sigma(i,j,2)
   endif

   if ((nwpo%rh_sigma(i,j,1) < 0.0) .or. (nwpo%rh_sigma(i,j,2) < 0.0)) then
      nwp26%rh_sigma(i,j) = missing_value_real4
   else
      nwp26%rh_sigma(i,j) = w1 * nwpo%rh_sigma(i,j,1)   + w2 * nwpo%rh_sigma(i,j,2)
   endif

   if ((nwpo%u_sigma(i,j,1) < 0.0) .or. (nwpo%u_sigma(i,j,2) < 0.0)) then
      nwp26%u_sigma(i,j) = missing_value_real4
   else
      nwp26%u_sigma(i,j) = w1 * nwpo%u_sigma(i,j,1)   + w2 * nwpo%u_sigma(i,j,2)
   endif

   if ((nwpo%v_sigma(i,j,1) < 0.0) .or. (nwpo%v_sigma(i,j,2) < 0.0)) then
      nwp26%v_sigma(i,j) = missing_value_real4
   else
      nwp26%v_sigma(i,j) = w1 * nwpo%v_sigma(i,j,1)   + w2 * nwpo%v_sigma(i,j,2)
   endif   

   if ((nwpo%tpw(i,j,1) < 0.0) .or. (nwpo%tpw(i,j,2) < 0.0)) then
      nwp26%tpw(i,j) = missing_value_real4
   else
      nwp26%tpw(i,j) = w1 * nwpo%tpw(i,j,1)   + w2 * nwpo%tpw(i,j,2)
   endif 

   if ((nwpo%weasd(i,j,1) < 0.0) .or. (nwpo%weasd(i,j,2) < 0.0)) then
      nwp26%weasd(i,j) = missing_value_real4
   else
      nwp26%weasd(i,j) = w1 * nwpo%weasd(i,j,1)   + w2 * nwpo%weasd(i,j,2)
   endif    

   if ((nwpo%o3col(i,j,1) < 0.0) .or. (nwpo%o3col(i,j,2) < 0.0)) then
      nwp26%o3col(i,j) = missing_value_real4
   else
      nwp26%o3col(i,j) = w1 * nwpo%o3col(i,j,1)   + w2 * nwpo%o3col(i,j,2)
   endif    
    
   if ((nwpo%ttropo(i,j,1) < 0.0) .or. (nwpo%ttropo(i,j,2) < 0.0)) then
      nwp26%ttropo(i,j) = missing_value_real4
   else
      nwp26%ttropo(i,j) = w1 * nwpo%ttropo(i,j,1)   + w2 * nwpo%ttropo(i,j,2)
   endif  

   !print*,nwpo%tlev(i,j,1:26,1)
   do k=1,nwp26%nlevels


       if ((nwpo%tlev(i,j,k,1) < 0.0) .or. (nwpo%tlev(i,j,k,2) < 0.0)) then
          nwp26%tlev(i,j,k) = missing_value_real4
       else
          nwp26%tlev(i,j,k) = w1 * nwpo%tlev(i,j,k,1)   + w2 * nwpo%tlev(i,j,k,2) 
       endif  

       if ((nwpo%zlev(i,j,k,1) < 0.0) .or. (nwpo%zlev(i,j,k,2) < 0.0)) then
          nwp26%zlev(i,j,k) = missing_value_real4
       else
          nwp26%zlev(i,j,k) = w1 * nwpo%zlev(i,j,k,1)   + w2 * nwpo%zlev(i,j,k,2)
       endif 

       if ((nwpo%o3lev(i,j,k,1) < 0.0) .or. (nwpo%o3lev(i,j,k,2) < 0.0)) then
          nwp26%o3lev(i,j,k) = missing_value_real4
       else
          nwp26%o3lev(i,j,k) = w1 * nwpo%o3lev(i,j,k,1)   + w2 * nwpo%o3lev(i,j,k,2)
       endif               

       if ((nwpo%rhlev(i,j,k,1) < 0.0) .or. (nwpo%rhlev(i,j,k,2) < 0.0)) then
          nwp26%rhlev(i,j,k) = missing_value_real4
       else
          nwp26%rhlev(i,j,k) = w1 * nwpo%rhlev(i,j,k,1)   + w2 * nwpo%rhlev(i,j,k,2)
       endif  
       
       if ((nwpo%clwlev(i,j,k,1) < 0.0) .or. (nwpo%clwlev(i,j,k,2) < 0.0)) then
          nwp26%clwlev(i,j,k) = missing_value_real4
       else
          nwp26%clwlev(i,j,k) = w1 * nwpo%clwlev(i,j,k,1)   + w2 * nwpo%clwlev(i,j,k,2)
       endif    
       

       ! convert rh to water vapor mix ration [g/kg]   
       if (nwp26%rhlev(i,j,k) > 0.0) then
            
            !mpav (3-28-2008) - For the GFS RH is with respect to water/ice
            ! depending on temperature.
            call rh_to_wv(nwp26%rhlev(i,j,k), nwp26%plev(i,j,k), &
                          nwp26%tlev(i,j,k), buf_r4)
                            
            nwp26%wlev(i,j,k) = buf_r4
       
       endif

       if ((nwpo%ulev(i,j,k,1) < -200.0) .or. (nwpo%ulev(i,j,k,2) < -200.0)) then
          nwp26%ulev(i,j,k) = missing_value_real4
       else
          nwp26%ulev(i,j,k) = w1 * nwpo%ulev(i,j,k,1)  + w2 * nwpo%ulev(i,j,k,2)
       endif  
       
       if ((nwpo%vlev(i,j,k,1) < -200.0) .or. (nwpo%vlev(i,j,k,2) < -200.0)) then
          nwp26%vlev(i,j,k) = missing_value_real4
       else
          nwp26%vlev(i,j,k) = w1 * nwpo%vlev(i,j,k,1)  + w2 * nwpo%vlev(i,j,k,2)
       endif

   end do ! k=1, nwp26%nlevels
      
   !------------------------------------------------------------------
   ! Extrapolate NWP profile from top level of RH profile to last 
   ! level of the pressure/temperature profile.
   !------------------------------------------------------------------
        
   wo = max(nwp26%wlev(i,j,kfirst_rh),0.0003)
   do k = kfirst_rh-1, 1, -1
      nwp26%wlev(i,j,k) = max((wo*(nwp26%plev(i,j,k)/nwp26%plev(i,j,kfirst_rh))**3),0.0003)
   end do
   
   do k=1, nwp26%nlevels
      if (nwp26%wlev(i,j,k) <= 0.0) then
          nwp26%wlev(i,j,k) = 0.0003
      endif
   end do
   
   
end do
end do

endif


! T639 nwp36
if (fylat_nwp_opt == 3) then   ! 3=T639

nwp36%lon  = nwpo%lon(:,:)
nwp36%lat  = nwpo%lat(:,:)

nwp36%wlev = 0.0
kfirst_rh = (nwp36%nlevels-nwp36%nlevels_rh) + 1

do i=1, nwp36%nlon
do j=1, nwp36%nlat


   nwp36%plev(i,j,:) = nwpo%plev_nointerp(:,1)
   
   if ((nwpo%psfc(i,j,1) < 0.0) .or. (nwpo%psfc(i,j,2) < 0.0)) then
      nwp36%psfc(i,j) = missing_value_real4
   else
      nwp36%psfc(i,j) = w1 * nwpo%psfc(i,j,1)   + w2 * nwpo%psfc(i,j,2)
   endif
   
   if ((nwpo%pmsl(i,j,1) < 0.0) .or. (nwpo%pmsl(i,j,2) < 0.0)) then
      nwp36%pmsl(i,j) = missing_value_real4
   else
       nwp36%pmsl(i,j) = w1 * nwpo%pmsl(i,j,1)   + w2 * nwpo%pmsl(i,j,2)
   endif

   if ((nwpo%tsfc(i,j,1) < 0.0) .or. (nwpo%tsfc(i,j,2) < 0.0)) then
      nwp36%tsfc(i,j) = missing_value_real4
   else
      nwp36%tsfc(i,j) = w1 * nwpo%tsfc(i,j,1)   + w2 * nwpo%tsfc(i,j,2)
   endif

   if ((nwpo%zsfc(i,j,1) < -10000.0) .or. (nwpo%zsfc(i,j,2) < -10000.0)) then
      nwp36%zsfc(i,j) = missing_value_real4
   else
      nwp36%zsfc(i,j) = w1 * nwpo%zsfc(i,j,1)   + w2 * nwpo%zsfc(i,j,2)
   endif
   
   if ((nwpo%albedo(i,j,1) < 0.0) .or. (nwpo%albedo(i,j,2) < 0.0)) then
      nwp36%albedo(i,j) = missing_value_real4
   else
      nwp36%albedo(i,j) = w1 * nwpo%albedo(i,j,1)   + w2 * nwpo%albedo(i,j,2)
   endif
   
   if ((nwpo%t_sigma(i,j,1) < 0.0) .or. (nwpo%t_sigma(i,j,2) < 0.0)) then
      nwp36%t_sigma(i,j) = missing_value_real4
   else
      nwp36%t_sigma(i,j) = w1 * nwpo%t_sigma(i,j,1) + w2 * nwpo%t_sigma(i,j,2)
   endif

   if ((nwpo%rh_sigma(i,j,1) < 0.0) .or. (nwpo%rh_sigma(i,j,2) < 0.0)) then
      nwp36%rh_sigma(i,j) = missing_value_real4
   else
      nwp36%rh_sigma(i,j) = w1 * nwpo%rh_sigma(i,j,1) + w2 * nwpo%rh_sigma(i,j,2)
   endif

   if ((nwpo%u_sigma(i,j,1) < 0.0) .or. (nwpo%u_sigma(i,j,2) < 0.0)) then
      nwp36%u_sigma(i,j) = missing_value_real4
   else
      nwp36%u_sigma(i,j) = w1 * nwpo%u_sigma(i,j,1) + w2 * nwpo%u_sigma(i,j,2)
   endif

   if ((nwpo%v_sigma(i,j,1) < 0.0) .or. (nwpo%v_sigma(i,j,2) < 0.0)) then
      nwp36%v_sigma(i,j) = missing_value_real4
   else
      nwp36%v_sigma(i,j) = w1 * nwpo%v_sigma(i,j,1) + w2 * nwpo%v_sigma(i,j,2)
   endif   

   if ((nwpo%tpw(i,j,1) < 0.0) .or. (nwpo%tpw(i,j,2) < 0.0)) then
      nwp36%tpw(i,j) = missing_value_real4
   else
      nwp36%tpw(i,j) = w1 * nwpo%tpw(i,j,1) + w2 * nwpo%tpw(i,j,2)
   endif 

   if ((nwpo%weasd(i,j,1) < 0.0) .or. (nwpo%weasd(i,j,2) < 0.0)) then
      nwp36%weasd(i,j) = missing_value_real4
   else
      nwp36%weasd(i,j) = w1 * nwpo%weasd(i,j,1) + w2 * nwpo%weasd(i,j,2)
   endif    

   if ((nwpo%o3col(i,j,1) < 0.0) .or. (nwpo%o3col(i,j,2) < 0.0)) then
      nwp36%o3col(i,j) = missing_value_real4
   else
      nwp36%o3col(i,j) = w1 * nwpo%o3col(i,j,1) + w2 * nwpo%o3col(i,j,2)
   endif    
    
   if ((nwpo%ttropo(i,j,1) < 0.0) .or. (nwpo%ttropo(i,j,2) < 0.0)) then
      nwp36%ttropo(i,j) = missing_value_real4
   else
      nwp36%ttropo(i,j) = w1 * nwpo%ttropo(i,j,1) + w2 * nwpo%ttropo(i,j,2)
   endif  


   wlayer = 1

   nw = 1

   do k=1,nwp36%nlevels

   if ( nwpo%tlev(i,j,k,1) < 400. .and. nwpo%tlev(i,j,k,2) < 400.0) then  ! good layer

       wlayer = k

       if ((nwpo%tlev(i,j,k,1) < 0.0) .or. (nwpo%tlev(i,j,k,2) < 0.0)) then
          nwp36%tlev(i,j,k) = missing_value_real4
       else
          nwp36%tlev(i,j,k) = w1 * nwpo%tlev(i,j,k,1)   + w2 * nwpo%tlev(i,j,k,2) 
       endif  

       if ((nwpo%zlev(i,j,k,1) < 0.0) .or. (nwpo%zlev(i,j,k,2) < 0.0)) then
          nwp36%zlev(i,j,k) = missing_value_real4
       else
          nwp36%zlev(i,j,k) = w1 * nwpo%zlev(i,j,k,1)   + w2 * nwpo%zlev(i,j,k,2)
       endif 

       if ((nwpo%o3lev(i,j,k,1) < 0.0) .or. (nwpo%o3lev(i,j,k,2) < 0.0)) then
          nwp36%o3lev(i,j,k) = missing_value_real4
       else
          nwp36%o3lev(i,j,k) = w1 * nwpo%o3lev(i,j,k,1)   + w2 * nwpo%o3lev(i,j,k,2)
       endif               

       if ((nwpo%rhlev(i,j,k,1) < 0.0) .or. (nwpo%rhlev(i,j,k,2) < 0.0)) then
          nwp36%rhlev(i,j,k) = missing_value_real4
       else
          nwp36%rhlev(i,j,k) = w1 * nwpo%rhlev(i,j,k,1)   + w2 * nwpo%rhlev(i,j,k,2)
       endif  
       
       if ((nwpo%clwlev(i,j,k,1) < 0.0) .or. (nwpo%clwlev(i,j,k,2) < 0.0)) then
          nwp36%clwlev(i,j,k) = missing_value_real4
       else
          nwp36%clwlev(i,j,k) = w1 * nwpo%clwlev(i,j,k,1)   + w2 * nwpo%clwlev(i,j,k,2)
       endif    
       

       ! convert rh to water vapor mix ration [g/kg]   
       if (nwp36%rhlev(i,j,k) > 0.0) then
            
            !mpav (3-28-2008) - For the GFS RH is with respect to water/ice
            ! depending on temperature.
            call rh_to_wv(nwp36%rhlev(i,j,k), nwp36%plev(i,j,k), &
                          nwp36%tlev(i,j,k), buf_r4)
                            
            nwp36%wlev(i,j,k) = buf_r4
       
       endif

       IF ((nwpo%ulev(i,j,k,1) < -200.0) .or. (nwpo%ulev(i,j,k,2) < -200.0)) THEN
          nwp36%ulev(i,j,k) = missing_value_real4
       ELSE
          nwp36%ulev(i,j,k) = w1 * nwpo%ulev(i,j,k,1)  + w2 * nwpo%ulev(i,j,k,2)
       ENDIF  
       
       IF ((nwpo%vlev(i,j,k,1) < -200.0) .or. (nwpo%vlev(i,j,k,2) < -200.0)) THEN
          nwp36%vlev(i,j,k) = missing_value_real4
       ELSE
          nwp36%vlev(i,j,k) = w1 * nwpo%vlev(i,j,k,1)  + w2 * nwpo%vlev(i,j,k,2)
       ENDIF 

   else  ! bad layer

       delta_p = log (nwp36%plev(i,j,k) / nwp36%plev(i,j,wlayer))

       delta_1 = (nwpo%tlev(i,j,wlayer,1)-nwpo%tlev(i,j,wlayer-1,1)) / log (nwp36%plev(i,j,wlayer) / nwp36%plev(i,j,wlayer-1))
       delta_2 = (nwpo%tlev(i,j,wlayer,2)-nwpo%tlev(i,j,wlayer-1,2)) / log (nwp36%plev(i,j,wlayer) / nwp36%plev(i,j,wlayer-1))
       nwp36%tlev(i,j,k) = w1 * (nwpo%tlev(i,j,wlayer,1) + delta_1*delta_p)    + w2 * (nwpo%tlev(i,j,wlayer,2) + delta_2*delta_p)

       delta_1 = (nwpo%zlev(i,j,wlayer,1)-nwpo%zlev(i,j,wlayer-1,1)) / log (nwp36%plev(i,j,wlayer) / nwp36%plev(i,j,wlayer-1))
       delta_2 = (nwpo%zlev(i,j,wlayer,2)-nwpo%zlev(i,j,wlayer-1,2)) / log (nwp36%plev(i,j,wlayer) / nwp36%plev(i,j,wlayer-1))
       nwp36%zlev(i,j,k) = w1 * (nwpo%zlev(i,j,wlayer,1) + delta_1*delta_p)    + w2 * (nwpo%zlev(i,j,wlayer,2) + delta_2*delta_p)

       delta_1 = (nwpo%rhlev(i,j,wlayer,1)-nwpo%rhlev(i,j,wlayer-1,1)) / log (nwp36%plev(i,j,wlayer) / nwp36%plev(i,j,wlayer-1))
       delta_2 = (nwpo%rhlev(i,j,wlayer,2)-nwpo%rhlev(i,j,wlayer-1,2)) / log (nwp36%plev(i,j,wlayer) / nwp36%plev(i,j,wlayer-1))
       nwp36%rhlev(i,j,k) = w1 * (nwpo%rhlev(i,j,wlayer,1) + delta_1*delta_p)    + w2 * (nwpo%rhlev(i,j,wlayer,2) + delta_2*delta_p)

       nw = nw + 1

       ! convert rh to water vapor mix ration [g/kg]
       if (nwp36%rhlev(i,j,k) > 0.0) then

           !mpav (3-28-2008) - For the GFS RH is with respect to water/ice
           ! depending on temperature.
           call rh_to_wv(nwp36%rhlev(i,j,k), nwp36%plev(i,j,k), &
                         nwp36%tlev(i,j,k), buf_r4)

           nwp36%wlev(i,j,k) = buf_r4

       endif

   endif

   end do ! k=1, nwp26%nlevels

   !------------------------------------------------------------------
   ! Extrapolate NWP profile from top level of RH profile to last
   ! level of the pressure/temperature profile.
   !------------------------------------------------------------------
        
   wo = max(nwp36%wlev(i,j,kfirst_rh),0.0003)
   do k = kfirst_rh-1, 1, -1
      nwp36%wlev(i,j,k) = max((wo*(nwp36%plev(i,j,k)/nwp36%plev(i,j,kfirst_rh))**3),0.0003)
   end do
   
   do k=1, nwp36%nlevels
      if (nwp36%wlev(i,j,k) <= 0.0) then
          nwp36%wlev(i,j,k) = 0.0003
      endif
   end do
   
end do
end do

endif


! grapes gfs 0.25*0.25
IF (fylat_nwp_opt == 6) THEN   ! 9=grapes gfs 0.25*0.25 grib2

nwp40%lon  = nwpo%lon(:,:)
nwp40%lat  = nwpo%lat(:,:)

nwp40%wlev = 0.0
kfirst_rh = (nwp40%nlevels-nwp40%nlevels_rh) + 1
print*,'w1,w2 =',w1,w2
DO i = 1, nwp40%nlon
DO j = 1, nwp40%nlat


   nwp40%plev(i,j,:) = nwpo%plev_nointerp(:,1)
   
   IF ((nwpo%psfc(i,j,1) < 0.0) .or. (nwpo%psfc(i,j,2) < 0.0)) THEN
      nwp40%psfc(i,j) = missing_value_real4
   ELSE
      nwp40%psfc(i,j) = w1 * nwpo%psfc(i,j,1)   + w2 * nwpo%psfc(i,j,2)
   ENDIF
   
   IF ((nwpo%pmsl(i,j,1) < 0.0) .or. (nwpo%pmsl(i,j,2) < 0.0)) THEN
      nwp40%pmsl(i,j) = missing_value_real4
   ELSE
       nwp40%pmsl(i,j) = w1 * nwpo%pmsl(i,j,1)   + w2 * nwpo%pmsl(i,j,2)
   ENDIF

   IF ((nwpo%tsfc(i,j,1) < 0.0) .or. (nwpo%tsfc(i,j,2) < 0.0)) THEN
      nwp40%tsfc(i,j) = missing_value_real4
   ELSE
      nwp40%tsfc(i,j) = w1 * nwpo%tsfc(i,j,1)   + w2 * nwpo%tsfc(i,j,2)
   ENDIF

   IF ((nwpo%zsfc(i,j,1) < -10000.0) .or. (nwpo%zsfc(i,j,2) < -10000.0)) THEN
      nwp40%zsfc(i,j) = missing_value_real4
   ELSE
      nwp40%zsfc(i,j) = w1 * nwpo%zsfc(i,j,1)   + w2 * nwpo%zsfc(i,j,2)
   ENDIF
   
   IF ((nwpo%albedo(i,j,1) < 0.0) .or. (nwpo%albedo(i,j,2) < 0.0)) THEN
      nwp40%albedo(i,j) = missing_value_real4
   ELSE
      nwp40%albedo(i,j) = w1 * nwpo%albedo(i,j,1)   + w2 * nwpo%albedo(i,j,2)
   ENDIF
   
   IF ((nwpo%t_sigma(i,j,1) < 0.0) .or. (nwpo%t_sigma(i,j,2) < 0.0)) THEN
      nwp40%t_sigma(i,j) = missing_value_real4
   ELSE
      nwp40%t_sigma(i,j) = w1 * nwpo%t_sigma(i,j,1) + w2 * nwpo%t_sigma(i,j,2)
   ENDIF

   IF ((nwpo%rh_sigma(i,j,1) < 0.0) .or. (nwpo%rh_sigma(i,j,2) < 0.0)) THEN
      nwp40%rh_sigma(i,j) = missing_value_real4
   ELSE
      nwp40%rh_sigma(i,j) = w1 * nwpo%rh_sigma(i,j,1) + w2 * nwpo%rh_sigma(i,j,2)
   ENDIF

   IF ((nwpo%u_sigma(i,j,1) < 0.0) .or. (nwpo%u_sigma(i,j,2) < 0.0)) THEN
      nwp40%u_sigma(i,j) = missing_value_real4
   ELSE
      nwp40%u_sigma(i,j) = w1 * nwpo%u_sigma(i,j,1) + w2 * nwpo%u_sigma(i,j,2)
   ENDIF

   IF ((nwpo%v_sigma(i,j,1) < 0.0) .or. (nwpo%v_sigma(i,j,2) < 0.0)) THEN
      nwp40%v_sigma(i,j) = missing_value_real4
   ELSE
      nwp40%v_sigma(i,j) = w1 * nwpo%v_sigma(i,j,1) + w2 * nwpo%v_sigma(i,j,2)
   ENDIF   

   IF ((nwpo%tpw(i,j,1) < 0.0) .or. (nwpo%tpw(i,j,2) < 0.0)) THEN
      nwp40%tpw(i,j) = missing_value_real4
   ELSE
      nwp40%tpw(i,j) = w1 * nwpo%tpw(i,j,1) + w2 * nwpo%tpw(i,j,2)
   ENDIF 

   IF ((nwpo%weasd(i,j,1) < 0.0) .or. (nwpo%weasd(i,j,2) < 0.0)) THEN
      nwp40%weasd(i,j) = missing_value_real4
   ELSE
      nwp40%weasd(i,j) = w1 * nwpo%weasd(i,j,1) + w2 * nwpo%weasd(i,j,2)
   ENDIF    

   IF ((nwpo%o3col(i,j,1) < 0.0) .or. (nwpo%o3col(i,j,2) < 0.0)) THEN
      nwp40%o3col(i,j) = missing_value_real4
   ELSE
      nwp40%o3col(i,j) = w1 * nwpo%o3col(i,j,1) + w2 * nwpo%o3col(i,j,2)
   ENDIF    
    
   IF ((nwpo%ttropo(i,j,1) < 0.0) .or. (nwpo%ttropo(i,j,2) < 0.0)) THEN
      nwp40%ttropo(i,j) = missing_value_real4
   ELSE
      nwp40%ttropo(i,j) = w1 * nwpo%ttropo(i,j,1) + w2 * nwpo%ttropo(i,j,2)
   ENDIF  

 !  IF ((nwpo%u10m(i,j,1) < -200.0) .or. (nwpo%u10m(i,j,2) < -200.0)) THEN
 !     nwp40%u10m(i,j) = missing_value_real4
 !  ELSE
 !     nwp40%u10m(i,j) = w1 * nwpo%u10m(i,j,1)   + w2 * nwpo%u10m(i,j,2)
 !  ENDIF    

 !  IF ((nwpo%v10m(i,j,1) < -200.0) .or. (nwpo%v10m(i,j,2) < -200.0)) THEN
 !     nwp40%v10m(i,j) = missing_value_real4
 !  ELSE
 !     nwp40%v10m(i,j) = w1 * nwpo%v10m(i,j,1)  + w2 * nwpo%v10m(i,j,2)
 !  ENDIF 

   wlayer = 1

   nw = 1

   DO k=1,nwp40%nlevels
      
      w11 = w1
      w22 = w2
      
      IF ( (nwpo%tlev(i,j,k,1) > 400. .or. nwpo%tlev(i,j,k,1) < 120.) .and. nwpo%tlev(i,j,k,2) < 400.) THEN  ! good layer
          w11 = 0.0
          w22 = 1.0
      ENDIF

      IF ( nwpo%tlev(i,j,k,1) < 400. .and. (nwpo%tlev(i,j,k,2) > 400. .or. nwpo%tlev(i,j,k,2) < 120.)) THEN  ! good layer
          w11 = 1.0
          w22 = 0.0
      ENDIF
      

       wlayer = k

       IF ((nwpo%tlev(i,j,k,1) < 0.0) .or. (nwpo%tlev(i,j,k,2) < 0.0)) THEN
          nwp40%tlev(i,j,k) = missing_value_real4
       ELSE
          nwp40%tlev(i,j,k) = w11 * nwpo%tlev(i,j,k,1)   + w22 * nwpo%tlev(i,j,k,2) 
       ENDIF  

       IF ((nwpo%zlev(i,j,k,1) < 0.0) .or. (nwpo%zlev(i,j,k,2) < 0.0)) THEN
          nwp40%zlev(i,j,k) = missing_value_real4
       ELSE
          nwp40%zlev(i,j,k) = w11 * nwpo%zlev(i,j,k,1)   + w22 * nwpo%zlev(i,j,k,2)
       ENDIF 

       IF ((nwpo%o3lev(i,j,k,1) < 0.0) .or. (nwpo%o3lev(i,j,k,2) < 0.0)) THEN
          nwp40%o3lev(i,j,k) = missing_value_real4
       ELSE
          nwp40%o3lev(i,j,k) = w11 * nwpo%o3lev(i,j,k,1)   + w22 * nwpo%o3lev(i,j,k,2)
       ENDIF               

       IF ((nwpo%rhlev(i,j,k,1) < 0.0) .or. (nwpo%rhlev(i,j,k,2) < 0.0)) THEN
          nwp40%rhlev(i,j,k) = missing_value_real4
       ELSE
          nwp40%rhlev(i,j,k) = w11 * nwpo%rhlev(i,j,k,1)   + w22 * nwpo%rhlev(i,j,k,2)
       ENDIF  
       
       IF ((nwpo%clwlev(i,j,k,1) < 0.0) .or. (nwpo%clwlev(i,j,k,2) < 0.0)) THEN
          nwp40%clwlev(i,j,k) = missing_value_real4
       ELSE
          nwp40%clwlev(i,j,k) = w11 * nwpo%clwlev(i,j,k,1)   + w22 * nwpo%clwlev(i,j,k,2)
       ENDIF    
       

       ! convert rh to water vapor mix ration [g/kg]   
       IF (nwp40%rhlev(i,j,k) > 0.0) THEN
            
            !mpav (3-28-2008) - For the GFS RH is with respect to water/ice
            ! depending on temperature.
            CALL rh_to_wv(nwp40%rhlev(i,j,k), nwp40%plev(i,j,k), &
                          nwp40%tlev(i,j,k), buf_r4)
                            
            nwp40%wlev(i,j,k) = buf_r4
       
       ENDIF
       
       IF ((nwpo%ulev(i,j,k,1) < -200.0) .or. (nwpo%ulev(i,j,k,2) < -200.0)) THEN
          nwp40%ulev(i,j,k) = missing_value_real4
       ELSE
          nwp40%ulev(i,j,k) = w1 * nwpo%ulev(i,j,k,1)  + w2 * nwpo%ulev(i,j,k,2)
       ENDIF  
       
       IF ((nwpo%vlev(i,j,k,1) < -200.0) .or. (nwpo%vlev(i,j,k,2) < -200.0)) THEN
          nwp40%vlev(i,j,k) = missing_value_real4
       ELSE
          nwp40%vlev(i,j,k) = w1 * nwpo%vlev(i,j,k,1)  + w2 * nwpo%vlev(i,j,k,2)
       ENDIF 

   !ENDIF

   END DO ! k=1, nwp40%nlevels

   !------------------------------------------------------------------
   ! Extrapolate NWP profile from top level of RH profile to last
   ! level of the pressure/temperature profile.
   !------------------------------------------------------------------
        
   wo = max(nwp40%wlev(i,j,kfirst_rh),0.0003)
   DO k = kfirst_rh-1, 1, -1
      nwp40%wlev(i,j,k) = max((wo*(nwp40%plev(i,j,k)/nwp40%plev(i,j,kfirst_rh))**3),0.0003)
   END DO
   
   !DO k=1, nwp40%nlevels
   !   IF (nwp40%wlev(i,j,k) <= 0.0) THEN
   !       nwp40%wlev(i,j,k) = 0.0003
   !   ENDIF
   !END DO
   DO k=1, nwp40%nlevels
      IF (k <= 14) THEN
          nwp40%wlev(i,j,k) = 0.0003
      ENDIF
      !IF (nwp40%wlev(i,j,k) <= 0.0) THEN
      IF (nwp40%wlev(i,j,k) <= 0.0003) THEN
          nwp40%wlev(i,j,k) = 0.0003
      ENDIF
   END DO
   
END DO
END DO

ENDIF ! nwp_opt = 6 grapes gfs 


IF (fylat_nwp_opt == 7 .or. fylat_nwp_opt == 8) THEN   ! 8=gfs 0.25 grib2

nwp31%lon  = nwpo%lon(:,:)
nwp31%lat  = nwpo%lat(:,:)

nwp31%wlev = 0.0
kfirst_rh = (nwp31%nlevels-nwp31%nlevels_rh) + 1

IF (fylat_nwp_opt == 7 .or. fylat_nwp_opt == 8) THEN
   xn = nwp31%nlon25
   yn = nwp31%nlat25
ENDIF
print*,'w1,w2 =',w1,w2

DO i = 1, xn !nwp26%nlon
DO j = 1, yn !nwp26%nlat

   nwp31%plev(i,j,:) = nwpo%plev_nointerp(:,1)
   
   !print*,'t1 = ',i,j,xn,yn,w1,w2,nwpo%tsfc(i,j,1:2),nwpo%psfc(i,j,1:2),nwpo%zsfc(i,j,1:2) 
  
   IF ((nwpo%psfc(i,j,1) < 0.0) .or. (nwpo%psfc(i,j,2) < 0.0)) THEN
      nwp31%psfc(i,j) = missing_value_real4
   ELSE
      nwp31%psfc(i,j) = w1 * nwpo%psfc(i,j,1)   + w2 * nwpo%psfc(i,j,2)
   ENDIF
 
   IF ((nwpo%pmsl(i,j,1) < 0.0) .or. (nwpo%pmsl(i,j,2) < 0.0)) THEN
      nwp31%pmsl(i,j) = missing_value_real4
   ELSE
      nwp31%pmsl(i,j) = w1 * nwpo%pmsl(i,j,1)   + w2 * nwpo%pmsl(i,j,2)
   ENDIF

   IF ((nwpo%tsfc(i,j,1) < 0.0) .or. (nwpo%tsfc(i,j,2) < 0.0)) THEN
      nwp31%tsfc(i,j) = missing_value_real4
   ELSE
      nwp31%tsfc(i,j) = w1 * nwpo%tsfc(i,j,1)   + w2 * nwpo%tsfc(i,j,2)
   ENDIF

   IF ((nwpo%zsfc(i,j,1) < 0.0) .or. (nwpo%zsfc(i,j,2) < 0.0)) THEN
      nwp31%zsfc(i,j) = missing_value_real4
   ELSE
      nwp31%zsfc(i,j) = w1 * nwpo%zsfc(i,j,1)   + w2 * nwpo%zsfc(i,j,2)
   ENDIF
 
   IF ((nwpo%albedo(i,j,1) < 0.0) .or. (nwpo%albedo(i,j,2) < 0.0)) THEN
      nwp31%albedo(i,j) = missing_value_real4
   ELSE
      nwp31%albedo(i,j) = w1 * nwpo%albedo(i,j,1)   + w2 * nwpo%albedo(i,j,2)
   ENDIF
   
   IF ((nwpo%t_sigma(i,j,1) < 0.0) .or. (nwpo%t_sigma(i,j,2) < 0.0)) THEN
      nwp31%t_sigma(i,j) = missing_value_real4
   ELSE
      nwp31%t_sigma(i,j) = w1 * nwpo%t_sigma(i,j,1)   + w2 * nwpo%t_sigma(i,j,2)
   ENDIF

   IF ((nwpo%rh_sigma(i,j,1) < 0.0) .or. (nwpo%rh_sigma(i,j,2) < 0.0)) THEN
      nwp31%rh_sigma(i,j) = missing_value_real4
   ELSE
      nwp31%rh_sigma(i,j) = w1 * nwpo%rh_sigma(i,j,1)   + w2 * nwpo%rh_sigma(i,j,2)
   ENDIF

   IF ((nwpo%u_sigma(i,j,1) < 0.0) .or. (nwpo%u_sigma(i,j,2) < 0.0)) THEN
      nwp31%u_sigma(i,j) = missing_value_real4
   ELSE
      nwp31%u_sigma(i,j) = w1 * nwpo%u_sigma(i,j,1)   + w2 * nwpo%u_sigma(i,j,2)
   ENDIF

   IF ((nwpo%v_sigma(i,j,1) < 0.0) .or. (nwpo%v_sigma(i,j,2) < 0.0)) THEN
      nwp31%v_sigma(i,j) = missing_value_real4
   ELSE
      nwp31%v_sigma(i,j) = w1 * nwpo%v_sigma(i,j,1)   + w2 * nwpo%v_sigma(i,j,2)
   ENDIF   

   IF ((nwpo%tpw(i,j,1) < 0.0) .or. (nwpo%tpw(i,j,2) < 0.0)) THEN
      nwp31%tpw(i,j) = missing_value_real4
   ELSE
      nwp31%tpw(i,j) = w1 * nwpo%tpw(i,j,1)   + w2 * nwpo%tpw(i,j,2)
   ENDIF 

   IF ((nwpo%weasd(i,j,1) < 0.0) .or. (nwpo%weasd(i,j,2) < 0.0)) THEN
      nwp31%weasd(i,j) = missing_value_real4
   ELSE
      nwp31%weasd(i,j) = w1 * nwpo%weasd(i,j,1)   + w2 * nwpo%weasd(i,j,2)
   ENDIF    

   IF ((nwpo%o3col(i,j,1) < 0.0) .or. (nwpo%o3col(i,j,2) < 0.0)) THEN
      nwp31%o3col(i,j) = missing_value_real4
   ELSE
      nwp31%o3col(i,j) = w1 * nwpo%o3col(i,j,1)   + w2 * nwpo%o3col(i,j,2)
   ENDIF    
    
   IF ((nwpo%ttropo(i,j,1) < 0.0) .or. (nwpo%ttropo(i,j,2) < 0.0)) THEN
      nwp31%ttropo(i,j) = missing_value_real4
   ELSE
      nwp31%ttropo(i,j) = w1 * nwpo%ttropo(i,j,1)   + w2 * nwpo%ttropo(i,j,2)
   ENDIF  
   !print*,'ttrop',nwpo%ttropo(i,j,1) ,nwpo%ttropo(i,j,2),nwpo%o3col(i,j,1:2),nwpo%tpw(i,j,1:2) 
 
 
   !print*,nwpo%tlev(i,j,1:26,1)
   DO k=1,nwp31%nlevels


       IF ((nwpo%tlev(i,j,k,1) < 0.0) .or. (nwpo%tlev(i,j,k,2) < 0.0)) THEN
          nwp31%tlev(i,j,k) = missing_value_real4
       ELSE
          nwp31%tlev(i,j,k) = w1 * nwpo%tlev(i,j,k,1)   + w2 * nwpo%tlev(i,j,k,2) 
       ENDIF  

       IF ((nwpo%zlev(i,j,k,1) < 0.0) .or. (nwpo%zlev(i,j,k,2) < 0.0)) THEN
          nwp31%zlev(i,j,k) = missing_value_real4
       ELSE
          nwp31%zlev(i,j,k) = w1 * nwpo%zlev(i,j,k,1)   + w2 * nwpo%zlev(i,j,k,2)
       ENDIF 

       IF ((nwpo%o3lev(i,j,k,1) < 0.0) .or. (nwpo%o3lev(i,j,k,2) < 0.0)) THEN
          nwp31%o3lev(i,j,k) = missing_value_real4
       ELSE
          nwp31%o3lev(i,j,k) = w1 * nwpo%o3lev(i,j,k,1)   + w2 * nwpo%o3lev(i,j,k,2)
       ENDIF               

       IF ((nwpo%rhlev(i,j,k,1) < 0.0) .or. (nwpo%rhlev(i,j,k,2) < 0.0)) THEN
          nwp31%rhlev(i,j,k) = missing_value_real4
       ELSE
          nwp31%rhlev(i,j,k) = w1 * nwpo%rhlev(i,j,k,1)   + w2 * nwpo%rhlev(i,j,k,2)
       ENDIF  
       
       IF ((nwpo%clwlev(i,j,k,1) < 0.0) .or. (nwpo%clwlev(i,j,k,2) < 0.0)) THEN
          nwp31%clwlev(i,j,k) = missing_value_real4
       ELSE
          nwp31%clwlev(i,j,k) = w1 * nwpo%clwlev(i,j,k,1)   + w2 * nwpo%clwlev(i,j,k,2)
       ENDIF    
       

       ! convert rh to water vapor mix ration [g/kg]   
       IF (nwp31%rhlev(i,j,k) > 0.0) THEN
            
            !mpav (3-28-2008) - For the GFS RH is with respect to water/ice
            ! depending on temperature.
            CALL rh_to_wv(nwp31%rhlev(i,j,k), nwp31%plev(i,j,k), &
                          nwp31%tlev(i,j,k), buf_r4)
                            
            nwp31%wlev(i,j,k) = buf_r4
       
       ENDIF
          !print*,'t2 = ',i,j,k,nwpo%tlev(i,j,k,1:2),nwpo%rhlev(i,j,k,1:2),nwpo%zlev(i,j,k,1:2)

       
      IF ((nwpo%ulev(i,j,k,1) < -200.0) .or. (nwpo%ulev(i,j,k,2) < -200.0)) THEN
         nwp31%ulev(i,j,k) = missing_value_real4
      ELSE
         nwp31%ulev(i,j,k) = w1 * nwpo%ulev(i,j,k,1)  + w2 * nwpo%ulev(i,j,k,2)
      ENDIF  
       
      IF ((nwpo%vlev(i,j,k,1) < -200.0) .or. (nwpo%vlev(i,j,k,2) < -200.0)) THEN
         nwp31%vlev(i,j,k) = missing_value_real4
      ELSE
         nwp31%vlev(i,j,k) = w1 * nwpo%vlev(i,j,k,1)  + w2 * nwpo%vlev(i,j,k,2)
      ENDIF 
         
   END DO ! k=1, nwp26%nlevels
       
   !------------------------------------------------------------------
   ! Extrapolate NWP profile from top level of RH profile to last 
   ! level of the pressure/temperature profile.
   !------------------------------------------------------------------
        
   wo = max(nwp31%wlev(i,j,kfirst_rh),0.0003)
   DO k = kfirst_rh-1, 1, -1
      nwp31%wlev(i,j,k) = max((wo*(nwp31%plev(i,j,k)/nwp31%plev(i,j,kfirst_rh))**3),0.0003)
   END DO
   
   DO k=1, nwp31%nlevels
      IF (nwp31%wlev(i,j,k) <= 0.0) THEN
          nwp31%wlev(i,j,k) = 0.0003
      ENDIF
   END DO
    
   
END DO
END DO

ENDIF


IF (fylat_nwp_opt == 9 .or. fylat_nwp_opt == 10) THEN   ! 8=gfs 0.25 grib2

nwp41%lon  = nwpo%lon(:,:)
nwp41%lat  = nwpo%lat(:,:)

nwp41%wlev = 0.0
kfirst_rh = (nwp41%nlevels-nwp41%nlevels_rh) + 1

IF (fylat_nwp_opt == 9) THEN
   xn = nwp41%nlon05
   yn = nwp41%nlat05
ENDIF

IF (fylat_nwp_opt == 10) THEN
   xn = nwp41%nlon25
   yn = nwp41%nlat25
ENDIF


DO i = 1, xn !nwp41%nlon
DO j = 1, yn !nwp41%nlat

   nwp41%plev(i,j,:) = nwpo%plev_nointerp(:,1)
   
   !print*,'t1 = ',i,j,xn,yn,w1,w2,nwpo%tsfc(i,j,1:2),nwpo%psfc(i,j,1:2),nwpo%zsfc(i,j,1:2) 
  
   IF ((nwpo%psfc(i,j,1) < 0.0) .or. (nwpo%psfc(i,j,2) < 0.0)) THEN
      nwp41%psfc(i,j) = missing_value_real4
   ELSE
      nwp41%psfc(i,j) = w1 * nwpo%psfc(i,j,1)   + w2 * nwpo%psfc(i,j,2)
   ENDIF
 
   IF ((nwpo%pmsl(i,j,1) < 0.0) .or. (nwpo%pmsl(i,j,2) < 0.0)) THEN
      nwp41%pmsl(i,j) = missing_value_real4
   ELSE
      nwp41%pmsl(i,j) = w1 * nwpo%pmsl(i,j,1)   + w2 * nwpo%pmsl(i,j,2)
   ENDIF

   IF ((nwpo%tsfc(i,j,1) < 0.0) .or. (nwpo%tsfc(i,j,2) < 0.0)) THEN
      nwp41%tsfc(i,j) = missing_value_real4
   ELSE
      nwp41%tsfc(i,j) = w1 * nwpo%tsfc(i,j,1)   + w2 * nwpo%tsfc(i,j,2)
   ENDIF

   IF ((nwpo%zsfc(i,j,1) < 0.0) .or. (nwpo%zsfc(i,j,2) < 0.0)) THEN
      nwp41%zsfc(i,j) = missing_value_real4
   ELSE
      nwp41%zsfc(i,j) = w1 * nwpo%zsfc(i,j,1)   + w2 * nwpo%zsfc(i,j,2)
   ENDIF
 
   IF ((nwpo%albedo(i,j,1) < 0.0) .or. (nwpo%albedo(i,j,2) < 0.0)) THEN
      nwp41%albedo(i,j) = missing_value_real4
   ELSE
      nwp41%albedo(i,j) = w1 * nwpo%albedo(i,j,1)   + w2 * nwpo%albedo(i,j,2)
   ENDIF
   
   IF ((nwpo%t_sigma(i,j,1) < 0.0) .or. (nwpo%t_sigma(i,j,2) < 0.0)) THEN
      nwp41%t_sigma(i,j) = missing_value_real4
   ELSE
      nwp41%t_sigma(i,j) = w1 * nwpo%t_sigma(i,j,1)   + w2 * nwpo%t_sigma(i,j,2)
   ENDIF

   IF ((nwpo%rh_sigma(i,j,1) < 0.0) .or. (nwpo%rh_sigma(i,j,2) < 0.0)) THEN
      nwp41%rh_sigma(i,j) = missing_value_real4
   ELSE
      nwp41%rh_sigma(i,j) = w1 * nwpo%rh_sigma(i,j,1)   + w2 * nwpo%rh_sigma(i,j,2)
   ENDIF

   IF ((nwpo%u_sigma(i,j,1) < 0.0) .or. (nwpo%u_sigma(i,j,2) < 0.0)) THEN
      nwp41%u_sigma(i,j) = missing_value_real4
   ELSE
      nwp41%u_sigma(i,j) = w1 * nwpo%u_sigma(i,j,1)   + w2 * nwpo%u_sigma(i,j,2)
   ENDIF

   IF ((nwpo%v_sigma(i,j,1) < 0.0) .or. (nwpo%v_sigma(i,j,2) < 0.0)) THEN
      nwp41%v_sigma(i,j) = missing_value_real4
   ELSE
      nwp41%v_sigma(i,j) = w1 * nwpo%v_sigma(i,j,1)   + w2 * nwpo%v_sigma(i,j,2)
   ENDIF   

   IF ((nwpo%tpw(i,j,1) < 0.0) .or. (nwpo%tpw(i,j,2) < 0.0)) THEN
      nwp41%tpw(i,j) = missing_value_real4
   ELSE
      nwp41%tpw(i,j) = w1 * nwpo%tpw(i,j,1)   + w2 * nwpo%tpw(i,j,2)
   ENDIF 

   IF ((nwpo%weasd(i,j,1) < 0.0) .or. (nwpo%weasd(i,j,2) < 0.0)) THEN
      nwp41%weasd(i,j) = missing_value_real4
   ELSE
      nwp41%weasd(i,j) = w1 * nwpo%weasd(i,j,1)   + w2 * nwpo%weasd(i,j,2)
   ENDIF    

   IF ((nwpo%o3col(i,j,1) < 0.0) .or. (nwpo%o3col(i,j,2) < 0.0)) THEN
      nwp41%o3col(i,j) = missing_value_real4
   ELSE
      nwp41%o3col(i,j) = w1 * nwpo%o3col(i,j,1)   + w2 * nwpo%o3col(i,j,2)
   ENDIF    
    
   IF ((nwpo%ttropo(i,j,1) < 0.0) .or. (nwpo%ttropo(i,j,2) < 0.0)) THEN
      nwp41%ttropo(i,j) = missing_value_real4
   ELSE
      nwp41%ttropo(i,j) = w1 * nwpo%ttropo(i,j,1)   + w2 * nwpo%ttropo(i,j,2)
   ENDIF  
   !print*,'ttrop',nwpo%ttropo(i,j,1) ,nwpo%ttropo(i,j,2),nwpo%o3col(i,j,1:2),nwpo%tpw(i,j,1:2) 
 
 
   !print*,nwpo%tlev(i,j,1:26,1)
   DO k=1,nwp41%nlevels


       IF ((nwpo%tlev(i,j,k,1) < 0.0) .or. (nwpo%tlev(i,j,k,2) < 0.0)) THEN
          nwp41%tlev(i,j,k) = missing_value_real4
       ELSE
          nwp41%tlev(i,j,k) = w1 * nwpo%tlev(i,j,k,1)   + w2 * nwpo%tlev(i,j,k,2) 
       ENDIF  

       IF ((nwpo%zlev(i,j,k,1) < 0.0) .or. (nwpo%zlev(i,j,k,2) < 0.0)) THEN
          nwp41%zlev(i,j,k) = missing_value_real4
       ELSE
          nwp41%zlev(i,j,k) = w1 * nwpo%zlev(i,j,k,1)   + w2 * nwpo%zlev(i,j,k,2)
       ENDIF 

       IF ((nwpo%o3lev(i,j,k,1) < 0.0) .or. (nwpo%o3lev(i,j,k,2) < 0.0)) THEN
          nwp41%o3lev(i,j,k) = missing_value_real4
       ELSE
          nwp41%o3lev(i,j,k) = w1 * nwpo%o3lev(i,j,k,1)   + w2 * nwpo%o3lev(i,j,k,2)
       ENDIF               

       IF ((nwpo%rhlev(i,j,k,1) < 0.0) .or. (nwpo%rhlev(i,j,k,2) < 0.0)) THEN
          nwp41%rhlev(i,j,k) = missing_value_real4
       ELSE
          nwp41%rhlev(i,j,k) = w1 * nwpo%rhlev(i,j,k,1)   + w2 * nwpo%rhlev(i,j,k,2)
       ENDIF  
       
       IF ((nwpo%clwlev(i,j,k,1) < 0.0) .or. (nwpo%clwlev(i,j,k,2) < 0.0)) THEN
          nwp41%clwlev(i,j,k) = missing_value_real4
       ELSE
          nwp41%clwlev(i,j,k) = w1 * nwpo%clwlev(i,j,k,1)   + w2 * nwpo%clwlev(i,j,k,2)
       ENDIF    
       

       ! convert rh to water vapor mix ration [g/kg]   
       IF (nwp41%rhlev(i,j,k) > 0.0) THEN
            
            !mpav (3-28-2008) - For the GFS RH is with respect to water/ice
            ! depending on temperature.
            CALL rh_to_wv(nwp41%rhlev(i,j,k), nwp41%plev(i,j,k), &
                          nwp41%tlev(i,j,k), buf_r4)
                            
            nwp41%wlev(i,j,k) = buf_r4
       
       ENDIF
          !print*,'t2 = ',i,j,k,nwpo%tlev(i,j,k,1:2),nwpo%rhlev(i,j,k,1:2),nwpo%zlev(i,j,k,1:2)
 
       
      IF ((nwpo%ulev(i,j,k,1) < -200.0) .or. (nwpo%ulev(i,j,k,2) < -200.0)) THEN
         nwp41%ulev(i,j,k) = missing_value_real4
      ELSE
         nwp41%ulev(i,j,k) = w1 * nwpo%ulev(i,j,k,1)  + w2 * nwpo%ulev(i,j,k,2)
      ENDIF  
       
      IF ((nwpo%vlev(i,j,k,1) < -200.0) .or. (nwpo%vlev(i,j,k,2) < -200.0)) THEN
         nwp41%vlev(i,j,k) = missing_value_real4
      ELSE
         nwp41%vlev(i,j,k) = w1 * nwpo%vlev(i,j,k,1)  + w2 * nwpo%vlev(i,j,k,2)
      ENDIF 
         
   END DO ! k=1, nwp26%nlevels
       
   !------------------------------------------------------------------
   ! Extrapolate NWP profile from top level of RH profile to last 
   ! level of the pressure/temperature profile.
   !------------------------------------------------------------------
        
   wo = max(nwp41%wlev(i,j,kfirst_rh),0.0003)
   DO k = kfirst_rh-1, 1, -1
      nwp41%wlev(i,j,k) = max((wo*(nwp41%plev(i,j,k)/nwp41%plev(i,j,kfirst_rh))**3),0.0003)
   END DO
   
   DO k=1, nwp41%nlevels
      IF (nwp41%wlev(i,j,k) <= 0.0) THEN
          nwp41%wlev(i,j,k) = 0.0003
      ENDIF
   END DO
    
   
END DO
END DO

ENDIF

!=== 2.5. deallocate nwpo arrays
call deallocate_nwpo_arrays

! 3. end subroutine   
end subroutine fylat_time_interp_nwp
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~ subroutine 2: cal_julian ~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine cal_julian

!-----------------------------------------------------------------------
! !F90 cal_julian
!
! !Description:
!    This program is to calculate julian times of satellite and nwp data.
!
! !Input  parameters:
!    none
!
! !Output parameters:
!    none
!
!-----------------------------------------------------------------------

implicit none

! 1. define variables
integer(kind=4) :: y1,m1,d1,h1
integer(kind=4) :: length, ierr
character(len=100) :: ins

!*******
! 2. begin program
! 2.1. satellite time
call julian (sat%year, sat%month, sat%day, sat%hour, sat%mint, jutime%sate)

! 2.2. convert char to INTEGER for nwp time and get nwp julian time
ins=nwptime%year(1)
call ICNVRT(1,y1,ins,length,ierr)
ins=nwptime%month(1)
call ICNVRT(1,m1,ins,length,ierr)
ins=nwptime%day(1)
call ICNVRT(1,d1,ins,length,ierr)
ins=nwptime%hour(1)
call ICNVRT(1,h1,ins,length,ierr)
call julian (y1, m1, d1, h1, 0, jutime%nwp1)

ins=nwptime%year(2)
call ICNVRT(1,y1,ins,length,ierr)
ins=nwptime%month(2)
call ICNVRT(1,m1,ins,length,ierr)
ins=nwptime%day(2)
call ICNVRT(1,d1,ins,length,ierr)
ins=nwptime%hour(2)
call ICNVRT(1,h1,ins,length,ierr)
call julian (y1, m1, d1, h1, 0, jutime%nwp2)


! 3. end subroutine   
end subroutine cal_julian
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~ subroutine 3: cal_weight ~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine cal_weight(w1,w2)

!-----------------------------------------------------------------------
! !F90 cal_weight
!
! !Description:
!    This program is to calculate weights of julian times.
!
! !Input  parameters:
!    none
!
! !Output parameters:
!    none
!
!-----------------------------------------------------------------------

implicit none

! 1. define variables
real(kind=8) :: t, t1, t2
real(kind=4) :: w1, w2

!*******
! 2. begin program
t  = jutime%sate
t1 = jutime%nwp1
t2 = jutime%nwp2

w1 = 1.D0 - (t-t1)/(t2-t1)
w2 = (t-t1)/(t2-t1)

! 3. end subroutine   
end subroutine cal_weight
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~ subroutine 4: interp101_NWP data ~~~~~~~~~~~~~~~~~~
subroutine fylat_interp101_nwp()

!-----------------------------------------------------------------------
! !F90 interp101_NWP
!
! !Description:
!    This program is to interpolate NWP arrays to 101 layers
!
! !Input  parameters:
!    none
!
! !Output parameters:
!    none
!
!-----------------------------------------------------------------------

implicit none
! 1. define variables
!real(kind=4), dimension(101) :: o3
real(kind=4),   dimension(:), pointer :: pp,tt,zz,oo,ww,tp,hh,uu,vv
integer(kind=4), dimension(:), pointer :: iv
integer(kind=4) :: ix, iy, i, j, k, xn, yn
integer(kind=1) :: status

!*******
! 2. begin program
  print*,'  ... interpolate nwp data to 101 layers'


!=== previous deal
call allocate_nwp_arrays1

if (fylat_nwp_opt == 1 .or. fylat_nwp_opt == 2 .or. fylat_nwp_opt == 4 .or. fylat_nwp_opt == 5) then   ! 1=ncep and 2=gfs

!=== 2.1. interp start
nwp%nlevels = NLEVELS_INTERP
nwp%dat%psfc = nwp26%psfc/100.
nwp%dat%pmsl = nwp26%pmsl/100.
nwp%dat%tsfc = nwp26%tsfc
nwp%dat%zsfc = nwp26%zsfc
nwp%dat%albedo = nwp26%albedo
nwp%dat%t_sigma = nwp26%t_sigma
nwp%dat%rh_sigma = nwp26%rh_sigma
nwp%dat%u_sigma = nwp26%u_sigma
nwp%dat%v_sigma = nwp26%v_sigma
nwp%dat%tpw     = nwp26%tpw
nwp%dat%weasd   = nwp26%weasd
nwp%dat%o3col   = nwp26%o3col
nwp%dat%ttropo  = nwp26%ttropo


!=== 2.2. rh and o3 calculation
!call rh_cal(rh)
!call o3_cal(o3)

xn = nwp%nlon
yn = nwp%nlat
if (fylat_nwp_opt == 5) then   
   xn = nwp%nlon05
   yn = nwp%nlat05 
endif


do ix = 1, xn
do iy = 1, yn

  nwp%dat(ix,iy)%lon     = nwp26%lon(ix,iy)
  nwp%dat(ix,iy)%lat     = nwp26%lat(ix,iy)
  
  call allocate_nwp_arrays2(ix,iy)
  
!=== 2.3. assign pointers
  pp => nwp%dat(ix,iy)%plev
  zz => nwp%dat(ix,iy)%zlev
  tt => nwp%dat(ix,iy)%tlev
  ww => nwp%dat(ix,iy)%wlev
  oo => nwp%dat(ix,iy)%o3lev
  tp => nwp%dat(ix,iy)%tpwlev
  iv => nwp%dat(ix,iy)%inversion_lev
  hh => nwp%dat(ix,iy)%rhlev
  uu => nwp%dat(ix,iy)%ulev
  vv => nwp%dat(ix,iy)%vlev
    
!=== 2.4. interp
  
  !------------------------------------------------------------------
  ! Create the higher resolution vertical pressure levels. 
  !------------------------------------------------------------------
  call make_profile_101(pp)

  !------------------------------------------------------------------
  ! Interpolate the current temperature and water vapor profiles to
  ! the higher resolution pressure levels. 
  !------------------------------------------------------------------
  call profile_to_101(nwp26%plev(ix,iy,:), &
                      nwp26%tlev(ix,iy,:), &
                      nwp26%wlev(ix,iy,:), &
                      nwp26%nlevels,       &
                      nwp%dat(ix,iy)%lat,  &
                      pp,                  &  !nwp%dat(ix,iy)%plev, &
                      tt,                  &  !nwp%dat(ix,iy)%tlev, &
                      ww,                  &  !nwp%dat(ix,iy)%wlev, &
                      status)         
  if (status /= 0) then
     print*, "(a,'Error interpolating profile in the vertical - aborting')"
     stop
  endif  

  !call deal101(i, j, pp, zz, tt, ww)
  
  !------------------------------------------------------------------
  ! read o3 profile from climatological results
  !------------------------------------------------------------------
  call climoz_101_o3(nwp%dat(ix,iy)%lat, sat%month, oo)
             
  !------------------------------------------------------------------
  ! Calculate a high resolution high profile. 
  !------------------------------------------------------------------
  call height_profile(pp, &
                      tt, &
                      ww, &
                      zz, &
                      nwp%nlevels,         &
                      nwp%dat(ix,iy)%pmsl)
                      
  !------------------------------------------------------------------
  ! Compute the RH profile, assuming a GFS liquid/ice weighting
  ! scheme - MPAV (3/21/2011
  !------------------------------------------------------------------
  
  call compute_rh_profile(pp, &
                          tt, &
                          ww, &
                          nwp%nlevels,         &
                          nwp%dat(ix,iy)%rhlev)

!--- obtain tpw profile/tropopause-stratopause_level of 101 layers

  call compute_tpw_profile(pp, ww, tp) 
   
  call compute_levels(ix, iy, iv)
  
  nwp%dat(ix,iy)%tsfc_uni = determine_nwp_quality_flag(nwp26%tsfc, &
                                                           ix,iy,  &
                                                              xn,  &
                                                              yn,  &
                                                            1, 1)
                  
  pp => null()
  tt => null()
  zz => null()
  ww => null()
  oo => null()
  tp => null()
  iv => null()
  uu => null()
  vv => null() 
     !print*,ix,iy, nwp%dat(ix,iy)%sfc_level,nwp%dat(ix,iy)%psfc

end do
end do
!stop
!=== 2.5. dellocate nwp26
!call deallocate_nwp26_arrays

endif  !if (fylat_nwp_opt == 1 .or. fylat_nwp_opt == 2) then   ! 1=ncep and 2=gfs


if (fylat_nwp_opt == 3) then   ! 3=T639

!=== 2.1. interp start
nwp%nlevels = NLEVELS_INTERP
nwp%dat%psfc = nwp36%psfc/100.
nwp%dat%pmsl = nwp36%pmsl/100.
nwp%dat%tsfc = nwp36%tsfc
nwp%dat%zsfc = nwp36%zsfc
nwp%dat%albedo = nwp36%albedo
nwp%dat%t_sigma = nwp36%t_sigma
nwp%dat%rh_sigma = nwp36%rh_sigma
nwp%dat%u_sigma = nwp36%u_sigma
nwp%dat%v_sigma = nwp36%v_sigma
nwp%dat%tpw     = nwp36%tpw
nwp%dat%weasd   = nwp36%weasd
nwp%dat%o3col   = nwp36%o3col
nwp%dat%ttropo  = nwp36%ttropo


!=== 2.2. rh and o3 calculation
!call rh_cal(rh)
!call o3_cal(o3)


do ix = 1, nwp%nlon_T639
do iy = 1, nwp%nlat_T639

  nwp%dat(ix,iy)%lon     = nwp36%lon(ix,iy)
  nwp%dat(ix,iy)%lat     = nwp36%lat(ix,iy)
  
  call allocate_nwp_arrays2(ix,iy)
  
!=== 2.3. assign pointers
  pp => nwp%dat(ix,iy)%plev
  zz => nwp%dat(ix,iy)%zlev
  tt => nwp%dat(ix,iy)%tlev
  ww => nwp%dat(ix,iy)%wlev
  oo => nwp%dat(ix,iy)%o3lev
  tp => nwp%dat(ix,iy)%tpwlev
  iv => nwp%dat(ix,iy)%inversion_lev
  hh => nwp%dat(ix,iy)%rhlev
  uu => nwp%dat(ix,iy)%ulev
  vv => nwp%dat(ix,iy)%vlev

!=== 2.4. interp
  
  !------------------------------------------------------------------
  ! Create the higher resolution vertical pressure levels. 
  !------------------------------------------------------------------
  call make_profile_101(pp)

  !------------------------------------------------------------------
  ! Interpolate the current temperature and water vapor profiles to
  ! the higher resolution pressure levels. 
  !------------------------------------------------------------------

  call profile_to_101_uv(nwp26%plev(ix,iy,:), &
                         nwp26%ulev(ix,iy,:), &
                         nwp26%vlev(ix,iy,:), &
                         nwp26%nlevels,       &
                         nwp%dat(ix,iy)%lat,  &
                         pp,                  &  !nwp%dat(ix,iy)%plev, &
                         uu,                  &  !nwp%dat(ix,iy)%tlev, &
                         vv,                  &  !nwp%dat(ix,iy)%wlev, &
                         status)                            
  if (status /= 0) then
     print*, "(a,'Error interpolating profile in the vertical - aborting')"
     stop
  endif  
  
  !call deal101(i, j, pp, zz, tt, ww)
  
  !------------------------------------------------------------------
  ! read o3 profile from climatological results
  !------------------------------------------------------------------
  call climoz_101_o3(nwp%dat(ix,iy)%lat, sat%month, oo)
             
  !------------------------------------------------------------------
  ! Calculate a high resolution high profile. 
  !------------------------------------------------------------------
  call height_profile(pp, &
                      tt, &
                      ww, &
                      zz, &
                      nwp%nlevels,         &
                      nwp%dat(ix,iy)%pmsl)

  !------------------------------------------------------------------
  ! Compute the RH profile, assuming a GFS liquid/ice weighting
  ! scheme - MPAV (3/21/2011
  !------------------------------------------------------------------
  
  call compute_rh_profile(pp, &
                          tt, &
                          ww, &
                          nwp%nlevels,         &
                          nwp%dat(ix,iy)%rhlev)

!--- obtain tpw profile/tropopause-stratopause_level of 101 layers

  call compute_tpw_profile(pp,  ww, tp) 
   
  call compute_levels(ix, iy, iv)
  
  nwp%dat(ix,iy)%tsfc_uni = determine_nwp_quality_flag(nwp36%tsfc, &
                                                           ix,iy,  &
                                                    nwp%nlon_T639, &
                                                    nwp%nlat_T639, &
                                                            1, 1)
                                   
  pp => null()
  tt => null()
  zz => null()
  ww => null()
  oo => null()
  tp => null()
  iv => null()
  uu => null()
  vv => null()   
 !print*,ix,iy, nwp%dat(ix,iy)%sfc_level,nwp%dat(ix,iy)%psfc

end do
end do
!stop
!=== 2.5. dellocate nwp36
!call deallocate_nwp36_arrays

endif  !if (fylat_nwp_opt == 3) then   ! 3=T639


IF (fylat_nwp_opt == 6) THEN   ! 6=grapes gfs

!=== 2.1. interp start
nwp%nlevels = NLEVELS_INTERP
nwp%dat%psfc = nwp40%psfc/100.
nwp%dat%pmsl = nwp40%pmsl/100.
nwp%dat%tsfc = nwp40%tsfc
nwp%dat%zsfc = nwp40%zsfc
nwp%dat%albedo = nwp40%albedo
nwp%dat%t_sigma = nwp40%t_sigma
nwp%dat%rh_sigma = nwp40%rh_sigma
nwp%dat%u_sigma = nwp40%u_sigma
nwp%dat%v_sigma = nwp40%v_sigma
nwp%dat%tpw     = nwp40%tpw
nwp%dat%weasd   = nwp40%weasd
nwp%dat%o3col   = nwp40%o3col
nwp%dat%ttropo  = nwp40%ttropo
!nwp%dat%u10m    = nwp40%u10m
!nwp%dat%v10m    = nwp40%v10m


!=== 2.2. rh and o3 calculation
!CALL rh_cal(rh)
!CALL o3_cal(o3)


DO ix = 1, nwp%nlon25
DO iy = 1, nwp%nlat25

  xn = nwp%nlon25
  yn = nwp%nlat25

  nwp%dat(ix,iy)%lon     = nwp40%lon(ix,iy)
  nwp%dat(ix,iy)%lat     = nwp40%lat(ix,iy)
  
  CALL allocate_nwp_arrays2(ix,iy)
  
!=== 2.3. assign pointers
  pp => nwp%dat(ix,iy)%plev
  zz => nwp%dat(ix,iy)%zlev
  tt => nwp%dat(ix,iy)%tlev
  ww => nwp%dat(ix,iy)%wlev
  oo => nwp%dat(ix,iy)%o3lev
  tp => nwp%dat(ix,iy)%tpwlev
  iv => nwp%dat(ix,iy)%inversion_lev
  hh => nwp%dat(ix,iy)%rhlev
  uu => nwp%dat(ix,iy)%ulev
  vv => nwp%dat(ix,iy)%vlev
    
!=== 2.4. interp
  
  !------------------------------------------------------------------
  ! Create the higher resolution vertical pressure levels. 
  !------------------------------------------------------------------
  CALL make_profile_101(pp)

  CALL profile_to_101_uv(nwp40%plev(ix,iy,11:40), &
                         nwp40%ulev(ix,iy,11:40), &
                         nwp40%vlev(ix,iy,11:40), &
                         nwp40%nlevels-10,        &
                         nwp%dat(ix,iy)%lat,  &
                         pp,                  &  !nwp%dat(ix,iy)%plev, &
                         uu,                  &  !nwp%dat(ix,iy)%tlev, &
                         vv,                  &  !nwp%dat(ix,iy)%wlev, &
                         status)         

  IF (status /= 0) THEN
     PRINT*, "(a,'Error interpolating profile in the vertical - aborting')"
     STOP
  ENDIF  
  !------------------------------------------------------------------
  ! Interpolate the current temperature and water vapor profiles to
  ! the higher resolution pressure levels. 
  !------------------------------------------------------------------
  CALL profile_to_101(nwp40%plev(ix,iy,13:40), &
                      nwp40%tlev(ix,iy,13:40), &
                      nwp40%wlev(ix,iy,13:40), &
                      nwp40%nlevels-12,        &
                      nwp%dat(ix,iy)%lat,  &
                      pp,                  &  !nwp%dat(ix,iy)%plev, &
                      tt,                  &  !nwp%dat(ix,iy)%tlev, &
                      ww,                  &  !nwp%dat(ix,iy)%wlev, &
                      status)
                 
  IF (status /= 0) THEN
     PRINT*, "(a,'Error interpolating profile in the vertical - aborting')"
     STOP
  ENDIF  

  !CALL deal101(i, j, pp, zz, tt, ww)
  
  !------------------------------------------------------------------
  ! read o3 profile from climatological results
  !------------------------------------------------------------------
  CALL climoz_101_o3(nwp%dat(ix,iy)%lat, sat%month, oo)
             
  !------------------------------------------------------------------
  ! Calculate a high resolution high profile. 
  !------------------------------------------------------------------
  CALL height_profile(pp, &
                      tt, &
                      ww, &
                      zz, &
                      nwp%nlevels,         &
                      nwp%dat(ix,iy)%pmsl)

  !------------------------------------------------------------------
  ! Compute the RH profile, assuming a GFS liquid/ice weighting
  ! scheme - MPAV (3/21/2011
  !------------------------------------------------------------------
  
  CALL compute_rh_profile(pp, &
                          tt, &
                          ww, &
                          nwp%nlevels,         &
                          nwp%dat(ix,iy)%rhlev)

!--- obtain tpw profile/tropopause-stratopause_level of 101 layers

  CALL compute_tpw_profile(pp,  ww, tp) 
   
  CALL compute_levels(ix, iy, iv)
  
  nwp%dat(ix,iy)%tsfc_uni = determine_nwp_quality_flag(nwp40%tsfc, &
                                                           ix,iy,  &
                                                       nwp%nlon25, &
                                                       nwp%nlat25, &
                                                            1, 1)
    !  print*,ix,iy,nwp40%wlev(ix,iy,1:40)
    !    print*,ix,iy,' grapes tt=', tt
    !    print*,ix,iy,' grapes ww=', ww
                        
  pp => null()
  tt => null()
  zz => null()
  ww => null()
  oo => null()
  tp => null()
  iv => null()
  uu => null()
  vv => null()  
 !print*,ix,iy, nwp%dat(ix,iy)%sfc_level,nwp%dat(ix,iy)%psfc
 
  ! add pwat
  !nwp%dat(ix,iy)%tpw = sum(nwp%dat(ix,iy)%tpwlev(:))
  nwp%dat(ix,iy)%tpw = nwp%dat(ix,iy)%tpwlev(NLEVELS_INTERP)
  !print*,ix,iy,'tpw',nwp%dat(ix,iy)%tpw
  
  ! add total ozone
  CALL O3_to_Dobson(nwp%dat(ix,iy)%plev,nwp%dat(ix,iy)%o3lev,nwp%dat(ix,iy)%o3col)
  !print*,ix,iy,'o3',nwp%dat(ix,iy)%o3col

END DO
END DO
!stop
!=== 2.5. dellocate nwp40
!CALL deallocate_nwp40_arrays

ENDIF  !IF (nwp_opt == 9) THEN   ! 9=grapes gfs


IF (fylat_nwp_opt == 7 .or. fylat_nwp_opt == 8) THEN   ! 8=grib2 gfs  0p25

!print*,nwp26%psfc(100,100),nwp26%pmsl(100,100)
!print*,nwp26%zsfc(100:120,100:120),nwp26%tsfc(100,100)
!=== 2.1. interp start
nwp%nlevels = NLEVELS_INTERP
nwp%dat%psfc = nwp31%psfc/100.
nwp%dat%pmsl = nwp31%pmsl/100.
nwp%dat%tsfc = nwp31%tsfc
nwp%dat%zsfc = nwp31%zsfc
nwp%dat%albedo = nwp31%albedo
nwp%dat%t_sigma = nwp31%t_sigma
nwp%dat%rh_sigma = nwp31%rh_sigma
nwp%dat%u_sigma = nwp31%u_sigma
nwp%dat%v_sigma = nwp31%v_sigma
nwp%dat%tpw     = nwp31%tpw
nwp%dat%weasd   = nwp31%weasd
nwp%dat%o3col   = nwp31%o3col
nwp%dat%ttropo  = nwp31%ttropo


!=== 2.2. rh and o3 calculation
!CALL rh_cal(rh)
!CALL o3_cal(o3)


IF (fylat_nwp_opt == 7 .or. fylat_nwp_opt == 8) THEN
   xn = nwp%nlon25_gfs
   yn = nwp%nlat25_gfs
ENDIF


DO ix = 1, xn !nwp%nlon
DO iy = 1, yn !nwp%nlat
  !print*,'ix ',ix, iy , nwp26%lon(ix,iy),nwp26%lat(ix,iy),nwp%dat(ix,iy)%plev
  nwp%dat(ix,iy)%lon     = nwp31%lon(ix,iy)
  nwp%dat(ix,iy)%lat     = nwp31%lat(ix,iy)

  CALL allocate_nwp_arrays2(ix,iy)
  
!=== 2.3. assign pointers
  pp => nwp%dat(ix,iy)%plev
  zz => nwp%dat(ix,iy)%zlev
  tt => nwp%dat(ix,iy)%tlev
  ww => nwp%dat(ix,iy)%wlev
  oo => nwp%dat(ix,iy)%o3lev
  tp => nwp%dat(ix,iy)%tpwlev
  iv => nwp%dat(ix,iy)%inversion_lev
  hh => nwp%dat(ix,iy)%rhlev
  uu => nwp%dat(ix,iy)%ulev
  vv => nwp%dat(ix,iy)%vlev
  
!=== 2.4. interp
  
  !------------------------------------------------------------------
  ! Create the higher resolution vertical pressure levels. 
  !------------------------------------------------------------------
  CALL make_profile_101(pp)

  CALL profile_to_101_uv(nwp31%plev(ix,iy,:), &
                         nwp31%ulev(ix,iy,:), &
                         nwp31%vlev(ix,iy,:), &
                         nwp31%nlevels,       &
                         nwp%dat(ix,iy)%lat,  &
                         pp,                  &  !nwp%dat(ix,iy)%plev, &
                         uu,                  &  !nwp%dat(ix,iy)%tlev, &
                         vv,                  &  !nwp%dat(ix,iy)%wlev, &
                         status)         
  IF (status /= 0) THEN
     PRINT*, "(a,'Error interpolating profile in the vertical - aborting')"
     STOP
  ENDIF    
  !print*,'uv', uu(80:90),vv(80:90)
  !------------------------------------------------------------------
  ! Interpolate the current temperature and water vapor profiles to
  ! the higher resolution pressure levels. 
  !------------------------------------------------------------------
  CALL profile_to_101(nwp31%plev(ix,iy,:), &
                      nwp31%tlev(ix,iy,:), &
                      nwp31%wlev(ix,iy,:), &
                      nwp31%nlevels,       &
                      nwp%dat(ix,iy)%lat,  &
                      pp,                  &  !nwp%dat(ix,iy)%plev, &
                      tt,                  &  !nwp%dat(ix,iy)%tlev, &
                      ww,                  &  !nwp%dat(ix,iy)%wlev, &
                      status)         
  IF (status /= 0) THEN
     PRINT*, "(a,'Error interpolating profile in the vertical - aborting')"
     STOP
  ENDIF  

  !CALL deal101(i, j, pp, zz, tt, ww)
  
  !------------------------------------------------------------------
  ! read o3 profile from climatological results
  !------------------------------------------------------------------
  CALL climoz_101_o3(nwp%dat(ix,iy)%lat, sat%month, oo)
             
  !------------------------------------------------------------------
  ! Calculate a high resolution high profile. 
  !------------------------------------------------------------------
  CALL height_profile(pp, &
                      tt, &
                      ww, &
                      zz, &
                      nwp%nlevels,         &
                      nwp%dat(ix,iy)%pmsl)
                      
  !------------------------------------------------------------------
  ! Compute the RH profile, assuming a GFS liquid/ice weighting
  ! scheme - MPAV (3/21/2011
  !------------------------------------------------------------------
  
  CALL compute_rh_profile(pp, &
                          tt, &
                          ww, &
                          nwp%nlevels,         &
                          nwp%dat(ix,iy)%rhlev)

!--- obtain tpw profile/tropopause-stratopause_level of 101 layers

  CALL compute_tpw_profile(pp,  ww, tp) 
   
  CALL compute_levels(ix, iy, iv)
  
  nwp%dat(ix,iy)%tsfc_uni = determine_nwp_quality_flag(nwp31%tsfc, &
                                                           ix,iy,  &
                                                  nwp%nlon25_gfs,  &
                                                  nwp%nlat25_gfs,  &
                                                            1, 1)
                               
  pp => null()
  tt => null()
  zz => null()
  ww => null()
  oo => null()
  tp => null()
  iv => null()
  uu => null()
  vv => null()
  !print*,ix,iy, nwp%dat(ix,iy)%sfc_level,nwp%dat(ix,iy)%psfc !,nwp%dat(ix,iy)%tlev
  !print*,ix,iy, nwp%dat(ix,iy)%plev,'fgg',nwp%dat(ix,iy)%tlev,'fee',nwp%dat(ix,iy)%zlev,'fzz',nwp%dat(ix,iy)%wlev,'fff',nwp%dat(ix,iy)%o3lev
  !stop
END DO
END DO
!stop

END IF !(gfxxgat_nwp_opt == 7 .or. gfxxgat_nwp_opt == 8)


IF (fylat_nwp_opt == 9 .or. fylat_nwp_opt == 10) THEN   ! 8=grib2 gfs  0p25

!print*,nwp41%psfc(100,100),nwp41%pmsl(100,100)
!print*,nwp26%zsfc(100:120,100:120),nwp26%tsfc(100,100)
!=== 2.1. interp start
nwp%nlevels = NLEVELS_INTERP
nwp%dat%psfc = nwp41%psfc/100.
nwp%dat%pmsl = nwp41%pmsl/100.
nwp%dat%tsfc = nwp41%tsfc
nwp%dat%zsfc = nwp41%zsfc
nwp%dat%albedo = nwp41%albedo
nwp%dat%t_sigma = nwp41%t_sigma
nwp%dat%rh_sigma = nwp41%rh_sigma
nwp%dat%u_sigma = nwp41%u_sigma
nwp%dat%v_sigma = nwp41%v_sigma
nwp%dat%tpw     = nwp41%tpw
nwp%dat%weasd   = nwp41%weasd
nwp%dat%o3col   = nwp41%o3col
nwp%dat%ttropo  = nwp41%ttropo


!=== 2.2. rh and o3 calculation
!CALL rh_cal(rh)
!CALL o3_cal(o3)


IF (fylat_nwp_opt == 9) THEN
   xn = nwp%nlon05
   yn = nwp%nlat05
ENDIF

IF (fylat_nwp_opt == 10) THEN
   xn = nwp%nlon25_gfs
   yn = nwp%nlat25_gfs
ENDIF


DO ix = 1, xn !nwp%nlon
DO iy = 1, yn !nwp%nlat
  !print*,'ix ',ix, iy , nwp41%lon(ix,iy),nwp41%lat(ix,iy),nwp%dat(ix,iy)%plev
  nwp%dat(ix,iy)%lon     = nwp41%lon(ix,iy)
  nwp%dat(ix,iy)%lat     = nwp41%lat(ix,iy)

  CALL allocate_nwp_arrays2(ix,iy)
  
!=== 2.3. assign pointers
  pp => nwp%dat(ix,iy)%plev
  zz => nwp%dat(ix,iy)%zlev
  tt => nwp%dat(ix,iy)%tlev
  ww => nwp%dat(ix,iy)%wlev
  oo => nwp%dat(ix,iy)%o3lev
  tp => nwp%dat(ix,iy)%tpwlev
  iv => nwp%dat(ix,iy)%inversion_lev
  hh => nwp%dat(ix,iy)%rhlev
  uu => nwp%dat(ix,iy)%ulev
  vv => nwp%dat(ix,iy)%vlev
  
!=== 2.4. interp
  
  !------------------------------------------------------------------
  ! Create the higher resolution vertical pressure levels. 
  !------------------------------------------------------------------
  CALL make_profile_101(pp)

  CALL profile_to_101_uv(nwp41%plev(ix,iy,:), &
                         nwp41%ulev(ix,iy,:), &
                         nwp41%vlev(ix,iy,:), &
                         nwp41%nlevels,       &
                         nwp%dat(ix,iy)%lat,  &
                         pp,                  &  !nwp%dat(ix,iy)%plev, &
                         uu,                  &  !nwp%dat(ix,iy)%tlev, &
                         vv,                  &  !nwp%dat(ix,iy)%wlev, &
                         status)         
  IF (status /= 0) THEN
     PRINT*, "(a,'Error interpolating profile in the vertical - aborting')"
     STOP
  ENDIF    
  !print*,'uv', uu(80:90),vv(80:90)
  !------------------------------------------------------------------
  ! Interpolate the current temperature and water vapor profiles to
  ! the higher resolution pressure levels. 
  !------------------------------------------------------------------
  CALL profile_to_101(nwp41%plev(ix,iy,:), &
                      nwp41%tlev(ix,iy,:), &
                      nwp41%wlev(ix,iy,:), &
                      nwp41%nlevels,       &
                      nwp%dat(ix,iy)%lat,  &
                      pp,                  &  !nwp%dat(ix,iy)%plev, &
                      tt,                  &  !nwp%dat(ix,iy)%tlev, &
                      ww,                  &  !nwp%dat(ix,iy)%wlev, &
                      status)         
  IF (status /= 0) THEN
     PRINT*, "(a,'Error interpolating profile in the vertical - aborting')"
     STOP
  ENDIF  

  !CALL deal101(i, j, pp, zz, tt, ww)
  
  !------------------------------------------------------------------
  ! read o3 profile from climatological results
  !------------------------------------------------------------------
  CALL climoz_101_o3(nwp%dat(ix,iy)%lat, sat%month, oo)
             
  !------------------------------------------------------------------
  ! Calculate a high resolution high profile. 
  !------------------------------------------------------------------
  CALL height_profile(pp, &
                      tt, &
                      ww, &
                      zz, &
                      nwp%nlevels,         &
                      nwp%dat(ix,iy)%pmsl)
                      
  !------------------------------------------------------------------
  ! Compute the RH profile, assuming a GFS liquid/ice weighting
  ! scheme - MPAV (3/21/2011
  !------------------------------------------------------------------
  
  CALL compute_rh_profile(pp, &
                          tt, &
                          ww, &
                          nwp%nlevels,         &
                          nwp%dat(ix,iy)%rhlev)

!--- obtain tpw profile/tropopause-stratopause_level of 101 layers

  CALL compute_tpw_profile(pp,  ww, tp) 
   
  CALL compute_levels(ix, iy, iv)
  
  nwp%dat(ix,iy)%tsfc_uni = determine_nwp_quality_flag(nwp41%tsfc, &
                                                           ix,iy,  &
                                                              xn,  &
                                                              yn,  &
                                                            1, 1)
                               
  pp => null()
  tt => null()
  zz => null()
  ww => null()
  oo => null()
  tp => null()
  iv => null()
  uu => null()
  vv => null()
  !print*,ix,iy, nwp%dat(ix,iy)%sfc_level,nwp%dat(ix,iy)%psfc !,nwp%dat(ix,iy)%tlev
  !print*,ix,iy, nwp%dat(ix,iy)%plev,'fgg',nwp%dat(ix,iy)%tlev,'fee',nwp%dat(ix,iy)%zlev,'fzz',nwp%dat(ix,iy)%wlev,'fff',nwp%dat(ix,iy)%o3lev
  !stop
END DO
END DO

!stop
!=== 2.5. dellocate nwp26
!CALL deallocate_nwp31_arrays

ENDIF  !IF (gfxxgat_nwp_opt == 9 .or. gfxxgat_nwp_opt == 10) THEN   ! 1=ncep and 2=gfs


! 3. end subroutine   
end subroutine fylat_interp101_nwp


!::::::::: sub subroutines for interp101_NWP::::::::::::::
!%%%%%%%%%%%%%%%%%%  ......
! sub function 1  make 101 profile
subroutine make_profile_101(pp)
! 1. define variables    
   integer(kind=4) ::  i, j
   real(kind=4)            :: l
   real(kind=4), parameter :: a = -1.550789414500298e-04
   real(kind=4), parameter :: b = -5.593654380586063e-02
   real(kind=4), parameter :: c =  7.451622227151780e+00
   real(kind=4), intent(out),dimension(:) :: pp
   
   l = 101.
   
   do i = 1, 101
      pp(i) = pow((a*l*l + b*l + c), (7./2.))
      l = l - 1.
   enddo
   
end subroutine make_profile_101

!%%%%%%%%%%%%%%%%%%  ......
! sub function 2  profile_to_101
subroutine profile_to_101(p, t, w, n, lat, &  ! in
                          pp,tt,ww,status)    ! out
! 1. define variables    
   integer(kind=4) :: i_str,i_end, iy, iz, i, j
   integer(kind=4), parameter :: nlx=101
   real(kind=4), parameter :: wmin = 0.0003
   real(kind=4) :: wmax, rlogp
   real(kind=4),dimension(2) :: pb, tb, wb
   real(kind=4) :: anum,aden,delt,delw
                         
   real(kind=4),dimension(n),intent(in) :: p, t, w
   integer(kind=4), intent(in) :: n
   real(kind=4), intent(in)    :: lat  
   real(kind=4), dimension(:),intent(in)  :: pp 
   real(kind=4), dimension(:),intent(out) :: tt,ww                      
   integer(kind=1),intent(out) :: status
   

!*******
! 2. begin program   

   status = 0

   if (n < 2) then
      print*,"number of levels must be at least 2"
      status = -1
   endif

   !*make_profile_101(pp)
   call int_levels_pp(p, t, w, n, pp, tt, ww, nlx, i_str, i_end)

   if (i_str >= 35) then
      print*,"ERROR: temperature profile doesn't go high enough"
      status = -1
   endif

   do i=1, i_str  ! assign values for upper layers
      tt(i) = -1.
      ww(i) = wmin
   enddo
   
   call extem101(tt, lat)
     
   iz = n ! 26

   if (p(iz) < pp(nlx)) then
   
      pb(1) = p(iz)
      tb(1) = t(iz)
      wb(1) = w(iz)

      pb(2) = pp(nlx)

      iy = n - 1

      anum = log(pb(2) / pb(1))
      aden = log(pb(1) / p(iy))

      rlogp = anum / aden

      delt = t(iz) - t(iy)
      tb(2) = t(iz) + delt * rlogp

      delw = w(iz) - w(iy)
      wb(2) = w(iz) + delw * rlogp

      if (wb(2) < wmin) then
          wb(2) = wmin
      else 
          wmax = satmix(pb(2), tb(2))

          if (wb(2) > wmax) then
              wb(2) = wmax
          endif
      endif

      i = i_end + 1

      call int_levels_pp(pb, tb, wb, 2, pp(i:nlx), tt(i:nlx), ww(i:nlx), nlx, i_str, i_end)
      
   endif  !if (p(iz) < pp(nlx)) then
   
end subroutine profile_to_101


subroutine profile_to_101_uv(p, t, w, n, lat, &  ! in
                             pp,tt,ww,status)    ! out
! 1. define variables    
   INTEGER(KIND=4) :: i_str,i_end, iy, iz,i
   INTEGER(KIND=4), PARAMETER :: nlx=101
   REAL(KIND=4), PARAMETER :: wmin = 0.0
   REAL(KIND=4) :: wmax
   REAL(KIND=4),DIMENSION(2) :: pb, tb, wb
   REAL(KIND=4) :: anum,aden,delt,delw
                         
   REAL(KIND=4),DIMENSION(n),INTENT(in) :: p, t, w
   INTEGER(KIND=4), INTENT(in) :: n
   REAL(KIND=4), INTENT(in)    :: lat  
   REAL(KIND=4), DIMENSION(:),INTENT(in)  :: pp 
   REAL(KIND=4), DIMENSION(:),INTENT(out) :: tt,ww                      
   INTEGER(KIND=1),INTENT(out) :: status
   

!*******
! 2. begin program   

   status = 0

   IF (n < 2) THEN
      PRINT*,"number of levels must be at least 2"
      status = -1
   ENDIF

   !*make_profile_101(pp)
   CALL int_levels_pp(p, t, w, n, pp, tt, ww, nlx, i_str, i_end)

   IF (i_str >= 35) THEN
      PRINT*,"ERROR: temperature profile doesn't go high enough"
      status = -1
   ENDIF

   DO i=1, i_str  ! assign values for upper layers
      tt(i) = wmin
      ww(i) = wmin
   ENDDO
   
  
end subroutine profile_to_101_uv

!%%%%%%%%%%%%%%%%%%  ......
! sub function 3  int_levels_pp
subroutine int_levels_pp(p, t, w, n, pp, tt, ww, nn, i_str, i_end) 
                         
   integer(kind=4) :: n, nn ! nn=101, n=26 or 36
   real(kind=4),dimension(n) :: p, t, w
   real(kind=4), dimension(:) :: pp, tt, ww  
   integer(kind=4) ::  i_str, i_end
   real(kind=4)    :: dl, dl1
   real(kind=4)    :: slope_t
   real(kind=4)    :: slope_w
   integer(kind=4) :: ii, i, j, k

   if (n < 2) then
      i_str =  1
      i_end =  0
   endif
   
   do i=1, nn
     if (pp(i)>=p(1)) EXIT ! p(1)=10 mb
     tt(i) = 0.
     ww(i) = 0.
   enddo
   
   i_str = i
   
   ii = 2

   dl      = log(p(ii) / p(ii-1))
   slope_t = (t(ii) - t(ii-1)) / dl  ! calculate diff = V20mb - V10mb
   slope_w = (w(ii) - w(ii-1)) / dl
   
   do j=1, nn+50
   
     if (pp(i) < p(ii)) then
     
       dl    = log(pp(i) / p(ii-1))
       tt(i) = t(ii-1) + slope_t * dl
       ww(i) = w(ii-1) + slope_w * dl

       i=i+1
       
     else
      
       ii = ii +1

       if (ii>n) then !> 1000mb
          if (p(ii-1) == pp(i)) then
              tt(i) = t(ii-1)
              ww(i) = w(ii-1)
          endif 
          EXIT
       endif
       
       dl     = log(p(ii) / p(ii-1))
       slope_t = (t(ii) - t(ii-1)) / dl
       slope_w = (w(ii) - w(ii-1)) / dl
       
     endif
     
   enddo
   
   i_end = i - 1     ! last layer in 26 layers nwp
   
   if (n>10) then 
   do i = i_end+1, nn  ! assign values for upper layers 
      tt(i) = 0.
      ww(i) = 0.  
   enddo
   endif
   
end subroutine int_levels_pp

!%%%%%%%%%%%%%%%%%%  ......
! sub function 4  generate height profile
subroutine height_profile(p, t, w, z, n, p0) 

  integer(kind=4), intent(in) :: n ! n=26、36
  real(kind=4),intent(in)     :: p0
  real(kind=4),dimension(:),intent(in)  :: p, t, w
  real(kind=4),dimension(:),intent(out) :: z
  integer(kind=4) :: i_up, i_dn, i, j, k
  real(kind=4)   :: g, Rd, Rv, epsilon, &
                        t0, w0, z_cur, &
                        p_last, t_last, w_last, &
                        e, rho, PP, TT, WW
                  
  g = 9.80665  !* m / s^2	*/

  Rd = 287.05  !* J / kg K	*/
  Rv = 461.51  !* J / kg K	*/	
  epsilon = Rd / Rv

  do i=1,n
     if (p(i)>p0) then
        EXIT
     endif
  enddo

  if (i == 1) then
      i_up = i
      i_dn = i+1
  else 
      i_up = i-1
      i_dn = i
  endif

  t0 = t(i_up) + lin_int(p(i_up), p(i_dn), p0, t(i_up), t(i_dn))
  w0 = w(i_up) + lin_int(p(i_up), p(i_dn), p0, w(i_up), w(i_dn))

  z_cur = 0

  p_last = p0
  t_last = t0
  w_last = w0
  
  do i = i_up,1,-1
  
     PP = (p_last + p(i)) / 2. * 100.
     TT = (t_last + t(i)) / 2.
     WW = (w_last + w(i)) / 2. / 1000.
     
     e  = WW / (WW + epsilon) * PP
     rho= (PP-e) / (Rd*TT) + e / (Rv*TT)
     
     z_cur = z_cur + (p_last-p(i))*100. / (g*rho)

     z(i) = z_cur

     p_last = p(i)
     t_last = t(i)
     w_last = w(i)
     
  end do
  
  z_cur = 0

  p_last = p0
  t_last = t0
  w_last = w0
  
  do i = i_dn, n
     PP = (p_last + p(i)) / 2. * 100.
     TT = (t_last + t(i)) / 2.
     WW = (w_last + w(i)) / 2. / 1000.

     e = WW / (WW + epsilon) * PP

     rho = (PP-e) / (Rd*TT) + e / (Rv*TT)

     z_cur = z_cur + (p_last-p(i))*100. / (g*rho)

     z(i) = z_cur

     p_last = p(i)
     t_last = t(i)
     w_last = w(i)  
  end do
     
end subroutine height_profile

!---------------------------------------------------------------------------
! Convert mass mixing ratio to RH
!---------------------------------------------------------------------------

subroutine compute_rh_profile(plev,    &
                              tlev,    &
                              wlev,    &
                              nlevels, &
                              rhlev)

 real(kind=4), intent(in), dimension(:):: plev, tlev, wlev
 integer(kind=4), intent(in):: nlevels
 real(kind=4), intent(out), dimension(:):: rhlev
 integer:: i, j, k

 level_loop: do k=1, nlevels
   call mr_to_rh_gfs(rhlev(k), plev(k), tlev(k), wlev(k))
 end do level_loop

end subroutine compute_rh_profile

!%%%%%%%%%%%%%%%%%%  ......
subroutine mr_to_rh_gfs(rh, p, t, w) 

    real(kind=4), intent(in)  :: p,t,w   
    real(kind=4), intent(out) :: rh  
    real(kind=4)              :: weight
    real(kind=4)              :: wsat
     
    if (t >= 273.16) then
       weight = 1.0
    else if (t <= 253.16) then
       weight = 0.0
    else if (t < 273.16 .and. t > 253.16) then
       weight = (t - 253.16)/(273.16 - 253.16)
    else
       weight = 1.0
    endif
     
     wsat = (weight*satmixwat(p, t) + (1.0 - weight)*satmixice(p, t));
     rh = (w/wsat)*100.0
     
end subroutine mr_to_rh_gfs

!%%%%%%%%%%%%%%%%%%  ......
! sub function 5 for subroutine height_profile
function lin_int(x1,x2,x,y1,y2) RESULT(res)

  real(kind=4) :: x1,x2,x,y1,y2,res
  res = (x - x1) / (x2 - x1) * (y2 - y1)
  
end function lin_int
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~ subroutine 5: rh cal ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine rh_to_wv(rh,p,t,wv)

!-----------------------------------------------------------------------
! !F90 rh_to_wv
!
! !Description:
!    This program is to calculate water vaper mix ratio [g/kg] profile of 101 
!    layers using relative humidity profile of 101 layers.
!
!
! !Input  parameters:
!    t,p,rh          = temperature, pressure, relative humidity
!
! !Output parameters:
!    wv              = water vapor mix ratio [g/kg]
!
!-----------------------------------------------------------------------

implicit none
! 1. define variables
real(kind=4), intent(in) :: t, p, rh
real(kind=4), intent(out):: wv
!real(kind=4) :: E, p1
real(kind=4) :: weight

!*******
! 2. begin program

! --- old program:rh to water vapor mix ratio---------
!r=wstd/1000.  ! convert g/kg to kg/kg
!p=pstd*100.   ! convert hPa to Pa
!do i = 1, 101
!  E=10**((10.286*tstd(i)-2148.4909)/(tstd(i)-35.58))
!  rh(i)=100*r(i)*p(i)/(E*(0.62198+r(i)))
!enddo
!------------------------------------------------------
!----- old 2 ------------------------------------------
!E  = 10**((10.286*t-2148.4909)/(t-35.58))
!p1 = p*100. ! convert hPa to Pa
!wv = 1000*rh*E*0.62198 / (100*p1-rh*E) 
!------------------------------------------------------

weight = 1.0
if (t >= 273.16) weight = 1.0
if (t <= 253.16) weight = 0.0
if ((t < 273.16) .and. (t > 253.16)) weight = (t - 253.16)/(273.16 - 253.16)

wv = rh * 0.01 * ( weight*satmixwat(p, t) + (1.0 - weight)*satmixice(p, t) )

! 3. end subroutine   
end subroutine rh_to_wv

!%%%%%%%%%%%%%%%%%%  ......
! sub functions 1  for water
function satmixwat(p,t) RESULT(res)

  real(kind=4) :: p,t,es,res
  
  es  = svpwat(t)
  res = 622. * es / p
  
end function satmixwat

!%%%%%%%%%%%%%%%%%%  ......
! sub functions 2  for water
function svpwat(t1) RESULT(res)

  real(kind=4) :: t1,res,t
  real(kind=4), parameter :: a0 =  0.999996876e0
  real(kind=4), parameter :: a1 = -0.9082695004e-2
  real(kind=4), parameter :: a2 =  0.7873616869e-4
  real(kind=4), parameter :: a3 = -0.6111795727e-6
  real(kind=4), parameter :: a4 =  0.4388418740e-8
  real(kind=4), parameter :: a5 = -0.2988388486e-10
  real(kind=4), parameter :: a6 =  0.2187442495e-12
  real(kind=4), parameter :: a7 = -0.1789232111e-14
  real(kind=4), parameter :: a8 =  0.1111201803e-16
  real(kind=4), parameter :: a9 = -0.3099457145e-19
  real(kind=4), parameter :: b  =  0.61078e+1
  
  t = t1-273.16
  res = b / pow(a0+t*(a1+t*(a2+t*(a3+t*(a4+t*(a5+t*(a6+t*(a7+t*(a8+t*a9)))))))), 8.)

end function svpwat

!%%%%%%%%%%%%%%%%%%  ......
! sub functions 3  for ice
function satmixice(p,t) RESULT(res)
  real(kind=4) :: p,t,es,res
  
  es  = svpice(t)
  res = 622. * es / p
  
end function satmixice

function svpice(t) RESULT(res)
  real(kind=4) :: a,res, t
  
  a = 273.16 / t
  res = pow(10., -9.09718 * (a - 1.) - 3.56654 * log10(a) + &
                0.876793 * (1. - 1./a) + log10(6.1071))
                
end function svpice

!%%%%%%%%%%%%%%%%%%  ......
! sub functions 4  for ice
function satmix(p,t) RESULT(res)
  real(kind=4) :: p,t,es,res
  
  if (t > 253.) then
     res = satmixwat(p, t)
  else
     res = satmixice(p, t)
  endif
  
end function satmix
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~ subroutine 6: compute tpw profile ~~~~~~~~~~~~~~~ 
subroutine compute_tpw_profile(plev, wlev, tpwlev)  

!-----------------------------------------------------------------------
! !F90 compute_tpw_profile
!
! !Description:
!    This program is to compute tpw profile
!
! !Input  parameters:
!    plev            = pressure profile  
!    wlev            = water vapor mix ratio profile
!       
! !Output parameters:
!    tpwlev          = tpw profile
!
!-----------------------------------------------------------------------

 real(kind=4), intent(in), dimension(:) :: plev, wlev
! integer(kind=4), intent(in) :: nlevels
 real(kind=4), intent(out), dimension(:):: tpwlev
 integer(kind=4):: k
 real(kind=4)  :: w_mean, tpw_layer


!--- assume no water above highest level
  tpwlev(1) = 0.0

!--- construct profile at each level of integrated tpw
  level_loop: do k = 2, nwp%nlevels

    w_mean = 0.5*(wlev(k)+wlev(k-1)) / 1000.0         !kg/kg
    tpw_layer = (10.0/STANDARD_GRAVITY)*(plev(k)-plev(k-1))*w_mean   !cm
    tpwlev(k) = tpwlev(k-1) + tpw_layer               !cm! 3. end subroutine

  end do level_loop

! 3. end subroutine
end subroutine compute_tpw_profile
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~ subroutine 7: compute some levels ~~~~~~~~~~~~~~~ 
subroutine compute_levels(ix,iy,inverlev)
 
!-----------------------------------------------------------------------
! !F90 compute_levels
!
! !Description:
!    This program is to find some levels [sfc/tropopause/stratopause
!    /inversion]
!
! !Input  parameters:
!    ix,iy
!       
! !Output parameters:
!    none
!
!-----------------------------------------------------------------------

real(kind=4)    :: tmin, tmax, a 
integer(kind=4) :: ix, iy, i, j, k
integer(kind=4) :: il, isfc
integer(kind=4) :: ninver
integer(kind=4), dimension(:), intent(out) :: inverlev

  !------------------------------------------------------------------
  ! Find the lowest valid level [sfc]. 
  !------------------------------------------------------------------
    nwp%dat(ix,iy)%sfc_level = nwp%nlevels
    do k=nwp%nlevels/2, nwp%nlevels
       if (nwp%dat(ix,iy)%plev(k) > nwp%dat(ix,iy)%psfc) then
          nwp%dat(ix,iy)%sfc_level = k
          EXIT
       endif
    end do
    
    isfc = nwp%dat(ix,iy)%sfc_level
    il = isfc - 1
    a = (nwp%dat(ix,iy)%psfc - nwp%dat(ix,iy)%plev(il))/ &
        (nwp%dat(ix,iy)%plev(il+1)  - nwp%dat(ix,iy)%plev(il))
    !a = (nwp%dat(ix,iy)%zsfc - nwp%dat(ix,iy)%zlev(il))/ &
    !    (nwp%dat(ix,iy)%zlev(il+1) - nwp%dat(ix,iy)%zlev(il))
    nwp%dat(ix,iy)%plev(isfc) = nwp%dat(ix,iy)%psfc
    nwp%dat(ix,iy)%tlev(isfc) = nwp%dat(ix,iy)%tlev(il) + &
      a*( nwp%dat(ix,iy)%tlev(il+1) - nwp%dat(ix,iy)%tlev(il) )
    nwp%dat(ix,iy)%a = a

  !------------------------------------------------------------------
  ! Find the tropopause. 
  !------------------------------------------------------------------
    
    tmin = 99999.0
    do k = nwp%nlevels-5, 1, -1
    
      if (nwp%dat(ix,iy)%tlev(k) < tmin .and. nwp%dat(ix,iy)%plev(k) >= PTOP .and. nwp%dat(ix,iy)%plev(k) <= PBOT) then
         tmin = nwp%dat(ix,iy)%tlev(k)
         nwp%dat(ix,iy)%tropo_level = k
         if (nwp%dat(ix,iy)%tlev(k) - nwp%dat(ix,iy)%tlev(k-1) < 0.5) EXIT
      endif
      
    end do

  !------------------------------------------------------------------
  ! Find the stratopause. 
  !------------------------------------------------------------------
    
    tmax = -99999.0
    do k=nwp%dat(ix,iy)%tropo_level-5, 1, -1
    
      if (nwp%dat(ix,iy)%tlev(k) > tmax) then
         tmax = nwp%dat(ix,iy)%tlev(k)
         nwp%dat(ix,iy)%strato_level = k
         !+1.0 threshold was chosen to account for possible near-isothermal stratosphere
         if (nwp%dat(ix,iy)%tlev(k) - nwp%dat(ix,iy)%tlev(k-1) > 1.0) EXIT
      endif
      
    end do
    
  !------------------------------------------------------------------
  ! Determine if there are any tropospheric temperature inversions. 
  !------------------------------------------------------------------
    ninver = 0
    inverlev(:) = sym%NO
    nwp%dat(ix,iy)%ninversion       = ninver

    do k=nwp%dat(ix,iy)%tropo_level, nwp%dat(ix,iy)%sfc_level-1
       if ((nwp%dat(ix,iy)%tlev(k) - nwp%dat(ix,iy)%tlev(k+1)) > 0.0) then
          if (inverlev(k-1) /= sym%YES) nwp%dat(ix,iy)%ninversion = nwp%dat(ix,iy)%ninversion+1
          inverlev(k) = sym%YES
       endif
    end do
    
    
    !nwp%dat(ix,iy)%inversion_lev(:) = sym%NO
    !nwp%dat(ix,iy)%ninversion       = ninver

    !do k=nwp%dat(ix,iy)%tropo_level, nwp%dat(ix,iy)%sfc_level-1
    !  if ((nwp%dat(ix,iy)%tlev(k) - nwp%dat(ix,iy)%tlev(k+1)) > 0.0) then
    !    if (nwp%dat(ix,iy)%inversion_lev(k-1) /= sym%YES) nwp%dat(ix,iy)%ninversion = nwp%dat(ix,iy)%ninversion+1
    !    nwp%dat(ix,iy)%inversion_lev(k) = sym%YES
    !  endif
    !end do
    
 
! 3. end subroutine
end subroutine compute_levels
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
!~~~~~~~~~~~~~~~~~~~~~ subroutine 8: deallocate nwp arrays ~~~~~~~~~~~~~ 
subroutine deallocate_nwp101_arrays()
 
!-----------------------------------------------------------------------
! !F90 deallocate_nwp101
!
! !Description:
!    This program is to deallocate nwp%dat().
!
! !Input  parameters:
!    ix,iy
!       
! !Output parameters:
!    none
!
!-----------------------------------------------------------------------

! 1. define variables
  integer(kind=4) :: astatus
  integer(kind=4) :: ix, iy
  
  do ix = 1,nwp%nlon
  do iy = 1,nwp%nlat 
  
     call deallocate_nwp_arrays2(ix,iy)
     
  end do
  end do
  
  call deallocate_nwp_arrays1
  
! 3. end subroutine
end subroutine deallocate_nwp101_arrays
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~ subroutine 9: allocate_nwp26_arrays ~~~~~~~~~~~~~~~~
subroutine allocate_nwp26_arrays
 
!-----------------------------------------------------------------------
! !F90 allocate_nwp26_arrays
!
! !Description:
!    This program is to allocate nwp data at satellite observation time.
!
! !Input  parameters:
!    none
!
! !Output parameters:
!    none
!
!-----------------------------------------------------------------------
 
implicit none
 
! 1. define variables
  integer(kind=4) :: astatus
  integer(kind=4) :: x,y,z
  
  if (fylat_nwp_opt == 1 .or. fylat_nwp_opt == 2 .or. fylat_nwp_opt == 4) then
     x = 360
     y = 181
     z = 26
  endif
  
  if (fylat_nwp_opt == 5) then  
     x = 720
     y = 361
     z = 26
  endif
  
! 2. begin program 
  allocate(nwp26%lon(x,y), nwp26%lat(x,y),           & 
           nwp26%plev(x,y,z),                        &
           nwp26%psfc(x,y),nwp26%pmsl(x,y),          &
           nwp26%tsfc(x,y),nwp26%zsfc(x,y),          &
           nwp26%albedo(x,y),nwp26%t_sigma(x,y),     &
           nwp26%rh_sigma(x,y),nwp26%u_sigma(x,y),   &
           nwp26%v_sigma(x,y),nwp26%tpw(x,y),        &
           nwp26%weasd(x,y),nwp26%o3col(x,y),        &
           nwp26%ttropo(x,y),                        &
           nwp26%tlev(x,y,z),                        &
           nwp26%zlev(x,y,z),                        &  
           nwp26%o3lev(x,y,z),                       &    
           nwp26%rhlev(x,y,z),                       &   
           nwp26%wlev(x,y,z),                        &   
           nwp26%clwlev(x,y,z),                      &
           nwp26%ulev(x,y,z),                        &   
           nwp26%vlev(x,y,z),                        &
           stat=astatus)
           
  if (astatus /= 0) then
     print *,"(a,'Not enough memory to allocate nwp26 data structure.')"
     stop
  endif
 
! 3. end subroutine
end subroutine allocate_nwp26_arrays
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~ subroutine 10: deallocate nwpo arrays ~~~~~~~~~~~~~ 
subroutine deallocate_nwp26_arrays
 
!-----------------------------------------------------------------------
! !F90 deallocate_nwp26_arrays
!
! !Description:
!    This program is to deallocate nwp data at satellite observation time.
!
! !Input  parameters::
!    none
!
! !Output parameters:
!    none
!
!-----------------------------------------------------------------------
 
! 1. define variables
  integer(kind=4) :: astatus
 
! 2. begin program 
  deallocate(nwp26%lon, nwp26%lat,          &
             nwp26%plev,                    &
             nwp26%psfc,nwp26%pmsl,         &
             nwp26%tsfc,nwp26%zsfc,         &
             nwp26%albedo,nwp26%t_sigma,    &
             nwp26%rh_sigma,nwp26%u_sigma,  &
             nwp26%v_sigma,nwp26%tpw,       &
             nwp26%weasd,nwp26%o3col,       &
             nwp26%ttropo,                  &
             nwp26%tlev,                    &
             nwp26%zlev,                    &
             nwp26%o3lev,                   &
             nwp26%rhlev,                   &
             nwp26%wlev,                    &   
             nwp26%clwlev,                  &
             nwp26%ulev,                    &   
             nwp26%vlev,                    &
             stat=astatus)
             
  if (astatus /= 0) then
     print *,"(a,'Error deallocating nwp26 data structure.')"
     stop
  endif
 
! 3. end subroutine
end subroutine deallocate_nwp26_arrays
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~ subroutine 9: allocate_nwp26_arrays ~~~~~~~~~~~~~~~~
subroutine allocate_nwp36_arrays
 
!-----------------------------------------------------------------------
! !F90 allocate_nwp26_arrays
!
! !Description:
!    This program is to allocate nwp data at satellite observation time.
!
! !Input  parameters:
!    none
!
! !Output parameters:
!    none
!
!-----------------------------------------------------------------------
 
implicit none
 
! 1. define variables
  integer(kind=4) :: astatus
  integer(kind=4) :: x=1280,y=641,z=36
 
! 2. begin program 
  allocate(nwp36%lon(x,y), nwp36%lat(x,y),           & 
           nwp36%plev(x,y,z),                        &
           nwp36%psfc(x,y),nwp36%pmsl(x,y),          &
           nwp36%tsfc(x,y),nwp36%zsfc(x,y),          &
           nwp36%albedo(x,y),nwp36%t_sigma(x,y),     &
           nwp36%rh_sigma(x,y),nwp36%u_sigma(x,y),   &
           nwp36%v_sigma(x,y),nwp36%tpw(x,y),        &
           nwp36%weasd(x,y),nwp36%o3col(x,y),        &
           nwp36%ttropo(x,y),                        &
           nwp36%tlev(x,y,z),                        &
           nwp36%zlev(x,y,z),                        &  
           nwp36%o3lev(x,y,z),                       &    
           nwp36%rhlev(x,y,z),                       &   
           nwp36%wlev(x,y,z),                        &   
           nwp36%clwlev(x,y,z),                      &
           nwp36%ulev(x,y,z),                        &   
           nwp36%vlev(x,y,z),                        &
           stat=astatus)
           
  if (astatus /= 0) then
     print *,"(a,'Not enough memory to allocate nwp36 data structure.')"
     stop
  endif
 
! 3. end subroutine
end subroutine allocate_nwp36_arrays
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~ subroutine 10: deallocate nwpo arrays ~~~~~~~~~~~~~ 
subroutine deallocate_nwp36_arrays
 
!-----------------------------------------------------------------------
! !F90 deallocate_nwp26_arrays
!
! !Description:
!    This program is to deallocate nwp data at satellite observation time.
!
! !Input  parameters::
!    none
!
! !Output parameters:
!    none
!
!-----------------------------------------------------------------------
 
! 1. define variables
  integer(kind=4) :: astatus
 
! 2. begin program 
  deallocate(nwp36%lon, nwp36%lat,          &
             nwp36%plev,                    &
             nwp36%psfc,nwp36%pmsl,         &
             nwp36%tsfc,nwp36%zsfc,         &
             nwp36%albedo,nwp36%t_sigma,    &
             nwp36%rh_sigma,nwp36%u_sigma,  &
             nwp36%v_sigma,nwp36%tpw,       &
             nwp36%weasd,nwp36%o3col,       &
             nwp36%ttropo,                  &
             nwp36%tlev,                    &
             nwp36%zlev,                    &
             nwp36%o3lev,                   &
             nwp36%rhlev,                   &
             nwp36%wlev,                    &   
             nwp36%clwlev,                  &
             nwp36%ulev,                    &   
             nwp36%vlev,                    &
             stat=astatus)
             
  if (astatus /= 0) then
     print *,"(a,'Error deallocating nwp36 data structure.')"
     stop
  endif
 
! 3. end subroutine
end subroutine deallocate_nwp36_arrays
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~ SUBROUTINE 9: allocate_nwp26_arrays ~~~~~~~~~~~~~~~~
SUBROUTINE allocate_nwp40_arrays
 
!-----------------------------------------------------------------------
! !F90 allocate_nwp26_arrays
!
! !Description:
!    This program is to allocate nwp data at satellite observation time.
!
! !Input  parameters:
!    none
!
! !Output parameters:
!    none
!
!-----------------------------------------------------------------------
 
IMPLICIT NONE
 
! 1. define variables
  INTEGER(kind=4) :: astatus
  INTEGER(kind=4) :: x=1440,y=720,z=40
 
! 2. begin program 
  ALLOCATE(nwp40%lon(x,y), nwp40%lat(x,y),           & 
           nwp40%plev(x,y,z),                        &
           nwp40%psfc(x,y),nwp40%pmsl(x,y),          &
           nwp40%tsfc(x,y),nwp40%zsfc(x,y),          &
           nwp40%albedo(x,y),nwp40%t_sigma(x,y),     &
           nwp40%rh_sigma(x,y),nwp40%u_sigma(x,y),   &
           nwp40%v_sigma(x,y),nwp40%tpw(x,y),        &
           nwp40%weasd(x,y),nwp40%o3col(x,y),        &
           nwp40%ttropo(x,y),                        &
           nwp40%tlev(x,y,z),                        &
           nwp40%zlev(x,y,z),                        &  
           nwp40%o3lev(x,y,z),                       &    
           nwp40%rhlev(x,y,z),                       &   
           nwp40%wlev(x,y,z),                        &   
           nwp40%clwlev(x,y,z),                      &
           nwp40%ulev(x,y,z),                        &
           nwp40%vlev(x,y,z),                        & 
           nwp40%u10m(x,y),                          &
           nwp40%v10m(x,y),                          &   
           stat=astatus)
           
  IF (astatus /= 0) THEN
     PRINT *,"(a,'Not enough memory to allocate nwp40 data structure.')"
     STOP
  ENDIF
 
! 3. END SUBROUTINE
END SUBROUTINE allocate_nwp40_arrays
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~ SUBROUTINE 10: deallocate nwpo arrays ~~~~~~~~~~~~~ 
SUBROUTINE deallocate_nwp40_arrays
 
!-----------------------------------------------------------------------
! !F90 deallocate_nwp26_arrays
!
! !Description:
!    This program is to deallocate nwp data at satellite observation time.
!
! !Input  parameters::
!    none
!
! !Output parameters:
!    none
!
!-----------------------------------------------------------------------
 
! 1. define variables
  INTEGER(kind=4) :: astatus
 
! 2. begin program 
  DEALLOCATE(nwp40%lon, nwp40%lat,          &
             nwp40%plev,                    &
             nwp40%psfc,nwp40%pmsl,         &
             nwp40%tsfc,nwp40%zsfc,         &
             nwp40%albedo,nwp40%t_sigma,    &
             nwp40%rh_sigma,nwp40%u_sigma,  &
             nwp40%v_sigma,nwp40%tpw,       &
             nwp40%weasd,nwp40%o3col,       &
             nwp40%ttropo,                  &
             nwp40%tlev,                    &
             nwp40%zlev,                    &
             nwp40%o3lev,                   &
             nwp40%rhlev,                   &
             nwp40%wlev,                    &   
             nwp40%clwlev,                  &
             nwp40%ulev,                    &   
             nwp40%vlev,                    &
             nwp40%u10m,                    &
             nwp40%v10m,                    &   
             stat=astatus)
             
  IF (astatus /= 0) THEN
     PRINT *,"(a,'Error deallocating nwp40 data structure.')"
     STOP
  ENDIF
 
! 3. END SUBROUTINE
END SUBROUTINE deallocate_nwp40_arrays
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~ SUBROUTINE 9: allocate_nwp26_arrays ~~~~~~~~~~~~~~~~
SUBROUTINE allocate_nwp31_arrays
 
!-----------------------------------------------------------------------
! !F90 allocate_nwp26_arrays
!
! !Description:
!    This program is to allocate nwp data at satellite observation time.
!
! !Input  parameters:
!    none
!
! !Output parameters:
!    none
!
!-----------------------------------------------------------------------
 
IMPLICIT NONE
 
! 1. define variables
  INTEGER(KIND=4) :: astatus
  INTEGER(KIND=4) :: x,y,z

  IF (fylat_nwp_opt == 7 .or. fylat_nwp_opt == 8) THEN ! 9
      x = 1440 !nwp%num_lon !1440
      y = 721 !nwp%num_lat !721
      z = 31
  ENDIF
   
! 2. begin program 
  ALLOCATE(nwp31%lon(x,y), nwp31%lat(x,y),           & 
           nwp31%plev(x,y,z),                        &
           nwp31%psfc(x,y),nwp31%pmsl(x,y),          &
           nwp31%tsfc(x,y),nwp31%zsfc(x,y),          &
           nwp31%albedo(x,y),nwp31%t_sigma(x,y),     &
           nwp31%rh_sigma(x,y),nwp31%u_sigma(x,y),   &
           nwp31%v_sigma(x,y),nwp31%tpw(x,y),        &
           nwp31%weasd(x,y),nwp31%o3col(x,y),        &
           nwp31%ttropo(x,y),                        &
           nwp31%tlev(x,y,z),                        &
           nwp31%zlev(x,y,z),                        &  
           nwp31%o3lev(x,y,z),                       &    
           nwp31%rhlev(x,y,z),                       &   
           nwp31%wlev(x,y,z),                        &   
           nwp31%clwlev(x,y,z),                      &
           nwp31%ulev(x,y,z),                        &
           nwp31%vlev(x,y,z),                        &  
           stat=astatus)
           
  IF (astatus /= 0) THEN
     PRINT *,"(a,'Not enough memory to allocate nwp31 data structure.')"
     STOP
  ENDIF
 
! 3. END SUBROUTINE
END SUBROUTINE allocate_nwp31_arrays
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~ SUBROUTINE 10: deallocate nwpo arrays ~~~~~~~~~~~~~ 
SUBROUTINE deallocate_nwp31_arrays
 
!-----------------------------------------------------------------------
! !F90 deallocate_nwp31_arrays
!
! !Description:
!    This program is to deallocate nwp data at satellite observation time.
!
! !Input  parameters::
!    none
!
! !Output parameters:
!    none
!
!-----------------------------------------------------------------------
 
! 1. define variables
  INTEGER(KIND=4) :: astatus
 
! 2. begin program 
  DEALLOCATE(nwp31%lon, nwp31%lat,          &
             nwp31%plev,                    &
             nwp31%psfc,nwp31%pmsl,         &
             nwp31%tsfc,nwp31%zsfc,         &
             nwp31%albedo,nwp31%t_sigma,    &
             nwp31%rh_sigma,nwp31%u_sigma,  &
             nwp31%v_sigma,nwp31%tpw,       &
             nwp31%weasd,nwp31%o3col,       &
             nwp31%ttropo,                  &
             nwp31%tlev,                    &
             nwp31%zlev,                    &
             nwp31%o3lev,                   &
             nwp31%rhlev,                   &
             nwp31%wlev,                    &   
             nwp31%clwlev,                  &
             nwp31%ulev,                    &
             nwp31%vlev,                    &
             stat=astatus)
             
  IF (astatus /= 0) THEN
     PRINT *,"(a,'Error deallocating nwp31 data structure.')"
     STOP
  ENDIF
 
! 3. END SUBROUTINE
END SUBROUTINE deallocate_nwp31_arrays
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!~~~~~~~~~~~~~~~~~~~ SUBROUTINE 9: allocate_nwp26_arrays ~~~~~~~~~~~~~~~~
SUBROUTINE allocate_nwp41_arrays
 
!-----------------------------------------------------------------------
! !F90 allocate_nwp41_arrays
!
! !Description:
!    This program is to allocate nwp data at satellite observation time.
!
! !Input  parameters:
!    none
!
! !Output parameters:
!    none
!
!-----------------------------------------------------------------------
 
IMPLICIT NONE
 
! 1. define variables
  INTEGER(KIND=4) :: astatus
  INTEGER(KIND=4) :: x,y,z

  IF (fylat_nwp_opt == 9) THEN ! 9
      x = 720 !nwp%num_lon !1440
      y = 361 !nwp%num_lat !721
      z = 41
  ENDIF

  IF (fylat_nwp_opt == 10) THEN ! 9
      x = 1440 !nwp%num_lon !1440
      y = 721 !nwp%num_lat !721
      z = 41
  ENDIF
   
! 2. begin program 
  ALLOCATE(nwp41%lon(x,y), nwp41%lat(x,y),           & 
           nwp41%plev(x,y,z),                        &
           nwp41%psfc(x,y),nwp41%pmsl(x,y),          &
           nwp41%tsfc(x,y),nwp41%zsfc(x,y),          &
           nwp41%albedo(x,y),nwp41%t_sigma(x,y),     &
           nwp41%rh_sigma(x,y),nwp41%u_sigma(x,y),   &
           nwp41%v_sigma(x,y),nwp41%tpw(x,y),        &
           nwp41%weasd(x,y),nwp41%o3col(x,y),        &
           nwp41%ttropo(x,y),                        &
           nwp41%tlev(x,y,z),                        &
           nwp41%zlev(x,y,z),                        &  
           nwp41%o3lev(x,y,z),                       &    
           nwp41%rhlev(x,y,z),                       &   
           nwp41%wlev(x,y,z),                        &   
           nwp41%clwlev(x,y,z),                      &
           nwp41%ulev(x,y,z),                        &
           nwp41%vlev(x,y,z),                        &  
           stat=astatus)
           
  IF (astatus /= 0) THEN
     PRINT *,"(a,'Not enough memory to allocate nwp41 data structure.')"
     STOP
  ENDIF
 
! 3. END SUBROUTINE
END SUBROUTINE allocate_nwp41_arrays
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~ SUBROUTINE 10: deallocate nwpo arrays ~~~~~~~~~~~~~ 
SUBROUTINE deallocate_nwp41_arrays
 
!-----------------------------------------------------------------------
! !F90 deallocate_nwp41_arrays
!
! !Description:
!    This program is to deallocate nwp data at satellite observation time.
!
! !Input  parameters::
!    none
!
! !Output parameters:
!    none
!
!-----------------------------------------------------------------------
 
! 1. define variables
  INTEGER(KIND=4) :: astatus
 
! 2. begin program 
  DEALLOCATE(nwp41%lon, nwp41%lat,          &
             nwp41%plev,                    &
             nwp41%psfc,nwp41%pmsl,         &
             nwp41%tsfc,nwp41%zsfc,         &
             nwp41%albedo,nwp41%t_sigma,    &
             nwp41%rh_sigma,nwp41%u_sigma,  &
             nwp41%v_sigma,nwp41%tpw,       &
             nwp41%weasd,nwp41%o3col,       &
             nwp41%ttropo,                  &
             nwp41%tlev,                    &
             nwp41%zlev,                    &
             nwp41%o3lev,                   &
             nwp41%rhlev,                   &
             nwp41%wlev,                    &   
             nwp41%clwlev,                  &
             nwp41%ulev,                    &
             nwp41%vlev,                    &
             stat=astatus)
             
  IF (astatus /= 0) THEN
     PRINT *,"(a,'Error deallocating nwp41 data structure.')"
     STOP
  ENDIF
 
! 3. END SUBROUTINE
END SUBROUTINE deallocate_nwp41_arrays
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~ subroutine 11: climoz_101_o3 ~~~~~~~~~~~~~~~~~~~~~~
subroutine climoz_101_o3(rlat, month, omix)
 
!-----------------------------------------------------------------------
! !F90 climoz_101_o3
!
! !Description:
!    This program is to get climatological o3 profile.
!c * Obtain climatological ozone mixing-ratio profile
!c     by interpolating in latitude and month amongst
!c     six LBLRTM model atmospheres.
!c .... version of 02.08.00
!c * Pressure coordinate is 101-level AIRS SPACECRAFT
!
! !Input  parameters:
!    rlat            = real latitude (deg +N,-S)
!    month           = INTEGER month (1,...,12)
!
! !Output parameters:
!    omix            = real ozone profile (ppmv)
!
!-----------------------------------------------------------------------

implicit none

real(kind=4),intent(in) :: rlat
integer(kind=4), intent(in) :: month    
integer(kind=4), parameter :: nl=101,ns=2,nz=3
real(kind=4) :: tlat(nz),ozmr(nl,nz,ns),omr(nl,ns)
real(kind=4), intent(out) :: omix(nl)
integer(kind=4) :: m,n, i, j, k
integer(kind=4) :: jl, kk, jl1, jl2, nmon, imon
real(kind=4)   :: wt1,wt2, or1,or2, alat
data tlat/15.,40.,65./

! * TROPICAL
  DATA (ozmr(m,1,1),m=1,nl)/0.47628,0.26429,0.26278,0.45250,0.75116,     &
       1.10347,1.55742,1.96852,2.55758,3.09206,3.88655,4.78962,5.83901,  &
       6.91912,7.84273,8.61191,9.17298,9.48285,9.75589,9.69629,9.62726,  &
       9.33416,8.94015,8.56266,7.99507,7.36964,6.70970,5.98273,5.25233,  &
       4.50810,3.92803,3.42163,2.93116,2.45571,1.99446,1.62320,1.40553,  &
       1.19122,0.97212,0.77005,0.60182,0.43786,0.33360,0.25389,0.19479,  &
       0.16141,0.13792,0.13077,0.12378,0.11626,0.10909,0.10378,0.09858,  &
       0.09349,0.08850,0.08364,0.07914,0.07473,0.07041,0.06617,0.06285,  &
       0.05982,0.05684,0.05392,0.05204,0.05020,0.04839,0.04689,0.04561,  &
       0.04435,0.04332,0.04256,0.04183,0.04111,0.04041,0.03973,0.03908,  &
       0.03844,0.03780,0.03719,0.03660,0.03602,0.03571,0.03544,0.03525,  &
       0.03503,0.03464,0.03426,0.03388,0.03341,0.03294,0.03246,0.03191,  &
       0.03137,0.03076,0.03005,0.02935,0.02867,0.02798,0.02780,0.02780/ 
       
! * MIDLATITUDE WINTER
  DATA (ozmr(m,2,1),m=1,nl)/0.51382,0.24549,0.31340,0.49141,0.80263,     &
       1.21291,1.76263,2.26416,3.07810,3.81664,4.71840,5.60591,6.36097,  &
       6.82852,7.09197,7.16580,7.11413,6.99021,6.86566,6.61458,6.37030,  &
       6.15518,5.95848,5.77003,5.60211,5.44577,5.28917,5.12866,4.92166,  &
       4.65215,4.37542,4.10128,3.83576,3.57837,3.32867,3.05455,2.72902,  & 
       2.41998,2.14892,1.89915,1.69178,1.48968,1.33095,1.19127,1.07451,  &
       0.98597,0.91318,0.86575,0.81938,0.78001,0.73926,0.68011,0.62221,  &
       0.56549,0.50992,0.45592,0.40760,0.36020,0.31372,0.26811,0.23664,  &
       0.20922,0.18231,0.15587,0.13995,0.12431,0.10893,0.09822,0.09083,  &
       0.08356,0.07689,0.07090,0.06551,0.06044,0.05575,0.05266,0.04953,  &
       0.04644,0.04340,0.04050,0.03794,0.03542,0.03428,0.03334,0.03241,  &
       0.03151,0.03063,0.02976,0.02891,0.02863,0.02840,0.02820,0.02811,  &
       0.02801,0.02794,0.02788,0.02783,0.02779,0.02774,0.02773,0.02773/ 
       
! * SUBARCTIC WINTER
  DATA (ozmr(m,3,1),m=1,nl)/0.66449,0.26978,0.48593,0.67688,0.95783,     &
       1.35182,1.91218,2.41775,3.09378,3.70718,4.53368,5.24354,5.74008,  &
       6.04019,6.19492,6.21557,6.15697,6.02288,5.88833,5.68564,5.48967,  & 
       5.30404,5.12722,4.95782,4.86320,4.79873,4.73519,4.67163,4.60720,  &
       4.54084,4.40541,4.24629,4.09218,3.94278,3.79785,3.59108,3.26675,  &
       2.94742,2.62079,2.32052,2.07317,1.83210,1.64421,1.47960,1.33092,  & 
       1.20137,1.07704,0.95937,0.84432,0.74373,0.65183,0.59004,0.52954,  &
       0.47028,0.41223,0.35871,0.33964,0.32094,0.30259,0.28460,0.25715,  &
       0.22765,0.19868,0.17023,0.14697,0.12411,0.10165,0.08725,0.07858,  &
       0.07006,0.06190,0.05415,0.04896,0.04499,0.04143,0.03971,0.03815,  &
       0.03665,0.03516,0.03373,0.03241,0.03110,0.02987,0.02866,0.02752,  &
       0.02643,0.02540,0.02439,0.02339,0.02274,0.02213,0.02152,0.02092,  &
       0.02033,0.01973,0.01915,0.01857,0.01800,0.01744,0.01729,0.01729/ 
       
! * TROPICAL ... identical to ozmr(i,1,1), repeated for symmetry
  DATA (ozmr(m,1,2),m=1,nl)/0.47628,0.26429,0.26278,0.45250,0.75116,     &
       1.10347,1.55742,1.96852,2.55758,3.09206,3.88655,4.78962,5.83901,  &
       6.91912,7.84273,8.61191,9.17298,9.48285,9.75589,9.69629,9.62726,  &
       9.33416,8.94015,8.56266,7.99507,7.36964,6.70970,5.98273,5.25233,  &
       4.50810,3.92803,3.42163,2.93116,2.45571,1.99446,1.62320,1.40553,  &
       1.19122,0.97212,0.77005,0.60182,0.43786,0.33360,0.25389,0.19479,  &
       0.16141,0.13792,0.13077,0.12378,0.11626,0.10909,0.10378,0.09858,  &
       0.09349,0.08850,0.08364,0.07914,0.07473,0.07041,0.06617,0.06285,  &
       0.05982,0.05684,0.05392,0.05204,0.05020,0.04839,0.04689,0.04561,  &
       0.04435,0.04332,0.04256,0.04183,0.04111,0.04041,0.03973,0.03908,  &
       0.03844,0.03780,0.03719,0.03660,0.03602,0.03571,0.03544,0.03525,  &
       0.03503,0.03464,0.03426,0.03388,0.03341,0.03294,0.03246,0.03191,  &
       0.03137,0.03076,0.03005,0.02935,0.02867,0.02798,0.02780,0.02780/ 
       
! * MIDLATITUDE SUMMER
  DATA (ozmr(m,2,2),m=1,nl)/0.52876,0.21784,0.28626,0.50723,0.82434,     &
       1.15971,1.49190,1.79766,2.36209,2.87423,3.60218,4.42858,5.42260,  &
       6.51285,7.45312,8.23323,8.70375,8.74203,8.74005,8.36522,7.99382,  &
       7.58434,7.16521,6.76365,6.39723,6.05246,5.70201,5.33505,4.93858,  &
       4.49900,4.10257,3.73019,3.36953,3.01991,2.68074,2.38304,2.15283,  &
       1.92043,1.66023,1.41371,1.19088,0.97372,0.83431,0.72694,0.64802,  &
       0.60467,0.56397,0.52707,0.49100,0.46648,0.44137,0.39970,0.35890,  &
       0.31893,0.27978,0.24321,0.22521,0.20755,0.19024,0.17325,0.15940,  &
       0.14655,0.13393,0.12153,0.11510,0.10879,0.10258,0.09697,0.09181,  &
       0.08674,0.08211,0.07799,0.07393,0.06992,0.06608,0.06282,0.06019,  &
       0.05773,0.05530,0.05302,0.05107,0.04914,0.04740,0.04571,0.04406,  &
       0.04250,0.04107,0.03966,0.03827,0.03717,0.03611,0.03507,0.03419,  &
       0.03332,0.03248,0.03169,0.03092,0.03015,0.02939,0.02918,0.02918/ 
       
! * SUBARCTIC SUMMER
  DATA (ozmr(m,3,2),m=1,nl)/0.60245,0.21936,0.27680,0.48346,0.76740,     &
       1.06602,1.37644,1.65892,2.09879,2.49791,3.19647,3.97613,4.89203,  &
       5.89915,6.74986,7.42939,7.75508,7.68046,7.57369,7.19015,6.81393,  & 
       6.42741,6.04403,5.67673,5.46090,5.30670,5.11691,4.86912,4.63322,  &
       4.40961,4.11937,3.80935,3.50909,3.21802,2.93565,2.63975,2.31163,  &
       2.01203,1.79760,1.59367,1.40729,1.22565,1.09532,0.98670,0.89692,  & 
       0.83042,0.76869,0.71395,0.66044,0.61732,0.57697,0.53986,0.50352,  &
       0.46793,0.43306,0.39843,0.35992,0.32216,0.28512,0.24878,0.22238,  &
       0.19888,0.17580,0.15314,0.13730,0.12173,0.10643,0.09639,0.09009,  &
       0.08390,0.07884,0.07512,0.07161,0.06823,0.06486,0.06128,0.05810,  &
       0.05504,0.05203,0.04924,0.04697,0.04474,0.04308,0.04151,0.03996,  &
       0.03849,0.03717,0.03586,0.03458,0.03342,0.03229,0.03116,0.03001,  &
       0.02887,0.02768,0.02642,0.02518,0.02395,0.02274,0.02241,0.02241/ 

alat=abs(rlat)

if (alat <= 15.) then

    jl=1
    kk=1
    do n=1,kk
    do m=1,nl
       omr(m,n)=ozmr(m,jl,n)
    end do
    end do

    do m=1,nl
       omix(m)=omr(m,1)
    end do
    
endif

if (alat > 15. .and. alat < 65.) then  !!!!!!!!!!!

    jl1=1
    if (alat > 40.) jl1=2
    jl2=jl1+1
    wt1=(tlat(jl2)-alat)/25.
    wt2=1.-wt1
    
    do n=1,ns
       do m=1,nl
          or1=ozmr(m,jl1,n)
          or2=ozmr(m,jl2,n)
          omr(m,n)=wt1*or1+wt2*or2
       end do
    end do

    nmon=month
    !c * shift 6 months for Southern Hemisphere
    if(rlat < 0.) nmon=nmon+6
    if(nmon > 12) nmon=nmon-12
    imon=abs(nmon-7)
    wt1=float(imon)/6.
    wt2=1.-wt1
    
    do m=1,nl
       omix(m)=wt1*omr(m,1)+wt2*omr(m,2)
    end do

endif

if (alat >= 65.) then  !!!!!!!!!!!!

    jl=3
    kk=2
    do n=1,kk
       do m=1,nl
          omr(m,n)=ozmr(m,jl,n)
       end do
    end do
    nmon=month
    !c * shift 6 months for Southern Hemisphere
    if(rlat < 0.) nmon=nmon+6
    if(nmon > 12) nmon=nmon-12
    
    imon=abs(nmon-7)
    wt1=float(imon)/6.
    wt2=1.-wt1
    
    do m=1,nl
       omix(m)=wt1*omr(m,1)+wt2*omr(m,2)
    end do

endif
        
! 3. end subroutine
end subroutine climoz_101_o3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~ function 12: determine_nwp_quality_flag ~~~~~~~~~~~
function determine_nwp_quality_flag(tsfc, ilon, ilat, nlon, nlat, &
                                    nlon_uni, nlat_uni) RESULT(tsfc_uni)

!-----------------------------------------------------------------------
! !F90 determine_nwp_quality_flag
!
! !Description:
!    This program is to determine nwp quality flag.
!
! !Input  parameters::
!    none
!
! !Output parameters:
!    none
!
!-----------------------------------------------------------------------

  real(kind=4), dimension(:,:), intent(in) :: tsfc
  integer(kind=4) :: ilat, ilon, nlat, nlon, nlat_uni, nlon_uni
  
  integer(kind=4) :: ilat_start, ilat_end, ilon_start, ilon_end
  real(kind=4) :: tsfc_min, tsfc_max, tsfc_uni
  
  ilon_start = max(1,ilon-nlon_uni)
  ilon_end = min(nlon,ilon+nlon_uni)
  ilat_start = max(1,ilat-nlat_uni)
  ilat_end = min(nlat,ilat+nlat_uni)
 
  tsfc_min = minval(tsfc(ilon_start:ilon_end,ilat_start:ilat_end))
  tsfc_max = maxval(tsfc(ilon_start:ilon_end,ilat_start:ilat_end))
  tsfc_uni = (tsfc_max - tsfc_min) / 3.0
  
  return

! 3. end subroutine
end function determine_nwp_quality_flag
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~ subroutine 13: allocate nwp arrays1 ~~~~~~~~~~~~~~~
subroutine allocate_nwp_arrays1
 
!-----------------------------------------------------------------------
! !F90 allocate_nwp_arrays1
!
! !Description:
!    This program is to allocate nwp data at satellite observation time.
!
! !Input  parameters:
!    none
!
! !Output parameters:
!    none
!
!-----------------------------------------------------------------------
 
implicit none
 
! 1. define variables
  integer(kind=4) :: astatus
 
! 2. begin program 
  if (fylat_nwp_opt == 1 .or. fylat_nwp_opt == 2 .or. fylat_nwp_opt == 4) then   ! 1=ncep and 2=gfs
  
    allocate(nwp%dat(nwp%nlon,nwp%nlat), stat=astatus)
  
    if (astatus /= 0) then
       print *,"(a,'Not enough memory to allocate nwp dat data structure.')"
       stop
    endif
    
  endif

  if ( fylat_nwp_opt == 5) then   ! 1=ncep and 2=gfs
  
    allocate(nwp%dat(nwp%nlon05,nwp%nlat05), stat=astatus)
  
    if (astatus /= 0) then
       print *,"(a,'Not enough memory to allocate nwp dat data structure.')"
       stop
    endif
    
  endif
 
  if (fylat_nwp_opt == 3) then   ! 3=T639
  
     allocate(nwp%dat(nwp%nlon_T639,nwp%nlat_T639), stat=astatus)
  
     if (astatus /= 0) then
        print *,"(a,'Not enough memory to allocate nwp (T639) dat data structure.')"
        stop
     endif     
      
  endif

  if (fylat_nwp_opt == 6) then   ! 6=grapes
  
     allocate(nwp%dat(nwp%nlon25,nwp%nlat25), stat=astatus)
  
     if (astatus /= 0) then
        print *,"(a,'Not enough memory to allocate nwp (grapes) dat data structure.')"
        stop
     endif     
      
  endif

  if (fylat_nwp_opt == 7 .or. fylat_nwp_opt == 8) then   ! 8=gfs0p25 grib2
  
    allocate(nwp%dat(nwp%nlon25_gfs,nwp%nlat25_gfs), stat=astatus)
  
    if (astatus /= 0) then
       print *,"(a,'Not enough memory to allocate nwp 0p50 dat data structure.')"
       stop
    endif
    
  endif

  if (fylat_nwp_opt == 9) then   ! 9=gfs0p50 grib2 @41 layers
  
    allocate(nwp%dat(nwp%nlon05,nwp%nlat05), stat=astatus)
  
    if (astatus /= 0) then
       print *,"(a,'Not enough memory to allocate nwp 0p50 dat data structure.')"
       stop
    endif
    
  endif  

  if (fylat_nwp_opt == 10) then   ! 10=gfs0p25 grib2 @41 layers
  
    allocate(nwp%dat(nwp%nlon25_gfs,nwp%nlat25_gfs), stat=astatus)
  
    if (astatus /= 0) then
       print *,"(a,'Not enough memory to allocate nwp 0p50 dat data structure.')"
       stop
    endif
    
  endif    
  
! 3. end subroutine
end subroutine allocate_nwp_arrays1
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~ subroutine 14: deallocate nwp arrays1 ~~~~~~~~~~~~~ 
subroutine deallocate_nwp_arrays1
 
!-----------------------------------------------------------------------
! !F90 deallocate_nwp_arrays
!
! !Description:
!    This program is to deallocate nwp data at satellite observation time.
!
! !Input  parameters::
!    none
!
! !Output parameters:
!    none
!
!-----------------------------------------------------------------------
 
! 1. define variables
  integer(kind=4) :: astatus
 
! 2. begin program 
  deallocate(nwp%dat, stat=astatus)
           
  if (astatus /= 0) then
     print *,"(a,'Error deallocating nwp dat data structure.')"
     stop
  endif
 
! 3. end subroutine
end subroutine deallocate_nwp_arrays1
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~ subroutine 15: allocate nwp arrays2 ~~~~~~~~~~~~~~~
subroutine allocate_nwp_arrays2(ix,iy)
 
!-----------------------------------------------------------------------
! !F90 allocate_nwp_arrays
!
! !Description:
!    This program is to allocate nwp data at satellite observation time.
!
! !Input  parameters:
!    none
!
! !Output parameters:
!    none
!
!-----------------------------------------------------------------------
 
implicit none
 
! 1. define variables
  integer(kind=4) :: astatus
  integer(kind=4) :: ix, iy
 
! 2. begin program 
  
  allocate(nwp%dat(ix,iy)%plev(nwp%nlevels),          & 
           nwp%dat(ix,iy)%tlev(nwp%nlevels),          &  
           nwp%dat(ix,iy)%zlev(nwp%nlevels),          &  
           nwp%dat(ix,iy)%o3lev(nwp%nlevels),         &  
           nwp%dat(ix,iy)%wlev(nwp%nlevels),          &  
           nwp%dat(ix,iy)%tpwlev(nwp%nlevels),        &  
           nwp%dat(ix,iy)%rhlev(nwp%nlevels),         & 
           nwp%dat(ix,iy)%ulev(nwp%nlevels),          &  
           nwp%dat(ix,iy)%vlev(nwp%nlevels),          & 
           nwp%dat(ix,iy)%inversion_lev(nwp%nlevels), & 
           stat=astatus)
           
  if (astatus /= 0) then
     print *,"(a,'Not enough memory to allocate nwp data structure.')"
     stop
  endif
  
! 3. end subroutine
end subroutine allocate_nwp_arrays2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~ subroutine 16: deallocate nwp arrays2 ~~~~~~~~~~~~~ 
subroutine deallocate_nwp_arrays2(ix, iy)
 
!-----------------------------------------------------------------------
! !F90 deallocate_nwp_arrays2
!
! !Description:
!    This program is to deallocate nwp data at satellite observation time.
!
! !Input  parameters::
!    none
!
! !Output parameters:
!    none
!
!-----------------------------------------------------------------------
 
! 1. define variables
  integer(kind=4) :: astatus
  integer(kind=4) :: ix, iy
   
! 2. begin program 
  deallocate(nwp%dat(ix,iy)%plev,          & 
             nwp%dat(ix,iy)%tlev,          &  
             nwp%dat(ix,iy)%zlev,          &  
             nwp%dat(ix,iy)%o3lev,         &  
             nwp%dat(ix,iy)%wlev,          &  
             nwp%dat(ix,iy)%tpwlev,        &  
             nwp%dat(ix,iy)%rhlev,         &
             nwp%dat(ix,iy)%ulev,          &  
             nwp%dat(ix,iy)%vlev,          &               
             nwp%dat(ix,iy)%inversion_lev, & 
             stat=astatus)
           
  if (astatus /= 0) then
     print *,"(a,'Error deallocating nwp data structure.')"
     stop
  endif
 
! 3. end subroutine
end subroutine deallocate_nwp_arrays2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE O3_to_Dobson(p,o3,dobson)
 
!-----------------------------------------------------------------------
! !F90 convert O3 ppmv to Dobson unit
!
! !Description:
! Purpose:
!   Compute column ozone ammount in dobson units (DU), from a
!   profile of O3 mixing ratio vs. pressure.
!   1 DU = 2.69e16 molecules cm-2
!
! !Input  parameters::
!    p.......pressure, mb
!    o3......ozone mixing ratio,  ppmv
!
! !Output parameters:
!    dobs....column o3 abundance, DU
!
!-----------------------------------------------------------------------
 
! 1. define variables
  INTEGER(kind=4) :: i,j,n
  REAL(kind=4),INTENT(in),DIMENSION(:) :: p
  REAL(kind=4),INTENT(in),DIMENSION(:) :: O3
  REAL(kind=4)  ::  dobson, term, g, const, mdry, avagadro
!- some constants

   g        = STANDARD_GRAVITY  !9.8         ; gravity, m/s2
   avagadro = AVOGADRO_CONSTANT != 6.02214199e+23_fp!6.02252e23  ; molec/mol
   mdry     = 0.028964    !; molec. wt. of dry air, kg/mol

   const    = 0.01 * avagadro / (g * mdry)

!   n        = n_elements(p)

   n = size(p)

   dobson     = 0.

!- sum o3 over height

   DO j = 1, n-1 

     term = 0.5*(o3(j) + o3(j-1)) *1.e-6 * abs(p(j-1) - p(j)) * const/2.69e16

     dobson = dobson + term 

   ENDDO
   
! 3. END SUBROUTINE
END SUBROUTINE O3_to_Dobson
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!+++++++++++++++++++++ step 3: end module++++++++++++++++++++++++++++++
end module nwp_utils_module
