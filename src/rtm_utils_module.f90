module rtm_utils_module

!-----------------------------------------------------------------------
! !F90                                                                  
!
! !Description: 
!    This module is to deal with rtm data.
!
! !Input parameters
!    none
! 
! !Output parameters
!    none
!
!  
! !end
!----------------------------------------------------------------------

use data_arrays_module
use constant
use planck_module
use names_module
use rtm_tran_module
use numerical

implicit none

!+++++++++++++++++++ step 1: define global variables +++++++++++++++++++
! |------|
! | none |
! |------|
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  
contains
!+++++++++++++++++++ step 2: subroutines +++++++++++++++++++++++++++++++

!~~~~~~~~~~~~~~~~~~~ subroutine 1: ir_rtm_driver ~~~~~~~~~~~~~~~~~~~~~~
subroutine ir_rtm_driver()

!-----------------------------------------------------------------------
! !F90 deal_rtm_begin
!
! !Description:
!    The begin program for deal with rtm data.[allocate variables]
!
! !Input  parameters:
!    none
!
! !Output parameters:
!    none
!
!-----------------------------------------------------------------------


! 1. define variables
real(kind=4) :: a, cos_satzen
integer:: pixel_rtm_option = 1
integer :: i,j,k
integer(kind=1) wflag

! 2. begin program
  print*,'  ... run ir fast rtm  '

!%%%%%%%%%%%%%%%
! STEP 1. 
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
  
!=== 2.1. allocate rtm
if (fylat_nwp_opt == 1 .or. fylat_nwp_opt == 2 .or. fylat_nwp_opt == 4) then   ! 1=ncep and 2=gfs
   call allocate_rtm(nwp%nlon,nwp%nlat)
endif

if (fylat_nwp_opt == 5) then   ! 5=gfs0p50 grib2
   call allocate_rtm(nwp%nlon05,nwp%nlat05)
endif

if (fylat_nwp_opt == 3) then   ! 3=T639
   call allocate_rtm(nwp%nlon_T639,nwp%nlat_T639)
endif

if (fylat_nwp_opt == 6) then   ! 6=grapes
   call allocate_rtm(nwp%nlon25,nwp%nlat25)
endif

if (fylat_nwp_opt == 7 .or. fylat_nwp_opt == 8) then   ! 8=gfs 0p25 grib2 1440*721
   call allocate_rtm(nwp%nlon25_gfs,nwp%nlat25_gfs)
endif

if (fylat_nwp_opt == 9) then   ! 9=gfs 0p50 grib2 720*361 @41 layers
   call allocate_rtm(nwp%nlon05,nwp%nlat05)
endif

if (fylat_nwp_opt == 10) then   ! 10=gfs 0p25 grib2 1440*721  @41 layers
   call allocate_rtm(nwp%nlon25_gfs,nwp%nlat25_gfs)
endif

!----------------------------------------------
! Determine the number of viewing angle bins.
!---------------------------------------------- 
nwp%rtm_nvzen = (int(ceiling(1.0 / RTM_VZA_BINSIZE))) + 1

!=== 2.2. cycle begin 

do i = 1, sat%nElem
do j = 1, sat%nLine

!do i = 20, 20
!do j = 490, 498
      
      cos_satzen = cos(geo%SensorZenith(i,j)*DTOR)
      sat%ivza(i,j) = int( (cos_satzen / RTM_VZA_BINSIZE) ) + 1
      !print*,i,j,'ddsdsd ===',cos_satzen, sat%ivza(i,j) ,geo%SensorZenith(i,j)
     !print*,geo%lon(i,j),geo%lat(i,j),geo%cos_satzen(i,j)
     !print*,i,j,rtm(sat%x_nwp(i,j),sat%y_nwp(i,j))%flag, sat%ivza(i,j),sat%x_nwp(i,j),sat%y_nwp(i,j),nwp%rtm_nvzen
      if (rtm(sat%x_nwp(i,j),sat%y_nwp(i,j))%flag == 0) then
         call allocate_rtm_cell(sat%x_nwp(i,j),sat%y_nwp(i,j),nwp%rtm_nvzen)
         !print*,'22',rtm(sat%x_nwp(i,j),sat%y_nwp(i,j))%d(:)%flag
      endif

!print*,sat%ivza(i,j)
!print*,rtm(sat%x_nwp(i,j),sat%y_nwp(i,j))%d(:)%flag
!print*,i,j,sat%x_nwp(i,j),sat%y_nwp(i,j),rtm(sat%x_nwp(i,j),sat%y_nwp(i,j))%d(sat%ivza(i,j))%flag
!rtm(sat%x_nwp(i,j),sat%y_nwp(i,j))%d(sat%ivza(i,j))%flag  = 0
!print*,'---'
!print*,'dddd'
!print*,sat%satzen(i,j) ,sat%ivza(i,j),sat%x_nwp(i,j),sat%y_nwp(i,j)
!print*,'xxxxxx',i,j,sat%zsfc(i,j), 'vvv',nwp%dat(266,18)%sfc_level,'v2',nwp%dat(i,j)%sfc_level
!print*,rtm(sat%x_nwp(i,j),sat%y_nwp(i,j))%d(sat%ivza(i,j))%satzen
! print*,'xxxxx11'
!      if (i == 1000 .and. j == 1000) then
!         print*,'ggg1-----',i,j,geo%lon(i,j),geo%lat(i,j),sat%x_nwp(i,j),sat%y_nwp(i,j),sat%ivza(i,j)
!         print*, 'tt',nwp%dat(sat%x_nwp(i,j),sat%y_nwp(i,j))%tlev(:)
!         print*,'ww',nwp%dat(sat%x_nwp(i,j),sat%y_nwp(i,j))%wlev(:)
!      endif
      !print*,'rtm0',rtm(sat%x_nwp(i,j),sat%y_nwp(i,j))%d(sat%ivza(i,j))%flag
      if (rtm(sat%x_nwp(i,j),sat%y_nwp(i,j))%d(sat%ivza(i,j))%flag == 0) then
      !print*,'rtm1',rtm(sat%x_nwp(i,j),sat%y_nwp(i,j))%d(sat%ivza(i,j))%flag
         call vza_mid(sat%ivza(i,j), RTM_VZA_BINSIZE, rtm(sat%x_nwp(i,j),sat%y_nwp(i,j))%d(sat%ivza(i,j))%satzen)   
         !print*,'rtm2'     
         ! call rtm
         if (fylat_rtm_opt == 1) then ! rtm = 1 means use PFAAST modle
         
            if (fylat_sensor_id == 1 .or. fylat_sensor_id == 2) then  ! Aqua-modis2mersiii 
              ! print*,'run plod +++++++++++++++'
               call run_plod_mersiII_modis(sat%x_nwp(i,j),sat%y_nwp(i,j),sat%ivza(i,j),nwp%nlevels,pixel_rtm_option)
            endif
            
            if (fylat_sensor_id == 21) then  ! real fy3d-mersiii rtm
               ! print*,'run plod +++++++++++++++'
               call run_plod_fy3d_mersi_ii(sat%x_nwp(i,j),sat%y_nwp(i,j),sat%ivza(i,j),nwp%nlevels,pixel_rtm_option)
               !print*,'rtm3'
            endif 
                       
         endif
         
      endif

      !---- note --------------------
      ! nwp%nlevels = 101 have beed defined in nwp_utils_modules and saved in data_arrays_moduls 
      ! sat%zsfc is elevation data
      !------------------------------
  !    print*,'11',geo%SolarZenith(i,j),sat%ivza(i,j),sat%isfc(i,j)
  !        integer(KIND=int1) :: SHALLOW_OCEAN         ! 0
  !  integer(KIND=int1) :: LAND                  ! 1
  !  integer(KIND=int1) :: COASTLINE             ! 2
  !  integer(KIND=int1) :: SHALLOW_INLAND_WATER  ! 3
  !  integer(KIND=int1) :: EPHEMERAL_WATER       ! 4
  !  integer(KIND=int1) :: DEEP_INLAND_WATER     ! 5
  !  integer(KIND=int1) :: MODERATE_OCEAN        ! 6
  !  integer(KIND=int1) :: DEEP_OCEAN            ! 7
      wflag = 1
      if (geo%lsm(i,j) == 1 .or. geo%lsm(i,j) == 2) then  ! find land/sea mask to use sst from oisst
         wflag = 0
      endif 
      !print*,'ggg1-----', 'h',sat%zsfc(i,j), sat%isfc(i,j)
      call get_pixel_sfc_level(nwp%nlevels, nwp%dat(sat%x_nwp(i,j),sat%y_nwp(i,j))%zlev(:), sat%zsfc(i,j), sat%isfc(i,j), a)
      !print*,'ggg2-----', 'dd',geo%lsm(i,j),sat%x_nwp(i,j),sat%y_nwp(i,j),sat%ivza(i,j),sat%isfc(i,j),i,j,a, wflag
      call get_pixel_clear_rad(sat%x_nwp(i,j),sat%y_nwp(i,j),sat%ivza(i,j),sat%isfc(i,j),i,j,a,sat%chan_flag, wflag)

enddo
enddo
!=== 2.3 

!print*,'RTM test', rtm(100,100)%d(20)%bt_clr38,rtm(100,100)%d(20)%bt_clr40, &
!                   rtm(100,100)%d(20)%bt_clr73,rtm(100,100)%d(20)%bt_clr86, &
!                   rtm(100,100)%d(20)%bt_clr11,rtm(100,100)%d(20)%bt_clr12


! 3. end subroutine   
end subroutine ir_rtm_driver
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~ subroutine 4: allocate rtm program ~~~~~~~~~~~~~~~~
subroutine allocate_rtm(nx,ny)

!-----------------------------------------------------------------------
! !F90 allocate rtm
!
! !Description:
!    The main program for allocating memory for nlon x nlat rtm structure.
!
! !Input  parameters:
!    filename        = L1b hdf file's name
!
! !Output parameters:
!    none
!-----------------------------------------------------------------------

IMPLICIT NONE

  integer(kind=4), intent(in) :: nx, ny
  integer:: astatus
  
  allocate(rtm(nx,ny),stat=astatus)
  
  if (astatus /= 0) then
     print *,"(a,'Not enough memory to allocate rtm_params structure.')"
     stop
  endif
  
! 3. end subroutine   
end subroutine allocate_rtm
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~ subroutine 5: deallocate rtm program ~~~~~~~~~~~~~~
subroutine deallocate_rtm()

!-----------------------------------------------------------------------
! !F90 deallocate_rtm
!
! !Description:
!    The main program for deallocating memory used for nlon x nlat rtm structure.
!
! !Input  parameters:
!    none
!
! !Output parameters:
!    none
!
!-----------------------------------------------------------------------
  
  integer:: astatus
  
  deallocate(rtm,stat=astatus)
  
  if (astatus /= 0) then
     print *,"(a,'Error deallocating rtm_params structure.')"
     stop
  endif
  
! 3. end subroutine   
end subroutine deallocate_rtm
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~ subroutine 6: allocate rtm cell program ~~~~~~~~~~~
subroutine allocate_rtm_cell(ilon,ilat,nvza)

!-----------------------------------------------------------------------
! !F90 allocate_rtm_cell
!
! !Description:
!    Thisprogram is to allocate memory for the RTM structure. 
!
! !Input  parameters:
!    ilon            = longitude position of satellite pixel in nwp data
!    ilat            = latitude position of satellite pixel in nwp data
!    nvza            = number of satellite zenith angle
!
! !Output parameters:
!    none
!
!-----------------------------------------------------------------------

  integer(kind=4), intent(in) :: ilon, ilat
  integer(kind=4), intent(in) :: nvza

  integer:: astatus
  !print*,ilon,ilat,nvza
  allocate(rtm(ilon, ilat)%d(nvza),stat=astatus)
  !print*,'dd',rtm(ilon, ilat)%d(:)%flag
  !print*,'astatus',astatus
  if (astatus /= 0) then
     print *,"(a,'Not enough memory to allocate rtm_params structure.')"
     stop
  endif
  !print*,rtm(ilon,ilat)%flag
  rtm(ilon,ilat)%flag = 1
    !print*,rtm(ilon,ilat)%flag
    !print*,rtm(ilon, ilat)%d(:)%flag

! 3. end subroutine   
end subroutine allocate_rtm_cell
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~ subroutine 7: deallocate rtm cell program ~~~~~~~~~~
subroutine deallocate_rtm_cell(ilon,ilat)

!-----------------------------------------------------------------------
! !F90 deallocate_rtm_cell
!
! !Description:
!    This program is to deallocate memory for the RTM structure. 
!
! !Input  parameters:
!    none
!
! !Output parameters:
!    none
!
!-----------------------------------------------------------------------

  integer(kind=4), intent(in) :: ilon, ilat

  integer:: astatus
  
  deallocate(rtm(ilon, ilat)%d,stat=astatus)

  if (astatus /= 0) then
     print *,"(a,'Error deallocating rtm_params structure.')"
     stop
  endif
  
  rtm(ilon,ilat)%flag = 0

! 3. end subroutine   
end subroutine deallocate_rtm_cell
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~ subroutine 8: calculate vza mid values ~~~~~~~~~~~~
subroutine vza_mid(bin, RTM_VZA_BINSIZE, vza)

!-----------------------------------------------------------------------
! !F90 vza_mid
!
! !Description:
!    This program is to calculate view zenith angle middle values 
!
! !Input  parameters:
!    bin
!    RTM_VZA_BINSIZE
!
! !Output parameters:
!    vza
!
!-----------------------------------------------------------------------

  integer(kind=4), intent(in)  :: bin
  real (kind=4), intent(in)    :: RTM_VZA_BINSIZE
  real (kind=4), intent(out)   :: vza
  real (kind=4) :: temp

  temp = (bin-1) * (RTM_VZA_BINSIZE) + (RTM_VZA_BINSIZE) / 2.
  if (temp > 1.0) temp = 1.0
  vza = acos(temp)/DTOR
  
! 3. end subroutine   
end subroutine vza_mid
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~ subroutine 9: determine surface level at pixel level
subroutine get_pixel_sfc_level(nlevels, zprof, zsfc, isfc, a)

!-----------------------------------------------------------------------
! !F90 get_pixel_sfc_level
!
! !Description:
!    This program is to determine surface level at pixel level using high
!    resolution surface elevation database. 
!
! !Input  parameters:
!
! !Output parameters
!
!-----------------------------------------------------------------------

  integer(kind=4), intent(in) :: nlevels
  real(kind=4), dimension(:), intent(in) :: zprof
  real(kind=4), intent(in) :: zsfc
  integer(kind=4), intent(out) :: isfc
  real(kind=4), intent(out) :: a
  
  integer:: k, il
  
  do k = nlevels/2, nlevels
     if (zprof(k) < zsfc) then
        isfc = k
        exit
     endif
  end do
  il = isfc - 1
  a = (zsfc - zprof(il))/(zprof(il+1) - zprof(il))

! 3. end subroutine   
end subroutine get_pixel_sfc_level
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~ SUBROUTINE 11:  compute clear sky atmospheric ~~~~~~
!                                    radiance profiles.
subroutine clear_radiance_prof(ichan, t_prof, tau_prof, rad_prof, cloud_prof)

!-----------------------------------------------------------------------
! !F90 clear_radiance_prof
!
! !Description:
!    This program is to compute clear sky atmospheric radiance profiles. 
! 
! !Input  parameters:
!    none
!
! !Output parameters:
!    none
!
!-----------------------------------------------------------------------

! 1. define variables
integer(kind=4), intent(in) :: ichan
real(kind=4), dimension(:), intent(in) :: t_prof, tau_prof
real(kind=4), dimension(:), intent(out) :: rad_prof, cloud_prof
integer(kind=4) :: nlev, ilev
  
!integer :: ilev, nlev
real (kind=4) :: dtrn, B1, B2
  
! 2. begin program    
nlev = size(t_prof,dim=1)
  
B1 = planck_rad_fast(ichan, t_prof(1))
rad_prof(1) = 0.0
cloud_prof(1) = B1*tau_prof(1)
  
do ilev=2, nlev
   B2 = planck_rad_fast(ichan, t_prof(ilev))
   !print*,'B2',ilev,B2
   dtrn = -(tau_prof(ilev) - tau_prof(ilev-1))
  ! print*,'dtrn',dtrn, tau_prof(ilev) , tau_prof(ilev-1), rad_prof(ilev-1)
   rad_prof(ilev) = rad_prof(ilev-1) + (B1+B2)/2.0 * dtrn
   !print*,'rad',ilev,rad_prof(ilev)
!   !Should subtract clear sky at the algorithm level since it may be a pixel resolution variable
   cloud_prof(ilev) = rad_prof(ilev) + B2*tau_prof(ilev)
   B1 = B2
end do

! 3. END SUBROUTINE   
end subroutine clear_radiance_prof
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~ subroutine 10:  calculate pixel-level clear-sky TOA ~
subroutine get_pixel_clear_rad(ilon, ilat, ivza, isfc, xpix, ypix, a, chflg, water_flag)

!-----------------------------------------------------------------------
! !F90 get_pixel_sfc_level
!
! !Description:
!    This program is to calculate pixel-level clear-sky TOA.  
!
! !Input  parameters:
!
! !Output parameters
!
!-----------------------------------------------------------------------

  integer(kind=4), intent(in) :: ilon, ilat, ivza, isfc, xpix, ypix
  integer(kind=1) water_flag
  real(kind=4), intent(in) :: a
  integer(kind=1), dimension(:), intent(in) :: chflg
  
  real (kind=4), dimension(:), pointer :: tau, rad, cloud_prof
  real (kind=4) :: tsfc, esfc
  !real (kind=4) :: lapse_rate
  integer:: il
  
  il = isfc - 1
  
  if (water_flag == 0) then
     tsfc = nwp%dat(ilon,ilat)%tsfc
  endif
  if (water_flag > 0) then
     tsfc = sat%sst(xpix, ypix)
     if (tsfc < 0.0 .or. tsfc >1000) then
         tsfc = nwp%dat(ilon,ilat)%tsfc
     endif
  endif

  if (tsfc < 0.0 .or. tsfc >1000) then
      tsfc = 273.15
  endif
  !print*,'tsfc',tsfc,water_flag,isfc ,il,ivza,a,ilon,ilat
  !if (sat%sfc_type(xpix,ypix) > 0) then
  !  lapse_rate = (nwp%dat(ilon,ilat)%tsfc - nwp%dat(ilon,ilat)%tlev(nwp%dat(ilon,ilat)%sfc_level))/ &
  !               (nwp%dat(ilon,ilat)%zsfc - nwp%dat(ilon,ilat)%zlev(nwp%dat(ilon,ilat)%sfc_level))
  ! 
  !  tsfc = nwp%dat(ilon,ilat)%tsfc + lapse_rate*(sat%zsfc(xpix,ypix) - nwp%dat(ilon,ilat)%zlev(isfc))
  !  !print*,sat%sfc_type(xpix,ypix),lapse_rate,tsfc,nwp%dat(ilon,ilat)%tsfc,sat%zsfc(xpix,ypix)
  !endif
  if (chflg(20) > 0) then 
    tau => rtm(ilon,ilat)%d(ivza)%trans_atm_clr38
    rad => rtm(ilon,ilat)%d(ivza)%rad_atm_clr38
    cloud_prof => rtm(ilon,ilat)%d(ivza)%cloud_prof38
    esfc = sat%sfc_emiss38(xpix,ypix)
    tau(isfc) = tau(il) + a * (tau(il+1) - tau(il))
    rad(isfc) = rad(il) + a * (rad(il+1) - rad(il))
    cloud_prof(isfc) = cloud_prof(il) + a * (cloud_prof(il+1) - cloud_prof(il))
    call clear_radiance_toa(1, 20, rad(isfc), tau(isfc), tsfc, esfc, &
                            sat%rad_clr38(xpix,ypix), sat%bt_clr38(xpix,ypix))
                         
    tau => null()
    rad => null()
    cloud_prof => null()
  endif

  if (chflg(21) > 0) then 
    tau => rtm(ilon,ilat)%d(ivza)%trans_atm_clr40
    rad => rtm(ilon,ilat)%d(ivza)%rad_atm_clr40
    cloud_prof => rtm(ilon,ilat)%d(ivza)%cloud_prof40
    esfc = sat%sfc_emiss40(xpix,ypix)
    tau(isfc) = tau(il) + a * (tau(il+1) - tau(il))
    rad(isfc) = rad(il) + a * (rad(il+1) - rad(il))
    cloud_prof(isfc) = cloud_prof(il) + a * (cloud_prof(il+1) - cloud_prof(il))
    
    call clear_radiance_toa(1, 21, rad(isfc), tau(isfc), tsfc, esfc, &
                            sat%rad_clr40(xpix,ypix), sat%bt_clr40(xpix,ypix))
    tau => null()
    rad => null()
    cloud_prof => null()
  endif
 
  if (chflg(22) > 0) then 
    tau => rtm(ilon,ilat)%d(ivza)%trans_atm_clr73
    rad => rtm(ilon,ilat)%d(ivza)%rad_atm_clr73
    cloud_prof => rtm(ilon,ilat)%d(ivza)%cloud_prof73
    esfc = sat%sfc_emiss73(xpix,ypix)
    tau(isfc) = tau(il) + a * (tau(il+1) - tau(il))
    rad(isfc) = rad(il) + a * (rad(il+1) - rad(il))
    cloud_prof(isfc) = cloud_prof(il) + a * (cloud_prof(il+1) - cloud_prof(il))
    
    call clear_radiance_toa(1, 22, rad(isfc), tau(isfc), tsfc, esfc, &
                            sat%rad_clr73(xpix,ypix), sat%bt_clr73(xpix,ypix))
    tau => null()
    rad => null()
    cloud_prof => null()
  endif

  if (chflg(23) > 0) then 
    tau => rtm(ilon,ilat)%d(ivza)%trans_atm_clr86
    rad => rtm(ilon,ilat)%d(ivza)%rad_atm_clr86
    cloud_prof => rtm(ilon,ilat)%d(ivza)%cloud_prof86
    esfc = sat%sfc_emiss86(xpix,ypix)
    tau(isfc) = tau(il) + a * (tau(il+1) - tau(il))
    rad(isfc) = rad(il) + a * (rad(il+1) - rad(il))
    cloud_prof(isfc) = cloud_prof(il) + a * (cloud_prof(il+1) - cloud_prof(il)) 
    
    call clear_radiance_toa(1, 23, rad(isfc), tau(isfc), tsfc, esfc, &
                            sat%rad_clr86(xpix,ypix), sat%bt_clr86(xpix,ypix))
    tau => null()
    rad => null()
    cloud_prof => null()
  endif
    
  if (chflg(24) > 0) then 
    tau => rtm(ilon,ilat)%d(ivza)%trans_atm_clr11
    rad => rtm(ilon,ilat)%d(ivza)%rad_atm_clr11
    cloud_prof => rtm(ilon,ilat)%d(ivza)%cloud_prof11
    esfc = sat%sfc_emiss11(xpix,ypix)
    tau(isfc) = tau(il) + a * (tau(il+1) - tau(il))
    rad(isfc) = rad(il) + a * (rad(il+1) - rad(il))
    cloud_prof(isfc) = cloud_prof(il) + a * (cloud_prof(il+1) - cloud_prof(il))
    
    call clear_radiance_toa(1, 24, rad(isfc), tau(isfc), tsfc, esfc, &
                            sat%rad_clr11(xpix,ypix), sat%bt_clr11(xpix,ypix))
  !  print*,'11.0',rad(isfc),tau(isfc),sat%rad_clr11(xpix,ypix), sat%bt_clr11(xpix,ypix)
    tau => null()
    rad => null()
    cloud_prof => null()
  endif
!print*,'25'
  if (chflg(25) > 0) then 
    tau => rtm(ilon,ilat)%d(ivza)%trans_atm_clr12
    rad => rtm(ilon,ilat)%d(ivza)%rad_atm_clr12
    cloud_prof => rtm(ilon,ilat)%d(ivza)%cloud_prof12
    esfc = sat%sfc_emiss12(xpix,ypix)
    tau(isfc) = tau(il) + a * (tau(il+1) - tau(il))
    rad(isfc) = rad(il) + a * (rad(il+1) - rad(il))
    cloud_prof(isfc) = cloud_prof(il) + a * (cloud_prof(il+1) - cloud_prof(il))  
    
    call clear_radiance_toa(1, 25, rad(isfc), tau(isfc), tsfc, esfc, &
                            sat%rad_clr12(xpix,ypix), sat%bt_clr12(xpix,ypix))
    tau => null()
    rad => null()
    cloud_prof => null()
  endif
  
!  print*,'2777'
! 3. end subroutine   
end subroutine get_pixel_clear_rad
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!~~~~~~~~~~~~~~~~~~~ subroutine 12:  compute clear sky brightness ~~~~~~~
subroutine clear_radiance_toa(option, ichan, rad_atm, tau_atm, tsfc, esfc, &
                              rad_clr, bt_clr)

!-----------------------------------------------------------------------
! !F90 clear_radiance_toa
!
! !Description:
!    This program is to compute clear sky brightness temperature and radiance.
! 
! !Input  parameters:
!    none
!
! !Output parameters:
!    none
!
!-----------------------------------------------------------------------

! 1. define variables
  integer(kind=4), intent(in) :: option, ichan
  real(kind=4), intent(in) :: rad_atm, tau_atm, tsfc, esfc
  real(kind=4), intent(out) :: rad_clr, bt_clr
  real(kind=4) esfc2
! 2. begin program  
  if (esfc<=0.) then
    esfc2 = 0.99
  else
    esfc2 = esfc
  endif
 
  if (option == 1) then
    !rad_clr = rad_atm + esfc2*planck_rad_fast(ichan, tsfc)*tau_atm
    rad_clr = rad_atm + esfc2*fylat_planck_tbb2rad(tsfc, ichan, 1)*tau_atm
!    print*,'rr',rad_clr,rad_atm,tau_atm,ichan,tsfc
    !bt_clr =  planck_temp_fast(ichan, rad_clr)
    !CALL  msg_planck_temp(ichan, rad_clr, bt_clr)
    ! CALL call_planck_temp_func(sc_ind, ichan, rad_clr, bt_clr)
    !CALL call_planck_temp_func(ichan, rad_clr, bt_clr)
    bt_clr = fylat_planck_rad2tbb(rad_clr, ichan, 1)
!    print*,'tbb',bt_clr
!    print*,'-----------'
  else
    rad_clr = 0.0
    bt_clr = 0.0
  endif

! 3. end subroutine   
end subroutine clear_radiance_toa
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~ subroutine run_3: run plod fy4 ~~~~~~~~~~~~~~~~~~~~
subroutine run_plod_mersiII_modis(ilon,ilat,ivza,nlevels,option)

!-----------------------------------------------------------------------
! !F90 deal_rtm_begin
!
! !Description:
!    The begin program for deal with rtm data.[allocate variables]
!
! !Input  parameters:
!    none
!
! !Output parameters:
!    none
!
!-----------------------------------------------------------------------
  integer (kind=4), intent(in) :: ilon, ilat
  integer (kind=4), intent(in) :: ivza, nlevels, option
  
  integer (kind=4) :: pnum, isfc
  real (kind=4) :: satzen, tsfc, esfc
  real (kind=4), dimension(:), pointer :: tlev,wlev,o3lev,tau,rad,cloud_prof,ax
  !real (kind=4), dimension(101) :: tlev,wlev,o3lev
  !real (kind=4), dimension(:), pointer :: tau,rad,cloud_prof
  integer :: astatus, il, status
  real (kind=4) :: a

  CHARACTER(LEN=6) :: Craft
  
  real(kind=4), parameter :: Co2_Mix = 380 !ppm
  
  Craft = 'AQUA'
  
  esfc = 0.99
 ! tlev = nwp%dat(ilon,ilat)%tlev(:)
 ! wlev = nwp%dat(ilon,ilat)%wlev(:)
 ! o3lev = nwp%dat(ilon,ilat)%o3lev(:)
  tlev => nwp%dat(ilon,ilat)%tlev(:)
  wlev => nwp%dat(ilon,ilat)%wlev(:)
  o3lev => nwp%dat(ilon,ilat)%o3lev(:)
  tsfc = nwp%dat(ilon,ilat)%tsfc
  satzen = rtm(ilon,ilat)%d(ivza)%satzen
  !pnum = int(scinfo%pnum,kind=4)
  if (option == 1) then
    isfc = nwp%dat(ilon,ilat)%sfc_level
  else
    isfc = nlevels
  endif
  a = nwp%dat(ilon,ilat)%a
  il = isfc - 1
  
  allocate(rtm(ilon,ilat)%d(ivza)%rtm_util(nlevels), stat=astatus)
  if (astatus /= 0) then
    print *,"(a,'Not enough memory to allocate rtm utility pointer.')"
    stop
  endif
  
  !------------------------------------------------------------------
  ! Channel 3.8um 
  !------------------------------------------------------------------ 
  tau => null()
  rad => null()
  cloud_prof => null()
  if (sat%chan_flag(20) > 0) then 
    allocate(rtm(ilon,ilat)%d(ivza)%trans_atm_clr38(nlevels),&
             rtm(ilon,ilat)%d(ivza)%rad_atm_clr38(nlevels),  &
             rtm(ilon,ilat)%d(ivza)%cloud_prof38(nlevels), stat=astatus)
    if (astatus /= 0) then
       print *,"(a,'Not enough memory to allocate channel 3.8um rtm profile.')"
       stop
    endif
    tau => rtm(ilon,ilat)%d(ivza)%trans_atm_clr38
    rad => rtm(ilon,ilat)%d(ivza)%rad_atm_clr38
    cloud_prof => rtm(ilon,ilat)%d(ivza)%cloud_prof38
   !tau = rtm(ilon,ilat)%d(ivza)%trans_atm_clr7
   !rad = rtm(ilon,ilat)%d(ivza)%rad_atm_clr7
   !cloud_prof = rtm(ilon,ilat)%d(ivza)%cloud_prof7
    call tranvmodisd101(tlev,                                    &
                        wlev,                                    &
                        o3lev,                                   &
                        satzen,                                  &
                        Co2_Mix,                                 &
                        Craft,                                   &
                        20,                                      &
                        0,                                       &
                        tau,                                     &
                        status)
    tau(isfc) = tau(il) + a * (tau(il+1) - tau(il))
    CALL clear_radiance_prof(20, tlev(1:isfc), tau(1:isfc), rad(1:isfc), cloud_prof(1:isfc))
    CALL clear_radiance_toa(option, 20, rad(isfc), tau(isfc), tsfc, esfc, &
                            rtm(ilon,ilat)%d(ivza)%rad_clr38, rtm(ilon,ilat)%d(ivza)%bt_clr38)
    tau => null()
    rad => null()
    cloud_prof => null()
   ! tau = 0.
   ! rad = 0.
   ! cloud_prof = 0.
  endif
  
  !------------------------------------------------------------------
  ! Channel 4.05 
  !------------------------------------------------------------------
  if (sat%chan_flag(21) > 0) then 
    allocate(rtm(ilon,ilat)%d(ivza)%trans_atm_clr40(nlevels),&
             rtm(ilon,ilat)%d(ivza)%rad_atm_clr40(nlevels),  &
             rtm(ilon,ilat)%d(ivza)%cloud_prof40(nlevels), stat=astatus)
    if (astatus /= 0) then
       print *,"(a,'Not enough memory to allocate channel 4.05um rtm profile.')"
       stop
    endif
    tau => rtm(ilon,ilat)%d(ivza)%trans_atm_clr40
    rad => rtm(ilon,ilat)%d(ivza)%rad_atm_clr40
    cloud_prof => rtm(ilon,ilat)%d(ivza)%cloud_prof40
    call tranvmodisd101(tlev,                                    &
                        wlev,                                    &
                        o3lev,                                   &
                        satzen,                                  &
                        Co2_Mix,                                 &
                        Craft,                                   &
                        23,                                      &
                        0,                                       &
                        tau,                                     &
                        status)
    tau(isfc) = tau(il) + a * (tau(il+1) - tau(il))
    CALL clear_radiance_prof(21, tlev(1:isfc), tau(1:isfc), rad(1:isfc), cloud_prof(1:isfc))
    CALL clear_radiance_toa(option, 21, rad(isfc), tau(isfc), tsfc, esfc, &
                            rtm(ilon,ilat)%d(ivza)%rad_clr40, rtm(ilon,ilat)%d(ivza)%bt_clr40)
    tau => null()
    rad => null()
    cloud_prof => null()
   ! tau = 0.
   ! rad = 0.
   ! cloud_prof = 0.
  endif

  !------------------------------------------------------------------
  ! Channel 7.3um
  !------------------------------------------------------------------
  if (sat%chan_flag(22) > 0) then
    allocate(rtm(ilon,ilat)%d(ivza)%trans_atm_clr73(nlevels),&
             rtm(ilon,ilat)%d(ivza)%rad_atm_clr73(nlevels),  &
             rtm(ilon,ilat)%d(ivza)%cloud_prof73(nlevels), stat=astatus)
    if (astatus /= 0) then
       print *,"(a,'Not enough memory to allocate channel 7.3um rtm profile.')"
       stop
    endif
    tau => rtm(ilon,ilat)%d(ivza)%trans_atm_clr73
    rad => rtm(ilon,ilat)%d(ivza)%rad_atm_clr73
    cloud_prof => rtm(ilon,ilat)%d(ivza)%cloud_prof73
    call tranvmodisd101(tlev,                                    &
                        wlev,                                    &
                        o3lev,                                   &
                        satzen,                                  &
                        Co2_Mix,                                 &
                        Craft,                                   &
                        28,                                      &
                        0,                                       &
                        tau,                                     &
                        status)
    tau(isfc) = tau(il) + a * (tau(il+1) - tau(il))
    CALL clear_radiance_prof(22, tlev(1:isfc), tau(1:isfc), rad(1:isfc), cloud_prof(1:isfc))
    CALL clear_radiance_toa(option, 22, rad(isfc), tau(isfc), tsfc, esfc, &
                            rtm(ilon,ilat)%d(ivza)%rad_clr73, rtm(ilon,ilat)%d(ivza)%bt_clr73)
    tau => null()
    rad => null()
    cloud_prof => null()
  endif
   
  !------------------------------------------------------------------
  ! Channel 8.6um
  !------------------------------------------------------------------
  if (sat%chan_flag(23) > 0) then 
    allocate(rtm(ilon,ilat)%d(ivza)%trans_atm_clr86(nlevels),&
             rtm(ilon,ilat)%d(ivza)%rad_atm_clr86(nlevels),  &
             rtm(ilon,ilat)%d(ivza)%cloud_prof86(nlevels), stat=astatus)
    if (astatus /= 0) then
       print *,"(a,'Not enough memory to allocate channel 8.6um rtm profile.')"
       stop
    endif
    tau => rtm(ilon,ilat)%d(ivza)%trans_atm_clr86
    rad => rtm(ilon,ilat)%d(ivza)%rad_atm_clr86
    cloud_prof => rtm(ilon,ilat)%d(ivza)%cloud_prof86
    call tranvmodisd101(tlev,                                    &
                        wlev,                                    &
                        o3lev,                                   &
                        satzen,                                  &
                        Co2_Mix,                                 &
                        Craft,                                   &
                        29,                                      &
                        0,                                       &
                        tau,                                     &
                        status)
    tau(isfc) = tau(il) + a * (tau(il+1) - tau(il))
    CALL clear_radiance_prof(23, tlev(1:isfc), tau(1:isfc), rad(1:isfc), cloud_prof(1:isfc))
    CALL clear_radiance_toa(option, 23, rad(isfc), tau(isfc), tsfc, esfc, &
                            rtm(ilon,ilat)%d(ivza)%rad_clr86, rtm(ilon,ilat)%d(ivza)%bt_clr86)
    tau => null()
    rad => null()
    cloud_prof => null()
  endif
  
  !------------------------------------------------------------------
  ! Channel 11. 
  !------------------------------------------------------------------
  
  if (sat%chan_flag(24) > 0) then 
    allocate(rtm(ilon,ilat)%d(ivza)%trans_atm_clr11(nlevels),&
             rtm(ilon,ilat)%d(ivza)%rad_atm_clr11(nlevels),&
             rtm(ilon,ilat)%d(ivza)%cloud_prof11(nlevels), stat=astatus)
    if (astatus /= 0) then
       print *,"(a,'Not enough memory to allocate channel 11um rtm profile.')"
       stop
    endif
    tau => rtm(ilon,ilat)%d(ivza)%trans_atm_clr11
    rad => rtm(ilon,ilat)%d(ivza)%rad_atm_clr11
    cloud_prof => rtm(ilon,ilat)%d(ivza)%cloud_prof11
    call tranvmodisd101(tlev,                                    &
                        wlev,                                    &
                        o3lev,                                   &
                        satzen,                                  &
                        Co2_Mix,                                 &
                        Craft,                                   &
                        31,                                      &
                        0,                                       &
                        tau,                                     &
                        status)
    tau(isfc) = tau(il) + a * (tau(il+1) - tau(il))
    CALL clear_radiance_prof(24, tlev(1:isfc), tau(1:isfc), rad(1:isfc), cloud_prof(1:isfc))
    CALL clear_radiance_toa(option, 24, rad(isfc), tau(isfc), tsfc, esfc, &
                            rtm(ilon,ilat)%d(ivza)%rad_clr11, rtm(ilon,ilat)%d(ivza)%bt_clr11)
!print*,'222-11',rad(isfc),tau(isfc), tsfc, esfc, isfc,rtm(ilon,ilat)%d(ivza)%rad_clr11, rtm(ilon,ilat)%d(ivza)%bt_clr11
 
    tau => null()
    rad => null()
    cloud_prof => null()
  endif
  
  !------------------------------------------------------------------
  ! Channel 12. 
  !------------------------------------------------------------------
  
  if (sat%chan_flag(25) > 0) then 
    allocate(rtm(ilon,ilat)%d(ivza)%trans_atm_clr12(nlevels),&
             rtm(ilon,ilat)%d(ivza)%rad_atm_clr12(nlevels),&
             rtm(ilon,ilat)%d(ivza)%cloud_prof12(nlevels), stat=astatus)
    if (astatus /= 0) then
       print *,"(a,'Not enough memory to allocate channel 12um rtm profile.')"
       stop
    endif
    tau => rtm(ilon,ilat)%d(ivza)%trans_atm_clr12
    rad => rtm(ilon,ilat)%d(ivza)%rad_atm_clr12
    cloud_prof => rtm(ilon,ilat)%d(ivza)%cloud_prof12
    call tranvmodisd101(tlev,                                    &
                        wlev,                                    &
                        o3lev,                                   &
                        satzen,                                  &
                        Co2_Mix,                                 &
                        Craft,                                   &
                        32,                                      &
                        0,                                       &
                        tau,                                     &
                        status)
    tau(isfc) = tau(il) + a * (tau(il+1) - tau(il))
    CALL clear_radiance_prof(25, tlev(1:isfc), tau(1:isfc), rad(1:isfc), cloud_prof(1:isfc))
    CALL clear_radiance_toa(option, 25, rad(isfc), tau(isfc), tsfc, esfc, &
                            rtm(ilon,ilat)%d(ivza)%rad_clr12, rtm(ilon,ilat)%d(ivza)%bt_clr12)
!print*,'333-12',rad(isfc),tau(isfc), tsfc, esfc, isfc,rtm(ilon,ilat)%d(ivza)%rad_clr12, rtm(ilon,ilat)%d(ivza)%bt_clr12

    tau => null()
    rad => null()
    cloud_prof => null()
  endif

  rtm(ilon,ilat)%d(ivza)%flag = 1

! 3. end subroutine  
end subroutine run_plod_mersiII_modis
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~ subroutine run_3: run plod fy3d mersi ii ~~~~~~~~~~
subroutine run_plod_fy3d_mersi_ii(ilon,ilat,ivza,nlevels,option)

!-----------------------------------------------------------------------
! !F90 deal_rtm_begin
!
! !Description:
!    The begin program for deal with rtm data.[allocate variables]
!
! !Input  parameters:
!    none
!
! !Output parameters:
!    none
!
!-----------------------------------------------------------------------
  integer (kind=4), intent(in) :: ilon, ilat
  integer (kind=4), intent(in) :: ivza, nlevels, option
  
  integer (kind=4) :: pnum, isfc
  real (kind=4) :: satzen, tsfc, esfc
  real (kind=4), dimension(:), pointer :: tlev,wlev,o3lev,tau,rad,cloud_prof,ax
  !real (kind=4), dimension(101) :: tlev,wlev,o3lev
  !real (kind=4), dimension(:), pointer :: tau,rad,cloud_prof
  integer :: astatus, il, status
  real (kind=4) :: a

  CHARACTER(LEN=6) :: Craft
  
  real(kind=4), parameter :: Co2_Mix = 380 !ppm
  
  !Craft = 'AQUA'
  
  esfc = 0.99
 ! tlev = nwp%dat(ilon,ilat)%tlev(:)
 ! wlev = nwp%dat(ilon,ilat)%wlev(:)
 ! o3lev = nwp%dat(ilon,ilat)%o3lev(:)
  tlev => nwp%dat(ilon,ilat)%tlev(:)
  wlev => nwp%dat(ilon,ilat)%wlev(:)
  o3lev => nwp%dat(ilon,ilat)%o3lev(:)
  tsfc = nwp%dat(ilon,ilat)%tsfc
  satzen = rtm(ilon,ilat)%d(ivza)%satzen
  !pnum = int(scinfo%pnum,kind=4)
  if (option == 1) then
    isfc = nwp%dat(ilon,ilat)%sfc_level
  else
    isfc = nlevels
  endif
  a = nwp%dat(ilon,ilat)%a
  il = isfc - 1
  
  allocate(rtm(ilon,ilat)%d(ivza)%rtm_util(nlevels), stat=astatus)
  if (astatus /= 0) then
     print *,"(a,'Not enough memory to allocate rtm utility pointer.')"
     stop
  endif
  
  !------------------------------------------------------------------
  ! Channel 3.8um 
  !------------------------------------------------------------------ 
  tau => null()
  rad => null()
  cloud_prof => null()
  if (sat%chan_flag(20) > 0) then 
    allocate(rtm(ilon,ilat)%d(ivza)%trans_atm_clr38(nlevels),&
             rtm(ilon,ilat)%d(ivza)%rad_atm_clr38(nlevels),  &
             rtm(ilon,ilat)%d(ivza)%cloud_prof38(nlevels), stat=astatus)
    if (astatus /= 0) then
       print *,"(a,'Not enough memory to allocate channel 3.8um rtm profile.')"
       stop
    endif
    tau => rtm(ilon,ilat)%d(ivza)%trans_atm_clr38
    rad => rtm(ilon,ilat)%d(ivza)%rad_atm_clr38
    cloud_prof => rtm(ilon,ilat)%d(ivza)%cloud_prof38
   !tau = rtm(ilon,ilat)%d(ivza)%trans_atm_clr7
   !rad = rtm(ilon,ilat)%d(ivza)%rad_atm_clr7
   !cloud_prof = rtm(ilon,ilat)%d(ivza)%cloud_prof7
    call fy3dtrn101  (  tlev,                                    &
                        wlev,                                    &
                        o3lev,                                   &
                        satzen,                                  &
!                        Co2_Mix,                                &
!                        Craft,                                  &
                         2,                                      &
!                        0,                                      &
                        tau,                                     &
                        status)
                        
    tau(isfc) = tau(il) + a * (tau(il+1) - tau(il))
    CALL clear_radiance_prof(20, tlev(1:isfc), tau(1:isfc), rad(1:isfc), cloud_prof(1:isfc))
    CALL clear_radiance_toa(option, 20, rad(isfc), tau(isfc), tsfc, esfc, &
                            rtm(ilon,ilat)%d(ivza)%rad_clr38, rtm(ilon,ilat)%d(ivza)%bt_clr38)
!print*,'333-10',rad(isfc),tau(isfc), tsfc, esfc, isfc,rtm(ilon,ilat)%d(ivza)%rad_clr38, rtm(ilon,ilat)%d(ivza)%bt_clr38

    tau => null()
    rad => null()
    cloud_prof => null()
   ! tau = 0.
   ! rad = 0.
   ! cloud_prof = 0.
  endif
  
  !------------------------------------------------------------------
  ! Channel 4.05 
  !------------------------------------------------------------------
  if (sat%chan_flag(21) > 0) then 
    allocate(rtm(ilon,ilat)%d(ivza)%trans_atm_clr40(nlevels),&
             rtm(ilon,ilat)%d(ivza)%rad_atm_clr40(nlevels),  &
             rtm(ilon,ilat)%d(ivza)%cloud_prof40(nlevels), stat=astatus)
    if (astatus /= 0) then
       print *,"(a,'Not enough memory to allocate channel 4.05um rtm profile.')"
       stop
    endif
    tau => rtm(ilon,ilat)%d(ivza)%trans_atm_clr40
    rad => rtm(ilon,ilat)%d(ivza)%rad_atm_clr40
    cloud_prof => rtm(ilon,ilat)%d(ivza)%cloud_prof40
    call fy3dtrn101  (  tlev,                                    &
                        wlev,                                    &
                        o3lev,                                   &
                        satzen,                                  &
!                        Co2_Mix,                                &
!                        Craft,                                  &
                         3,                                      &
!                        0,                                      &
                        tau,                                     &
                        status)
    tau(isfc) = tau(il) + a * (tau(il+1) - tau(il))
    CALL clear_radiance_prof(21, tlev(1:isfc), tau(1:isfc), rad(1:isfc), cloud_prof(1:isfc))
    CALL clear_radiance_toa(option, 21, rad(isfc), tau(isfc), tsfc, esfc, &
                            rtm(ilon,ilat)%d(ivza)%rad_clr40, rtm(ilon,ilat)%d(ivza)%bt_clr40)
    tau => null()
    rad => null()
    cloud_prof => null()
   ! tau = 0.
   ! rad = 0.
   ! cloud_prof = 0.
  endif

  !------------------------------------------------------------------
  ! Channel 7.3um
  !------------------------------------------------------------------
  if (sat%chan_flag(22) > 0) then
    allocate(rtm(ilon,ilat)%d(ivza)%trans_atm_clr73(nlevels),&
             rtm(ilon,ilat)%d(ivza)%rad_atm_clr73(nlevels),  &
             rtm(ilon,ilat)%d(ivza)%cloud_prof73(nlevels), stat=astatus)
    if (astatus /= 0) then
       print *,"(a,'Not enough memory to allocate channel 7.3um rtm profile.')"
       stop
    endif
    tau => rtm(ilon,ilat)%d(ivza)%trans_atm_clr73
    rad => rtm(ilon,ilat)%d(ivza)%rad_atm_clr73
    cloud_prof => rtm(ilon,ilat)%d(ivza)%cloud_prof73
    call fy3dtrn101  (  tlev,                                    &
                        wlev,                                    &
                        o3lev,                                   &
                        satzen,                                  &
!                        Co2_Mix,                                &
!                        Craft,                                  &
                         4,                                      &
!                        0,                                      &
                        tau,                                     &
                        status)
    tau(isfc) = tau(il) + a * (tau(il+1) - tau(il))
    CALL clear_radiance_prof(22, tlev(1:isfc), tau(1:isfc), rad(1:isfc), cloud_prof(1:isfc))
    CALL clear_radiance_toa(option, 22, rad(isfc), tau(isfc), tsfc, esfc, &
                            rtm(ilon,ilat)%d(ivza)%rad_clr73, rtm(ilon,ilat)%d(ivza)%bt_clr73)
    tau => null()
    rad => null()
    cloud_prof => null()
  endif
   
  !------------------------------------------------------------------
  ! Channel 8.6um
  !------------------------------------------------------------------
  if (sat%chan_flag(23) > 0) then 
    allocate(rtm(ilon,ilat)%d(ivza)%trans_atm_clr86(nlevels),&
             rtm(ilon,ilat)%d(ivza)%rad_atm_clr86(nlevels),  &
             rtm(ilon,ilat)%d(ivza)%cloud_prof86(nlevels), stat=astatus)
    if (astatus /= 0) then
       print *,"(a,'Not enough memory to allocate channel 8.6um rtm profile.')"
       stop
    endif
    tau => rtm(ilon,ilat)%d(ivza)%trans_atm_clr86
    rad => rtm(ilon,ilat)%d(ivza)%rad_atm_clr86
    cloud_prof => rtm(ilon,ilat)%d(ivza)%cloud_prof86
    call fy3dtrn101  (  tlev,                                    &
                        wlev,                                    &
                        o3lev,                                   &
                        satzen,                                  &
!                        Co2_Mix,                                &
!                        Craft,                                  &
                          5,                                     &
!                        0,                                      &
                        tau,                                     &
                        status)
    tau(isfc) = tau(il) + a * (tau(il+1) - tau(il))
    CALL clear_radiance_prof(23, tlev(1:isfc), tau(1:isfc), rad(1:isfc), cloud_prof(1:isfc))
    CALL clear_radiance_toa(option, 23, rad(isfc), tau(isfc), tsfc, esfc, &
                            rtm(ilon,ilat)%d(ivza)%rad_clr86, rtm(ilon,ilat)%d(ivza)%bt_clr86)
    tau => null()
    rad => null()
    cloud_prof => null()
  endif
  
  !------------------------------------------------------------------
  ! Channel 11. 
  !------------------------------------------------------------------
  if (sat%chan_flag(24) > 0) then 
    allocate(rtm(ilon,ilat)%d(ivza)%trans_atm_clr11(nlevels),&
             rtm(ilon,ilat)%d(ivza)%rad_atm_clr11(nlevels),&
             rtm(ilon,ilat)%d(ivza)%cloud_prof11(nlevels), stat=astatus)
    if (astatus /= 0) then
       print *,"(a,'Not enough memory to allocate channel 11um rtm profile.')"
       stop
    endif
    tau => rtm(ilon,ilat)%d(ivza)%trans_atm_clr11
    rad => rtm(ilon,ilat)%d(ivza)%rad_atm_clr11
    cloud_prof => rtm(ilon,ilat)%d(ivza)%cloud_prof11
    call fy3dtrn101  (  tlev,                                    &
                        wlev,                                    &
                        o3lev,                                   &
                        satzen,                                  &
!                        Co2_Mix,                                &
!                        Craft,                                  &
                         6,                                      &
!                        0,                                      &
                        tau,                                     &
                        status)
    tau(isfc) = tau(il) + a * (tau(il+1) - tau(il))
    CALL clear_radiance_prof(24, tlev(1:isfc), tau(1:isfc), rad(1:isfc), cloud_prof(1:isfc))
    CALL clear_radiance_toa(option, 24, rad(isfc), tau(isfc), tsfc, esfc, &
                            rtm(ilon,ilat)%d(ivza)%rad_clr11, rtm(ilon,ilat)%d(ivza)%bt_clr11)
!print*,'222-11',rad(isfc),tau(isfc), tsfc, esfc, isfc,rtm(ilon,ilat)%d(ivza)%rad_clr11, rtm(ilon,ilat)%d(ivza)%bt_clr11
 
    tau => null()
    rad => null()
    cloud_prof => null()
  endif
  
  !------------------------------------------------------------------
  ! Channel 12. 
  !------------------------------------------------------------------
  
  if (sat%chan_flag(25) > 0) then 
    allocate(rtm(ilon,ilat)%d(ivza)%trans_atm_clr12(nlevels),&
             rtm(ilon,ilat)%d(ivza)%rad_atm_clr12(nlevels),&
             rtm(ilon,ilat)%d(ivza)%cloud_prof12(nlevels), stat=astatus)
    if (astatus /= 0) then
       print *,"(a,'Not enough memory to allocate channel 12um rtm profile.')"
       stop
    endif
    tau => rtm(ilon,ilat)%d(ivza)%trans_atm_clr12
    rad => rtm(ilon,ilat)%d(ivza)%rad_atm_clr12
    cloud_prof => rtm(ilon,ilat)%d(ivza)%cloud_prof12
    call fy3dtrn101  (  tlev,                                    &
                        wlev,                                    &
                        o3lev,                                   &
                        satzen,                                  &
!                        Co2_Mix,                                &
!                        Craft,                                  &
                         7,                                      &
!                        0,                                      &
                        tau,                                     &
                        status)
    tau(isfc) = tau(il) + a * (tau(il+1) - tau(il))
    CALL clear_radiance_prof(25, tlev(1:isfc), tau(1:isfc), rad(1:isfc), cloud_prof(1:isfc))
    CALL clear_radiance_toa(option, 25, rad(isfc), tau(isfc), tsfc, esfc, &
                            rtm(ilon,ilat)%d(ivza)%rad_clr12, rtm(ilon,ilat)%d(ivza)%bt_clr12)
!print*,'333-12',rad(isfc),tau(isfc), tsfc, esfc, isfc,rtm(ilon,ilat)%d(ivza)%rad_clr12, rtm(ilon,ilat)%d(ivza)%bt_clr12

    tau => null()
    rad => null()
    cloud_prof => null()
  endif
!print*,'tttt',ilon,ilat,ivza 562         144          38
if (ilon == 562 .and.  ilat == 144 .and. ivza==38) then
print*,'RTM test', rtm(ilon,ilat)%d(ivza)%bt_clr38,rtm(ilon,ilat)%d(ivza)%bt_clr40, &
                   rtm(ilon,ilat)%d(ivza)%bt_clr73,rtm(ilon,ilat)%d(ivza)%bt_clr86, &
                   rtm(ilon,ilat)%d(ivza)%bt_clr11,rtm(ilon,ilat)%d(ivza)%bt_clr12, &
                   rtm(ilon,ilat)%d(ivza)%rad_clr38,rtm(ilon,ilat)%d(ivza)%rad_clr40,  &
                   rtm(ilon,ilat)%d(ivza)%rad_clr73,rtm(ilon,ilat)%d(ivza)%rad_clr86,  &
                   rtm(ilon,ilat)%d(ivza)%rad_clr11,rtm(ilon,ilat)%d(ivza)%rad_clr12
endif
  rtm(ilon,ilat)%d(ivza)%flag = 1

! 3. end subroutine  
end subroutine run_plod_fy3d_mersi_ii
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
end module rtm_utils_module
