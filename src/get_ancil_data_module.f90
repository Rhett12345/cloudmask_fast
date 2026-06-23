module get_ancil_data_module


!C-----------------------------------------------------------------------
!C !F90                                                                  
!C
!C !description: 
!C    This module is to read ancillary data.
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
!C !end
!C----------------------------------------------------------------------

use names_module
use constant
use message_module
use numerical
use data_arrays_module
use hdf4
use hdf5

implicit none

! sfc emiss
real(kind=4), parameter, dimension(16), public :: emiss_water_sfc = (/ &
                                                0.000, 0.000, 0.000, 0.000, 0.000, 0.000, &
                                                0.978, 0.979, 0.979, 0.979, 0.983, 0.988, 0.990, 0.993, 0.986, 0.968/)
    
integer(kind=4), parameter, private :: num_lat_emiss = 3600
integer(kind=4), parameter, private :: num_lon_emiss = 7200
real(kind=4),    parameter, private :: first_lat_emiss = 89.9750, last_lat_emiss = -89.9750
real(kind=4),    parameter, private :: first_lon_emiss = -179.975, last_lon_emiss = 179.975
real(kind=4),    parameter, private :: del_lat_emiss = 0.05
real(kind=4),    parameter, private :: del_lon_emiss = 0.05

integer(kind=4), parameter, private :: num_lat_smk = 720
integer(kind=4), parameter, private :: num_lon_smk = 1440
integer(kind=4), parameter, private :: num_lat_nise = 721
integer(kind=4), parameter, private :: num_lon_nise = 721
real(kind=4), parameter, private    :: first_lat_smk = 89.875, last_lat_smk = -89.875
real(kind=4), parameter, private    :: first_lon_smk = -179.875, last_lon_smk = 179.875
real(kind=4), parameter, private    :: del_lat_smk = 0.25
real(kind=4), parameter, private    :: del_lon_smk = 0.25
integer(kind=4), parameter, private :: num_lon_oisst = 1440
integer(kind=4), parameter, private :: num_lat_oisst = 720

! sfc alb
integer(kind=4),parameter, dimension(23) :: Alb_day = (/ &
       1, 17, 33, 49, 65, 81, 97, 113, 129, 145, 161, 177, 193,  &
       209, 225, 241, 257, 273, 289, 305, 321, 337, 353/)
       
integer(kind=4), parameter, private :: num_lat_alb = 2700
integer(kind=4), parameter, private :: num_lon_alb = 5400
real(kind=4), parameter, private   :: first_lat_alb = 89.9600, last_lat_alb = -89.9600
real(kind=4), parameter, private   :: first_lon_alb = -179.9600, last_lon_alb = 179.9600
real(kind=4), parameter, private   :: del_lat_alb = 0.06666
real(kind=4), parameter, private   :: del_lon_alb = 0.06666

! 1KM ECOSYSTEM 
! --- Fortran OPEN 
!integer(kind=4), parameter      :: ECO_1KM_UNIT = 10
!character(len= 3), parameter	:: ECO_1KM_STATUS = 'OLD'
!character(len=11), parameter	:: ECO_1KM_FORM = 'UNFORMATTED'
!character(len= 6),parameter		:: ECO_1KM_ACCESS = 'DIRECT'
!integer(kind=4), parameter      :: ECO_1KM_RECL = 512



contains
!+++++++++++++++++++ step 2: subroutine s +++++++++++++++++++++++++++++++
!~~~~~~~~~~~~~~~~ subroutine   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine  allocate_fylat_ancil_data()

integer :: error

allocate (sat%snow_mask(sat%nElem,sat%nLine),   & 
          sat%eco(sat%nElem,sat%nLine),         &      
          sat%sst(sat%nElem,sat%nLine),         &         	         
          stat=error )
            
if (error /= 0) then
    print *,"(ERROR: 'Not enough memory to allocate fylat ancillary data arrays.')"
    stop
endif

end subroutine  allocate_fylat_ancil_data
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~ subroutine   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine  deallocate_fylat_ancil_data()

integer :: error

deallocate (sat%snow_mask,     & 
            sat%eco,           &
            sat%sst,           &
            stat=error)
            
if (error /= 0) then
    print *,"(ERROR: 'Not enough memory to deallocate fylat ancillary data arrays.')"
    stop
endif

call allocate_emiss_arrays()

end subroutine  deallocate_fylat_ancil_data
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~ subroutine  :  read_sfc_snow_ice_mask ~~~~~~~~~~~~~~
subroutine  read_sfc_snow_ice_mask(smonth)

!-----------------------------------------------------------------------
! !F90 read_snow_ice_mask
!
! !description:
!    This is a main program for reading surface snow mask file.
!
! !Input  parameters:
!    smonth             =  month for corresponding day data (4)
!          
! !Output parameters:
!    none
!
!-----------------------------------------------------------------------

! 1. define variables
  integer(kind=4) :: smonth
  character(len=2) :: month
  integer(kind=4) :: id_snow
  integer length, ierr
! 2. begin program
  print*,'  ... read NISE sfc snow and ice mask  '


!=== 2.1.get month 
call ICNVRT(0,smonth,month,length,ierr)
if (smonth<10) month='0'//month

!=== 2.2. open snow maks file
!call open_smk_file(month, id_snow)
call open_smk_nise_file(month, id_snow)

!=== 2.3. read albedo file
!call read_smk_file(id_snow)  
call read_smk_nise_file(id_snow) 

!=== 2.4. close albedo file
!call close_smk_file(id_snow)
call close_smk_nise_file(id_snow)

!print*,sat%snow_mask(1000,1000)
! 3. end subroutine      
end subroutine  read_sfc_snow_ice_mask
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~ subroutine  17:  open emiss file ~~~~~~~~~~~~~~~~~~~
subroutine  open_smk_nise_file(month, id)

!-----------------------------------------------------------------------
! !F90 open_smk_nise_file
!
! !description:
!    This program is to open the surface snow mask file.
!
! !Input  parameters:
!    month           =  month for corresponding emissivity data
!          
! !Output parameters:
!    none
!
!-----------------------------------------------------------------------

! 1. define variables
  character(len=2), intent(in)    :: month
  integer(kind=4), intent(out) :: id
  
  character(len=200) :: data_dir
  character(len=256) :: filename

  LOGICAL :: file_exists
  
! 2. begin program  
  data_dir = trim(code_root_path)//"coeff/sfc_snow_ice"
  filename = trim(data_dir)//"/NISE_SSMIF13_EASEGRID_M"//trim(month)//".HDF"

  inquire(file = filename, exist = file_exists)
  if (.not. file_exists) then
     print *,"(/,a,'NISE Surface snow mask file, ',a,' does not exist.')"
     stop
  endif
  
  id = sfstart(trim(filename), DFACC_READ)
  if (id == FAIL) then
     print *,"(/,a,'Failed to open, ',a)",trim(filename)
     stop
  endif
    
! 3. end subroutine      
end subroutine  open_smk_nise_file
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

!~~~~~~~~~~~~~~~~~~~ subroutine  18:  read smk file ~~~~~~~~~~~~~~~~~~~~~
subroutine  read_smk_nise_file(id)

!-----------------------------------------------------------------------
! !F90 read_smk_nise_file
!
! !description:
!    This program is to read the surface snow mask file.
!
! !Input  parameters:
!    id1/2           =  hdf file id
!          
! !Output parameters:
!    emiss           =  emissivity data
!
!-----------------------------------------------------------------------

! 1. define variables
  integer(kind=4), intent(in)                   :: id
!  real(kind=4),allocatable, dimension(:,:)  :: lat, lon
  integer(kind=1), allocatable, dimension(:,:)  :: space_mask
  integer(kind=1), allocatable, dimension(:,:)  :: smk_n
  integer(kind=1), allocatable, dimension(:,:)  :: smk_s
!  real(kind=4), dimension(num_lon_alb,num_lat_alb) :: alb_org2, alb_org3, alb_org4, alb_org5, alb_org6
 ! integer(kind=1), allocatable, dimension(:,:)  :: smk_org
  integer(kind=4) ::  nx, ny, i, j
  
  integer(kind=4) :: astatus

  character(len=14) :: sds_name_s
  character(len=14) :: sds_name_n
  character(len=1)  :: gridname

  
! 2. begin program

!=== 2.1. allocate variables
allocate(smk_n(num_lon_nise,num_lat_nise), &
         smk_s(num_lon_nise,num_lat_nise))

!=== 2.2. assign data
!lon        = sat%lon
!lat        = sat%lat
!space_mask = sat%space_mask
!smk = 1
sat%snow_mask = 1
nx            = missing_value_int4
ny            = missing_value_int4

!=== 2.3. read start
! read channel 2 and 5

  sds_name_n = 'NL_NISE_Extent'
  sds_name_s = 'SL_NISE_Extent'
  
  call read_smk_nise_hdf(id, trim(sds_name_n), smk_n) 
  call read_smk_nise_hdf(id, trim(sds_name_s), smk_s) 
    
  ! match lon lat 
  do i = 1, sat%nElem
  do j = 1, sat%nLine
  
     !if (geo%lon(i,j) >= -180.0 .and. geo%lon(i,j) <= 180.0 .and. &
     !    geo%lat(i,j) >= -90.0 .and. geo%lat(i,j) <= 90.0) then
     ! 
     !    nx = abs(INT((geo%lon(i,j) - first_lon_smk)/del_lon_smk))+1
     !    ny = abs(INT((first_lat_smk - geo%lat(i,j))/del_lat_smk))+1
     !   
     !    sat%snow_mask(i,j) = smk(nx,ny)
     !
     !endif
     if  (geo%lon(i,j) >= -180.0 .and. geo%lon(i,j) <= 180.0 .and. &
          geo%lat(i,j) >= 0.0 .and. geo%lat(i,j) <= 90.0) then
          gridname = 'N'
          call nise_ezlh_convert (gridname, geo%lon(i,j), geo%lat(i,j), nx, ny)
          sat%snow_mask(i,j) = smk_n(nx,ny)
     endif
     if  (geo%lon(i,j) >= -180.0 .and. geo%lon(i,j) <= 180.0 .and. &
          geo%lat(i,j) >= -90.0 .and. geo%lat(i,j) < 0.0) then
          gridname = 'S'
          call nise_ezlh_convert (gridname, geo%lon(i,j), geo%lat(i,j), nx, ny)
          sat%snow_mask(i,j) = smk_s(nx,ny)
     endif     
  enddo
  enddo

  
!=== 2.4. deallocate variables
  deallocate(smk_n,        &
             smk_s,        &
             stat=astatus)
  if (astatus /= 0) then
     print *,"(a,'Error deallocating NISE surface snow mask data file.')"
     stop
  endif
  
! 3. end subroutine   
end subroutine  read_smk_nise_file
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

!~~~~~~~~~~~~~~~~~~~ subroutine  19:  read alb hdf file ~~~~~~~~~~~~~~~~~
subroutine  read_smk_nise_hdf(sd_id, sds_name, smk_org)

!-----------------------------------------------------------------------
! !F90 read_smk_nise_hdf
!
! !description:
!    This program is to read the surface albedo hdf file.
!
! !Input  parameters:
!    sd_id           =  hdf file id
!    sds_name        =  satellite channel number name
!          
! !Output parameters:
!    alb_org           =  albedo data
!
!-----------------------------------------------------------------------

! 1. define variables
!===== 1.1.hdf4 Function declaration.
character(*), intent(in) :: sds_name

integer(kind=4) sd_id, sds_id, attr_id
integer(kind=4), dimension(2) :: start, stride, edges

integer(kind=1),allocatable, dimension(:,:) :: smk_org1
!real(kind=4), dimension(num_lon_alb,num_lat_alb) :: alb_org
integer(kind=1),dimension(:,:), intent(OUT) :: smk_org

!===== 1.2. other
integer(kind=4) :: status  
character*100 message

!*******
! 2. begin program
  !===== 2.1. initialize
  allocate(smk_org1(num_lon_nise,num_lat_nise))
  start    = 0
  stride   = 1
  edges(1) = num_lon_nise
  edges(2) = num_lat_nise ! 16 column * 1 row
  status   = 0
  message = 'NISE Surface snow mask'

  ! read albedo
  sds_id = sfselect(sd_id, sfn2index(sd_id,sds_name))   
  status = sfrdata (sds_id, start, stride, edges, smk_org1)
  call hdf_info_message(message, sds_name, status) 
  
  status = sfendacc(sds_id)
  
  smk_org = smk_org1

  deallocate(smk_org1)

! 3. end subroutine     
end subroutine  read_smk_nise_hdf
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

!~~~~~~~~~~~~~~~~~~~ subroutine  20:  close smk hdf file ~~~~~~~~~~~~~~~~
subroutine  close_smk_nise_file(id)

!-----------------------------------------------------------------------
! !F90 close_smk_nise_file
!
! !description:
!    This program is to close the surface albedo hdf file.
!
! !Input  parameters:
!    id1/2           =  hdf file id
!          
! !Output parameters:
!    none
!
!-----------------------------------------------------------------------

  integer(kind=4), intent(in) :: id
  
  integer(kind=4) :: istatus
  integer(kind=4) :: sfend
  
  istatus = sfend(id)
  if (istatus /= 0) then
     print *,"(/,a,'Error closing NISE surface snow mask hdf file.')"
     stop
  endif

end subroutine  close_smk_nise_file
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
subroutine  nise_ezlh_convert (grid, lon, lat, col, row)

!--------------------------------------------------------------------------
!	character*(*) grid
!	real lat, lon, r, s

!	convert geographic coordinates (spherical earth) to 
!	azimuthal equal area or equal area cylindrical grid coordinates
!
!	status = ezlh_convert (grid, lat, lon, r, s)
!
!	input : grid - projection name '[NSM][lh]'
!               where l = "low"  = 25km resolution
!                     h = "high" = 12.5km resolution
!		lat, lon - geo. coords. (decimal degrees)
!
!	output: r, s - column, row coordinates
!
!	result: status = 0 indicates normal successful completion
!			-1 indicates error status (point not on grid)
!
!--------------------------------------------------------------------------
!   integer cols, rows, scale
!	real Rg, phi, lam, rho
    character(len=1), intent(in) :: grid
    real    :: RE_km, CELL_km, COS_PHI1, pai, Rg, r0, s0
    real    :: phi, lam, rho, r, s, lon, lat
    integer :: ezlh_convert, scale
    integer :: col, row,  cols, rows
    
!	radius of the earth (km), authalic sphere based on International datum 
	RE_km = 6371.228
!	nominal cell size in kilometers
	CELL_km = 25.067525

!	scale factor for standard paralles at +/-30.00 degrees
	COS_PHI1 = 0.866025403

	!pai = 3.141592653589793
	pai = PI
    !rad = inline('t*pai./180','pai','t');
	!deg = inline('t*180./pai','pai','t');
   
    ezlh_convert = -1

	!if (strcmp(grid(1:1), 'N')==1 | strcmp(grid(1:1),'S')==1) 
	!  cols = 721;
	!  rows = 721;
	!elseif (strcmp(grid(1:1),'M')==1) 
	!  cols = 1383;
	!  rows = 586;
	!else
	!  fprintf('!c20 !c20','ezlh_convert: unknown projection: ', grid);
	!  return;
    !end
    cols = num_lat_nise
    rows = num_lon_nise

	!if (strcmp(grid(2:2),'l')==1)
	!  scale = 1;
	!elseif (strcmp(grid(2:2),'h')==1) then
	!  scale = 2;
	!else
	!  fprintf('!c20 !c20','ezlh_convert: unknown projection: ', grid);
	!  return;
    !end
    scale = 1

    Rg = scale * RE_km/CELL_km

!
!	r0,s0 are defined such that cells at all scales
!	have coincident center points
!
    r0 = (cols-1)/2. * scale
    s0 = (rows-1)/2. * scale

    phi = lat*pai/180. !rad(lat,pai)
    lam = lon*pai/180. !rad(lon,pai)

    if (grid=='N') then 
	   rho = 2 * Rg * sin(PI/4. - phi/2.)
	   r = r0 + rho * sin(lam)
	   s = s0 + rho * cos(lam)
    endif
    if (grid=='S') then 
	   rho = 2 * Rg * cos(PI/4. - phi/2.)
	   r = r0 + rho * sin(lam)
	   s = s0 - rho * cos(lam)
    endif
    row = max(1,min((nint(s)+1), rows))
    col = max(1,min((nint(r)+1), cols))
    ezlh_convert = -1
  
end subroutine  nise_ezlh_convert
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

!~~~~~~~~~~~~~~~~~~~ subroutine  17:  open emiss file ~~~~~~~~~~~~~~~~~~~
subroutine  open_smk_file(month, id)

!-----------------------------------------------------------------------
! !F90 open_smk_file
!
! !description:
!    This program is to open the surface snow mask file.
!
! !Input  parameters:
!    month           =  month for corresponding emissivity data
!          
! !Output parameters:
!    none
!
!-----------------------------------------------------------------------

! 1. define variables
  character(len=2), intent(in)    :: month
  integer(kind=4), intent(out) :: id
  
  character(len=200) :: data_dir
  character(len=256) :: filename

  LOGICAL :: file_exists
  
! 2. begin program  
  data_dir = trim(code_root_path)//"coeff/sfc_snow_ice"
  filename = trim(data_dir)//"/NISE_SSMif13_GLL025Deg_M"//trim(month)//".hdf"

  inquire(file = filename, exist = file_exists)
  if (.not. file_exists) then
     print *,"(/,a,'Surface snow mask file, ',a,' does not exist.')"
     stop
  endif
  
  id = sfstart(trim(filename), DFACC_READ)
  if (id == FAIL) then
     print *,"(/,a,'Failed to open, ',a)",trim(filename)
     stop
  endif
    
! 3. end subroutine      
end subroutine  open_smk_file
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

!~~~~~~~~~~~~~~~~~~~ subroutine  18:  read smk file ~~~~~~~~~~~~~~~~~~~~~
subroutine  read_smk_file(id)

!-----------------------------------------------------------------------
! !F90 read_alb_file
!
! !description:
!    This program is to read the surface snow mask file.
!
! !Input  parameters:
!    id1/2           =  hdf file id
!          
! !Output parameters:
!    emiss           =  emissivity data
!
!-----------------------------------------------------------------------

! 1. define variables
  integer(kind=4), intent(in)                   :: id
!  real(kind=4),allocatable, dimension(:,:)  :: lat, lon
  integer(kind=1), allocatable, dimension(:,:)  :: space_mask
  integer(kind=1), allocatable, dimension(:,:)  :: smk
!  real(kind=4), dimension(num_lon_alb,num_lat_alb) :: alb_org2, alb_org3, alb_org4, alb_org5, alb_org6
 ! integer(kind=1), allocatable, dimension(:,:)  :: smk_org
  integer(kind=4) ::  nx, ny, i, j
  
  integer(kind=4) :: astatus

  character(len=16) :: sds_name

  
! 2. begin program

!=== 2.1. allocate variables
allocate(smk(num_lon_smk,num_lat_smk))

!=== 2.2. assign data
!lon        = sat%lon
!lat        = sat%lat
!space_mask = sat%space_mask
!smk = 1
sat%snow_mask = 1
nx            = missing_value_int4
ny            = missing_value_int4

!=== 2.3. read start
! read channel 2 and 5

  sds_name = 'snow_mask'

  call read_smk_hdf(id, trim(sds_name), smk) 

    
  ! match lon lat 
  do i = 1, sat%nElem
  do j = 1, sat%nLine
  
     if (geo%lon(i,j) >= -180.0 .and. geo%lon(i,j) <= 180.0 .and. &
         geo%lat(i,j) >= -90.0 .and. geo%lat(i,j) <= 90.0) then
     
         nx = abs(INT((geo%lon(i,j) - first_lon_smk)/del_lon_smk))+1
         ny = abs(INT((first_lat_smk - geo%lat(i,j))/del_lat_smk))+1
        
         sat%snow_mask(i,j) = smk(nx,ny)
   
     endif
     
  enddo
  enddo

  
!=== 2.4. deallocate variables
  deallocate(smk,        &
             stat=astatus)
  if (astatus /= 0) then
     print *,"(a,'Error deallocating surface snow mask data file.')"
     stop
  endif
  
! 3. end subroutine   
end subroutine  read_smk_file
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 


!~~~~~~~~~~~~~~~~~~~ subroutine  19:  read alb hdf file ~~~~~~~~~~~~~~~~~
subroutine  read_smk_hdf(sd_id, sds_name, smk_org)

!-----------------------------------------------------------------------
! !F90 read_alb_hdf
!
! !description:
!    This program is to read the surface albedo hdf file.
!
! !Input  parameters:
!    sd_id           =  hdf file id
!    sds_name        =  satellite channel number name
!          
! !Output parameters:
!    alb_org           =  albedo data
!
!-----------------------------------------------------------------------

! 1. define variables
!===== 1.1.hdf4 Function declaration.
character(*), intent(in) :: sds_name

integer(kind=4) sd_id, sds_id, attr_id
integer(kind=4), dimension(2) :: start, stride, edges

integer(kind=1),allocatable, dimension(:,:) :: smk_org1
!real(kind=4), dimension(num_lon_alb,num_lat_alb) :: alb_org
integer(kind=1),dimension(:,:), intent(OUT) :: smk_org

!===== 1.2. other
integer(kind=4) :: status  
character*100 message

!*******
! 2. begin program
  !===== 2.1. initialize
  allocate(smk_org1(num_lon_smk,num_lat_smk))
  start    = 0
  stride   = 1
  edges(1) = num_lon_smk
  edges(2) = num_lat_smk ! 16 column * 1 row
  status   = 0
  message = 'Surface snow mask'

  ! read albedo
  sds_id = sfselect(sd_id, sfn2index(sd_id,sds_name))   
  status = sfrdata (sds_id, start, stride, edges, smk_org1)
  call hdf_info_message(message, sds_name, status) 
  
  status = sfendacc(sds_id)
  
  smk_org = smk_org1

  deallocate(smk_org1)

! 3. end subroutine     
end subroutine  read_smk_hdf
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 


!~~~~~~~~~~~~~~~~~~~ subroutine  20:  close smk hdf file ~~~~~~~~~~~~~~~~
subroutine  close_smk_file(id)

!-----------------------------------------------------------------------
! !F90 close_alb_file
!
! !description:
!    This program is to close the surface albedo hdf file.
!
! !Input  parameters:
!    id1/2           =  hdf file id
!          
! !Output parameters:
!    none
!
!-----------------------------------------------------------------------

  integer(kind=4), intent(in) :: id
  
  integer(kind=4) :: istatus
  integer(kind=4) :: sfend
  
  istatus = sfend(id)
  if (istatus /= 0) then
     print *,"(/,a,'Error closing surface snow mask hdf file.')"
     stop
  endif

end subroutine  close_smk_file
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~ subroutine  20:  close smk hdf file ~~~~~~~~~~~~~~~~
subroutine  read_ecosystem_file(fname)

!-----------------------------------------------------------------------
! !F90 read_ecosystem_file
!
! !description:
!    This program is to close the surface albedo hdf file.
!
! !Input  parameters:
!    npixel     = sat%nElem
!    scans_cube = sat%nLine/10
!          
! !Output parameters:
!    cube_eco
!
!-----------------------------------------------------------------------
! 1. define variables
    character*(*)             :: fname
    character (len = 200)             :: eco_mapFN
    integer                           :: status
    !Stores the resolution (in degrees) of the global map
    real     (kind = 4)               :: Map_Resolution
    !Indices for section of global map that is needed for supplying values for granule:
    !Global map dims and maximum number of albedo wavelengths:
    integer ,parameter                :: NumMapCols        = 43200, &
                                         NumMapRows        = 21600

    !Number of Ecosystem classifications
    integer , parameter               :: NumEcosystems = 18
    
! --- constants
!      integer, parameter :: max_row = 17347
!      integer, parameter :: max_col = 40031
!      integer, parameter :: recsize = 512
!      integer, parameter :: max_count = 3000000    
!      integer, parameter :: max_rec = (max_row*max_col )/recsize
!      integer(KIND=4) :: io_err,i,j
!      integer(kind=1), allocatable,dimension(:,:) :: dd       
    
    ! 2. begin program
    print*,'  ... read ecosystem map information  '
    
    eco_mapFN = trim(code_root_path)//'coeff/'//trim(fname)

!  read IGBPmap
    !*****************************************************************************
    ! Set up the global albedo/ecosystem map resolution and determine the start/stop
    !  points for the portion of the global map to be read in.  This is done by
    !  computing the map indices from the min/max lat/lon values.
    !*****************************************************************************
    !Determine the map resolution:
    Map_Resolution = 360.000d0 / float(NumMapCols)            !???????360.000d0 

    allocate( cube_eco (1:NumMapCols, 1:NumMapRows) )
    !*****************************************************************************
    ! Read in the portion of the ecosystem map from the corresponding file:
    !*****************************************************************************
    call ReadEcosystemStats(eco_mapFN, NumMapCols, NumMapRows, status)


end subroutine  read_ecosystem_file
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~ subroutine  20:  ReadEcosystemStats ~~~~~~~~~~~~~~~~
subroutine  ReadEcosystemStats(EcosystemFN, NumMapCols, NumMapRows, status)

    character (len = *),                intent ( in)  :: EcosystemFN
    character (len = 25)                              :: SDSName
    integer,                            intent (out)  :: status
    !File IDs
    integer                  :: EcoMapFID
    !SDS IDs
    integer                  :: EcoMapSDSID
    !HDF IDs
    integer                  :: HDFstatus
    integer  , intent(in)    :: NumMapCols,      &
                                NumMapRows
    integer                  :: hdfStart(2), hdfStride(2), hdfEdge(2)
    ! match
    integer :: i,j,nlat,nlon
    real    :: dlat,dlon,first_lat,first_lon,rlat,rlon, mx, my
    integer :: factor

    !begin program
    !Define the SDS Name:
    hdfStart   = 0
    hdfStride  = 1
    hdfEdge(1) = NumMapCols
    hdfEdge(2) = NumMapRows

    SDSName = "ecosystem_type"

    ! Open the HDF file for reading:
    EcoMapFID = Sfstart(trim(EcosystemFN), Dfacc_Read)    !?????????Sfstart 
    if (EcoMapFID == Fail) then
       status = failure
       call local_message_handler('Problem detected ecosysyem file sfstart',status,'getAlbedoEco')
       return
    end if

    ! Obtain the data set ID's:
    EcoMapSDSID   = SFselect(EcoMapFID, SFn2index(EcoMapFID, trim(SDSName)))   !?????????SFselect( ≤√¥∫Ø ˝
    if (EcoMapSDSID == Fail) then
       status = failure
       call local_message_handler('Problem detected ecosystem file sfselect',status,'getAlbedoEco')
       return
    end if

    ! Read in the data:
    HDFstatus = SFrdata(EcoMapSDSID, hdfStart, hdfStride, hdfEdge, cube_eco)
    !Error checking:
    if (HDFstatus == FAIL) then
       status = failure
       call local_message_handler('Problem detected ecosystem file sfrdata',status,'getAlbedoEco')
       return
    end if

    !end Access to the datasets and files:
    HDFstatus = SFendacc( EcoMapSDSID )
    if (HDFstatus == FAIL) then
       status = failure
       call local_message_handler('Problem detected ecosystem file sfendacc',status,'getAlbedoEco')
       return
    end if
    
    HDFstatus = SFend ( EcoMapFID   )
    if (HDFstatus == FAIL) then
       status = failure
       call local_message_handler('Problem detected ecossystem file sfend',status,'getAlbedoEco')
       return
    end if

! match    
nlat = NumMapRows
nlon = NumMapCols
dlat = 180.0 / nlat
dlon = 360. / nlon
first_lat = 89.99583
first_lon = -179.99583
factor = 1
if (first_lat > 0.0 ) factor = -1

do j = 1, sat%nLine
   do i = 1, sat%nElem
   
      rlon = geo%lon(i,j)
      rlat = geo%lat(i,j)
    
      my = max(1, min(nlat, int((rlat-first_lat+0.5*dlat*factor)/(dlat*factor) +1)))
      mx = max(1, min(nlon, int((rlon-first_lon+0.5*dlon)/(dlon) +1)))
    
      sat%eco(i,j) = cube_eco(mx,my)
      !print*,i,j,rlon,rlat,sat%eco(i,j) 
   enddo
enddo

! deallocate IGBP data arrays
call deallocate_eco_arrays()

end subroutine  ReadEcosystemStats
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~ subroutine  20:  close smk hdf file ~~~~~~~~~~~~~~~~
subroutine  deallocate_eco_arrays()

  integer(kind=4) :: astatus
  
  deallocate(cube_eco,      &
             stat=astatus)
             
  if (astatus /= 0) then
     print *,"(a,'Error deallocating IGBP data file.')"
     stop
  endif

end subroutine  deallocate_eco_arrays
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~ subroutine  20:  close smk hdf file ~~~~~~~~~~~~~~~~
subroutine  read_oisst_file()
!-----------------------------------------------------------------------
! !F90 read_oisst_data_main
!
! !description:
!    The main program for reading nwp DATA.
!    Here, this program reads two neighbouring nwp DATA
!    The NCEP DATA are used as nwp DATA.  
!
! !Input  parameters:
!    none
!
! !Output parameters:
!    none
!
!-----------------------------------------------------------------------

!use names_module ,ONLY: nwp_data_path, code_root_path, fylat_nwp_opt


! 1. define variables
!===== 1.1.input variables
integer(kind=4) :: n  ! 1 is temporal / 2 is real

!===== 1.2.middle variables  
character*300 bin_name
integer sd_id, sds_id, attr_id, status
integer, dimension(2) :: start1, stride1
integer (HID_T)       :: sd_id5, sds_id5, attr_id5
integer (HSIZE_T), dimension(2) ::  edges25
integer (HSIZE_T), dimension(2) :: dims_geo2  
real(kind=4),dimension(:,:),allocatable :: sst0, sst1

    ! match
    integer :: i,j,nlat,nlon
    real    :: dlat,dlon,first_lat,first_lon,rlat,rlon, mx, my
    integer :: factor
    
! 2. begin program

  print*,'  ... read oisst hdf5 (daily) data'
  
  allocate(sst0(num_lon_oisst,num_lat_oisst), sst1(num_lon_oisst,num_lat_oisst))

  start1   = 0
  stride1  = 1
  !edges25(1)= num_lon_oisst
  !edges25(2)= num_lat_oisst
  
  !print*,trim(oisst_data)
  call h5open_f(status)
  call h5fopen_f(trim(oisst_data), H5F_ACC_RDONLY_F, sd_id5, status)
  
  call h5dopen_f(sd_id5, "sst", sds_id5, status)
  call h5dread_f(sds_id5, H5T_NATIVE_real, sst0, edges25, status)
  if (status /= 0) then 
     print*,'ERROR: oisst hdf5 read failed !'
     stop
  endif
  call h5dclose_f(sds_id5, status)
  
  ! Terminate access to the SD interface and close the file. 
  call h5fclose_f(sd_id5, status)  ! close 
  call h5close_f(status)

  !===== 2.3. find x, y nwp position
sst1(1:int(num_lon_oisst/2),:) = sst0(int(num_lon_oisst/2)+1:num_lon_oisst,:)
sst1(int(num_lon_oisst/2)+1:num_lon_oisst,:) = sst0(1:int(num_lon_oisst/2),:)

! match    
nlat = num_lat_oisst
nlon = num_lon_oisst
dlat = 180.0 / nlat
dlon = 360. / nlon
!first_lat = 89.75
first_lat = -89.75   ! revised by minmin 20221014
first_lon = -179.75
factor = 1
if (first_lat > 0.0 ) factor = -1

do j = 1, sat%nLine
   do i = 1, sat%nElem
   
      rlon = geo%lon(i,j)
      rlat = geo%lat(i,j)
    
      my = max(1, min(nlat, int((rlat-first_lat+0.5*dlat*factor)/(dlat*factor) +1)))
      mx = max(1, min(nlon, int((rlon-first_lon+0.5*dlon)/(dlon) +1)))
      
      sat%sst(i,j) = -999.0
      if (sst1(mx,my) > -20) then
         sat%sst(i,j) = sst1(mx,my)+273.15
      endif
      !print*,i,j,rlon,rlat,sat%eco(i,j) 
   enddo
enddo

deallocate(sst0, sst1)
! 3. end subroutine    
end subroutine  read_oisst_file
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~ subroutine  1:  read_emissivity_main ~~~~~~~~~~~~~~~
subroutine  read_emissivity_data(month)

!-----------------------------------------------------------------------
! !F90 read_emissivity_data_main
!
! !Description:
!    This is a main program for reading surface emissivity file.
!
! !Input  parameters:
!    month           =  month for corresponding emissivity data
!          
! !Output parameters:
!    none
!
!-----------------------------------------------------------------------

! 1. define variables
  integer(kind=4), intent(in) :: month
  integer(kind=4) :: id
  
! 2. begin program
  print*,'  ... read IR emissivity hdf4 (monthly) data'


!=== 2.1. deal internal memory
call allocate_emiss_arrays

!=== 2.2. open emissivity file
call open_emiss_file(month, id)

!=== 2.3. read emissivity file
call read_emiss_file(id)  

!=== 2.4. close emissivity file
call close_emiss_file(id)

!month_id = month

! 3. end subroutine      
end subroutine  read_emissivity_data
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

!~~~~~~~~~~~~~~~~~~~ subroutine  2:  open emiss file ~~~~~~~~~~~~~~~~~~~~
subroutine  open_emiss_file(month, id)

!-----------------------------------------------------------------------
! !F90 open_emiss_file
!
! !Description:
!    This program is to open the surface emissivity file.
!
! !Input  parameters:
!    month           =  month for corresponding emissivity data
!          
! !Output parameters:
!    none
!
!-----------------------------------------------------------------------

! 1. define variables
  integer(kind=4), intent(in) :: month
  integer(kind=4), intent(out) :: id
  integer(kind=4) :: i, j, k
  
  character(len=200) :: data_dir
  character(len=256) :: filename
  character(len=3)   :: jday_str
  character(len=4)   :: year_str
  
  logical :: file_exists

! 2. begin program  
  year_str = "2005"
  
  select case (month)
  case (1)
    jday_str = "001"
  case (2)
    jday_str = "032"
  case (3)
    jday_str = "060"
  case (4)
    jday_str = "091"
  case (5)
    jday_str = "121"
  case (6)
    jday_str = "152"
  case (7)
    jday_str = "182"
  case (8)
    jday_str = "213"
  case (9)
    jday_str = "244"
  case (10)
    jday_str = "274"
  case (11)
    jday_str = "305"
  case (12)
    jday_str = "335"
  end select
  
  data_dir = trim(code_root_path)//"/coeff/sfc_emiss"
  filename = trim(data_dir)//"/global_emiss_intABI_"//trim(year_str)//trim(jday_str)//".hdf"
  
  inquire(file = filename, exist = file_exists)
  if (.not. file_exists) then
     print *,"(/,a,'Surface emissivity file, ',a,' does not exist.')"
     stop
  endif
  
  id = sfstart(trim(filename), DFACC_READ)
  if (id == FAIL) then
     print *,"(/,a,'Failed to open, ',a)",trim(filename)
     stop
  endif
  
! 3. end subroutine      
end subroutine  open_emiss_file
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

!~~~~~~~~~~~~~~~~~~~ subroutine  3:  read emiss file ~~~~~~~~~~~~~~~~~~~~
subroutine  read_emiss_file(id)

!-----------------------------------------------------------------------
! !F90 read_emiss_file
!
! !Description:
!    This program is to read the surface emissivity file.
!
! !Input  parameters:
!    id              =  hdf file id
!          
! !Output parameters:
!    emiss           =  emissivity data
!
!-----------------------------------------------------------------------

! 1. define variables
  integer(kind=4), intent(in)                  :: id
  real(kind=4),allocatable, dimension(:,:)     :: lat, lon
  integer(kind=1), allocatable, dimension(:,:) :: space_mask
  real(kind=4), allocatable, dimension(:,:)    :: emiss
!  real(kind=4), dimension(num_lon_emiss,num_lat_emiss) :: emiss_org
  real(kind=4), allocatable, dimension(:,:) :: emiss_org
  integer(kind=4) ::  nx, ny, i, j ,k 
  
  integer :: astatus
  character(len=2) :: ichan
  character(len=10) :: sds_name
  
! 2. begin program

!=== 2.1. allocate variables
allocate(lat(sat%nElem,sat%nLine),lon(sat%nElem,sat%nLine), &
         emiss(sat%nElem,sat%nLine),          &
         emiss_org(num_lon_emiss,num_lat_emiss))

!=== 2.2. assign data
lon        = geo%lon
lat        = geo%lat
emiss      = missing_value_real4
nx         = missing_value_real4
ny         = missing_value_real4

!=== 2.3. read start
! read channel 7 to 16
do k=7,16

  if (k<10) then
     write(ichan,'(I1)')k
  else 
     write(ichan,'(I2)')k
  endif

  sds_name = 'emiss'//trim(ichan)
  call read_emiss_hdf(id, trim(sds_name), emiss_org) 

  ! match lon lat 
  do i=1,sat%nElem
  do j=1,sat%nLine

 !    if (space_mask(i,j) /= sym%SPACE ) then
     emiss(i,j) = 1.0
     if (lon(i,j)<=180 .and. lon(i,j)>=-180) then
        nx = INT((lon(i,j) - first_lon_emiss)/del_lon_emiss)+1
        ny = INT((first_lat_emiss - lat(i,j))/del_lat_emiss)+1
        
        if (emiss_org(nx,ny) > 1. .or. emiss_org(nx,ny) <= 0.) then
            emiss_org(nx,ny) = emiss_water_sfc(k)
        endif
        
        emiss(i,j) = emiss_org(nx,ny)
     endif
!     endif
     
  end do
  end do

  select case (k)
  case (7)
    sat%sfc_emiss38  = emiss
    sat%sfc_emiss40  = emiss
  case (8)
   ! sat%sfc_emiss8  = emiss
  case (9)
   ! sat%sfc_emiss9  = emiss
  case (10)
    sat%sfc_emiss73 = emiss
!    print*,'emiss73', sat%sfc_emiss73(1000,1000:1011)
  case (11)
    sat%sfc_emiss86 = emiss
  case (12)
   ! sat%sfc_emiss12 = emiss
  case (13)
   ! sat%sfc_emiss13 = emiss
  case (14)
    sat%sfc_emiss11 = emiss
  case (15)
    sat%sfc_emiss12 = emiss
  case (16)
   ! sat%sfc_emiss16 = emiss
  end select
  
enddo ! do k=7,16
  
!=== 2.4. deallocate variables
deallocate(lat,lon,emiss,emiss_org, stat=astatus)
if (astatus /= 0) then
   print *,"(a,'Error deallocating emissivity data file.')"
   stop
endif
  
! 3. end subroutine   
end subroutine  read_emiss_file
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

!~~~~~~~~~~~~~~~~~~~ subroutine  4:  read emiss hdf file ~~~~~~~~~~~~~~~~
subroutine  read_emiss_hdf(sd_id, sds_name, emiss_org)

!-----------------------------------------------------------------------
! !F90 read_emiss_hdf
!
! !Description:
!    This program is to read the surface emissivity hdf file.
!
! !Input  parameters:
!    sd_id           =  hdf file id
!    sds_name        =  satellite channel number name
!          
! !Output parameters:
!    emiss           =  emissivity data
!
!-----------------------------------------------------------------------

! 1. define variables
!===== 1.1.hdf4 Function declaration.
character(*), intent(in) :: sds_name

integer(kind=4) :: sd_id, sds_id, attr_id
integer(kind=4), dimension(2) :: start, stride, edges

real(kind=4) :: slope, offset
integer(kind=2),allocatable, dimension(:,:) :: emiss_org1
!real(kind=4), dimension(num_lon_emiss,num_lat_emiss) :: emiss_org
real(kind=4), dimension(:,:),intent(OUT) :: emiss_org

!===== 1.2. other
integer(kind=4) :: status  
character*100 message

!*******
! 2. begin program
  !===== 2.1. initialize
  !print*,size(emiss_org)
  allocate(emiss_org1(num_lon_emiss,num_lat_emiss))
  
  start    = 0
  stride   = 1
  edges(1) = num_lon_emiss
  edges(2) = num_lat_emiss ! 16 column * 1 row
  status   = 0
  message = 'Surface Emissivity'

  ! read emissivity
  sds_id = sfselect(sd_id, sfn2index(sd_id,sds_name))   
  status = sfrdata (sds_id, start, stride, edges, emiss_org1)
  call hdf_info_message(message, sds_name, status) 
  
  attr_id = sffattr(sds_id, "scale_factor")
  status  = sfrattr(sds_id, attr_id, slope)
  call hdf_info_message(message, sds_name//" scale_factor", status)
  
  attr_id = sffattr(sds_id, "add_offset")
  status  = sfrattr(sds_id, attr_id, offset)
  call hdf_info_message(message, sds_name//" add_offset", status)
  
  status = sfendacc(sds_id)
  
  emiss_org = emiss_org1*slope + offset

  deallocate(emiss_org1)

! 3. end subroutine     
end subroutine  read_emiss_hdf
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

!~~~~~~~~~~~~~~~~~~~ subroutine  5:  close emiss hdf file ~~~~~~~~~~~~~~~
subroutine  close_emiss_file(id)

!-----------------------------------------------------------------------
! !F90 close_emiss_file
!
! !Description:
!    This program is to close the surface emissivity hdf file.
!
! !Input  parameters:
!    id           =  hdf file id
!          
! !Output parameters:
!    none
!
!-----------------------------------------------------------------------

  integer(kind=4), intent(in) :: id
  
  integer(kind=4) :: istatus
  
  istatus = sfend(id)
  if (istatus /= 0) then
     print *,"(/,a,'Error closing surface emissivity hdf file.')"
     stop
  endif

end subroutine  close_emiss_file
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	

!~~~~~~~~~~~~~~~~~~~ subroutine  6: allocate_emiss_data ~~~~~~~~~~~~~~~~~
subroutine  allocate_emiss_arrays

!-----------------------------------------------------------------------
! !F90 allocate_emiss_arrays
!
! !Description:
!    This program is to allocate emiss data.
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
!  allocate(sat%sfc_emiss7(sat%nElem,sat%nLine), sat%sfc_emiss8(sat%nElem,sat%nLine),   &
!           sat%sfc_emiss9(sat%nElem,sat%nLine), sat%sfc_emiss10(sat%nElem,sat%nLine),  &
!           sat%sfc_emiss11(sat%nElem,sat%nLine),sat%sfc_emiss12(sat%nElem,sat%nLine),  &
!           sat%sfc_emiss13(sat%nElem,sat%nLine),sat%sfc_emiss14(sat%nElem,sat%nLine),  &
!           sat%sfc_emiss15(sat%nElem,sat%nLine),sat%sfc_emiss16(sat%nElem,sat%nLine),  &
!           stat=astatus)

  if (sat%chan_flag(20) > 0) then 
     allocate(sat%sfc_emiss38(sat%nElem,sat%nLine), &
              stat=astatus)
  endif
  if (sat%chan_flag(21) > 0) then 
     allocate(sat%sfc_emiss40(sat%nElem,sat%nLine), &
              stat=astatus)
  endif
  if (sat%chan_flag(22) > 0) then 
     allocate(sat%sfc_emiss73(sat%nElem,sat%nLine), &
              stat=astatus)
  endif
  if (sat%chan_flag(23) > 0) then 
     allocate(sat%sfc_emiss86(sat%nElem,sat%nLine), &
              stat=astatus)
  endif
  if (sat%chan_flag(24) > 0) then 
     allocate(sat%sfc_emiss11(sat%nElem,sat%nLine), &
              stat=astatus)
  endif
  if (sat%chan_flag(25) > 0) then 
     allocate(sat%sfc_emiss12(sat%nElem,sat%nLine), &
              stat=astatus)
  endif
           
  if (astatus /= 0) then
     print *,"(a,'Not enough memory to allocate emiss data structure.')"
     stop
  endif

! 3. end subroutine 
end subroutine  allocate_emiss_arrays
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~ subroutine  7: deallocate smiss data ~~~~~~~~~~~~~~~
subroutine  deallocate_emiss_arrays()
 
!-----------------------------------------------------------------------
! !F90 deallocate_emiss_data
!
! !Description:
!    This program is to deallocate emiss data.
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
!  deallocate(sat%sfc_emiss7, sat%sfc_emiss8,   &
!             sat%sfc_emiss9, sat%sfc_emiss10,  &
!             sat%sfc_emiss11,sat%sfc_emiss12,  &
!             sat%sfc_emiss13,sat%sfc_emiss14,  &
!             sat%sfc_emiss15,sat%sfc_emiss16,  &
!             stat=astatus)

  if (sat%chan_flag(20) > 0) then 
     deallocate(sat%sfc_emiss38, &
                stat=astatus)
  endif
  if (sat%chan_flag(21) > 0) then 
     deallocate(sat%sfc_emiss40, &
                stat=astatus)
  endif
  if (sat%chan_flag(22) > 0) then 
     deallocate(sat%sfc_emiss73, &
               stat=astatus)
  endif
  if (sat%chan_flag(23) > 0) then 
     deallocate(sat%sfc_emiss86, &
                stat=astatus)
  endif
  if (sat%chan_flag(24) > 0) then 
     deallocate(sat%sfc_emiss11, &
                stat=astatus)
  endif
  if (sat%chan_flag(25) > 0) then 
     deallocate(sat%sfc_emiss12, &
                stat=astatus)
  endif
               
  if (astatus /= 0) then
     print *,"(a,'Error deallocate emiss data structure.')"
     stop
  endif
 
! 3. end subroutine 
end subroutine  deallocate_emiss_arrays
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~ subroutine  9:  read_alb_main ~~~~~~~~~~~~~~~~~~~~~~
subroutine  read_albedo_data(daynum)

!-----------------------------------------------------------------------
! !F90 read_AlbMap_main
!
! !Description:
!    This is a main program for reading surface albedo file.
!
! !Input  parameters:
!    day             =  day for corresponding day data
!          
! !Output parameters:
!    none
!
!-----------------------------------------------------------------------

! 1. define variables
  integer(kind=4), intent(in) :: daynum
  integer(kind=4) :: day_alb_id
  integer(kind=4) :: id2, id3, id4, id5, id6

! 2. begin program
    print*,'  ... read ws albedo information  '


call deal_alb_dayid(daynum, day_alb_id)

!=== 2.1. deal internal memory
call allocate_alb_arrays

!=== 2.2. open albedo file
call open_alb_file(day_alb_id, id2, id3, id4, id5, id6)

!=== 2.3. read albedo file
call read_alb_file(id2, id3, id4, id5, id6)  

!=== 2.4. close albedo file
call close_alb_file(id2, id3, id4, id5, id6)

!day_id = day_alb_id


! 3. end subroutine      
end subroutine  read_albedo_data
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~ subroutine  10: allocate_alb_data ~~~~~~~~~~~~~~~~~~
subroutine  allocate_alb_arrays()

!-----------------------------------------------------------------------
! !F90 allocate_alb_arrays
!
! !Description:
!    This program is to allocate alb data.
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
  !allocate(sat%ws_albedo66(xr,yr),   &
  !         sat%ws_albedo87(xr,yr),   &
  !         sat%ws_albedo124(xr,yr),   &
  !         sat%ws_albedo164(xr,yr),   &
  !         sat%ws_albedo213(xr,yr),   &
  !         stat=astatus)

  if (sat%chan_flag(3) > 0) then 
     allocate(sat%ws_albedo66(sat%nElem,sat%nLine), &
              stat=astatus)
  endif
  if (sat%chan_flag(4) > 0) then 
     allocate(sat%ws_albedo87(sat%nElem,sat%nLine), &
              stat=astatus)
  endif
  if (sat%chan_flag(5) > 0) then 
     allocate(sat%ws_albedo124(sat%nElem,sat%nLine), &
              stat=astatus)
  endif
  if (sat%chan_flag(6) > 0) then 
     allocate(sat%ws_albedo164(sat%nElem,sat%nLine), &
              stat=astatus)
  endif
  if (sat%chan_flag(7) > 0) then 
     allocate(sat%ws_albedo213(sat%nElem,sat%nLine), &
              stat=astatus)
  endif
           
  if (astatus /= 0) then
     print *,"(a,'Not enough memory to allocate ws albedo data structure.')"
     stop
  endif

! 3. end subroutine 
end subroutine  allocate_alb_arrays
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~ subroutine  11: deallocate alb data ~~~~~~~~~~~~~~~~
subroutine  deallocate_alb_arrays()
 
!-----------------------------------------------------------------------
! !F90 deallocate_alb_data
!
! !Description:
!    This program is to deallocate alb data.
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
 ! deallocate(sat%ws_albedo2,   &
 !            sat%ws_albedo3,   &
 !            sat%ws_albedo4,   &
 !            sat%ws_albedo5,   &
 !            sat%ws_albedo6,   &
 !            stat=astatus)

  if (sat%chan_flag(3) > 0) then 
     deallocate(sat%ws_albedo66, &
                stat=astatus)
  endif
  if (sat%chan_flag(4) > 0) then 
     deallocate(sat%ws_albedo87, &
                stat=astatus)
  endif
  if (sat%chan_flag(5) > 0) then 
     deallocate(sat%ws_albedo124, &
                stat=astatus)
  endif
  if (sat%chan_flag(6) > 0) then 
     deallocate(sat%ws_albedo164, &
                stat=astatus)
  endif
  if (sat%chan_flag(7) > 0) then 
     deallocate(sat%ws_albedo213, &
                stat=astatus)
  endif  
             
  if (astatus /= 0) then
     print *,"(a,'Error deallocate ws albedo data structure.')"
     stop
  endif
 
! 3. end subroutine 
end subroutine  deallocate_alb_arrays
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

!~~~~~~~~~~~~~~~~~~~ subroutine 8:  deal_AlbMap_dayid ~~~~~~~~~~~~~~~~~~
subroutine deal_alb_dayid(day, day_alb_id)

!-----------------------------------------------------------------------
! !F90 deal_alb_dayid
!
! !Description:
!    This is a main program to deal AlbMap daytime.
!
! !Input  parameters:
!    day             =  day for corresponding day data
!          
! !Output parameters:
!    none
!
!-----------------------------------------------------------------------

! 1. define variables
  integer(kind=4), intent(in) :: day
  integer(kind=4), intent(out):: day_alb_id 
  integer(kind=4) :: dif, i

! 2. begin program

do i=1,23

   dif = day-Alb_day(i)
   
   if (dif >= 0) then
       day_alb_id = i
       cycle
   endif
   
enddo

! 3. end subroutine     
end subroutine deal_alb_dayid
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~ subroutine 12:  open emiss file ~~~~~~~~~~~~~~~~~~~
subroutine open_alb_file(day_id, id2, id3, id4, id5, id6)

!-----------------------------------------------------------------------
! !F90 open_alb_file
!
! !Description:
!    This program is to open the surface emissivity file.
!
! !Input  parameters:
!    month           =  month for corresponding emissivity data
!          
! !Output parameters:
!    none
!
!-----------------------------------------------------------------------

! 1. define variables
  integer(kind=4), intent(in) :: day_id
  integer(kind=4), intent(out) :: id2, id3, id4, id5, id6
  
  character(len=200) :: data_dir
  character(len=256) :: filename2, filename3, filename4, filename5, filename6
  character(len=4),parameter,dimension(23) :: jday_str = (/    &
       "001", "017", "033", "049", "065", "081", "097", "113", &
       "129", "145", "161", "177", "193", "209", "225", "241", &
       "257", "273", "289", "305", "321", "337", "353"/)
       
  character(len=5)   :: year_str
  
  logical :: file_exists
  
! 2. begin program  
  year_str = "2004."
  
  data_dir = trim(code_root_path)//"/coeff/sfc_albedo"
  !filename2 = trim(data_dir)//"/AlbMap.WS.c004.v2.0."//trim(year_str)//trim(jday_str(day_id))//".0.659_x4.hdf"
  !filename5 = trim(data_dir)//"/AlbMap.WS.c004.v2.0."//trim(year_str)//trim(jday_str(day_id))//".1.64_x4.hdf"
  filename2 = trim(data_dir)//"/AlbMap.WS.c004.v2.0.00-04."//trim(jday_str(day_id))//".0.659_x4.hdf"
  filename3 = trim(data_dir)//"/AlbMap.WS.c004.v2.0.00-04."//trim(jday_str(day_id))//".0.858_x4.hdf"
  filename4 = trim(data_dir)//"/AlbMap.WS.c004.v2.0.00-04."//trim(jday_str(day_id))//".1.24_x4.hdf"
  filename5 = trim(data_dir)//"/AlbMap.WS.c004.v2.0.00-04."//trim(jday_str(day_id))//".1.64_x4.hdf"   
  filename6 = trim(data_dir)//"/AlbMap.WS.c004.v2.0.00-04."//trim(jday_str(day_id))//".2.13_x4.hdf"   

if (sat%chan_flag(3) > 0) then     
  inquire(file = filename2, exist = file_exists)
  if (.not. file_exists) then
     print *,"(/,a,'Surface alb 0.66 file, ',a,' does not exist.')"
     stop
  endif
  
  id2 = sfstart(trim(filename2), DFACC_READ)
  if (id2 == FAIL) then
     print *,"(/,a,'Failed to open, ',a)",trim(filename2)
     stop
  endif
endif

if (sat%chan_flag(4) > 0) then 
  inquire(file = filename3, exist = file_exists)
  if (.not. file_exists) then
     print *,"(/,a,'Surface alb 0.858 file, ',a,' does not exist.')"
     stop
  endif
  
  id3 = sfstart(trim(filename3), DFACC_READ)
  if (id3 == FAIL) then
     print *,"(/,a,'Failed to open, ',a)",trim(filename3)
     stop
  endif
endif

if (sat%chan_flag(5) > 0) then 
  inquire(file = filename4, exist = file_exists)
  if (.not. file_exists) then
     print *,"(/,a,'Surface alb 1.24 file, ',a,' does not exist.')"
     stop
  endif
  
  id4 = sfstart(trim(filename4), DFACC_READ)
  if (id4 == FAIL) then
     print *,"(/,a,'Failed to open, ',a)",trim(filename4)
     stop
  endif
endif

if (sat%chan_flag(6) > 0) then 
  inquire(file = filename5, exist = file_exists)
  if (.not. file_exists) then
     print *,"(/,a,'Surface alb 1.64 file, ',a,' does not exist.')"
     stop
  endif
  
  id5 = sfstart(trim(filename5), DFACC_READ)
  if (id5 == FAIL) then
     print *,"(/,a,'Failed to open, ',a)",trim(filename5)
     stop
  endif
endif

if (sat%chan_flag(7) > 0) then 
  inquire(file = filename6, exist = file_exists)
  if (.not. file_exists) then
     print *,"(/,a,'Surface alb 2.13 file, ',a,' does not exist.')"
     stop
  endif
  
  id6 = sfstart(trim(filename6), DFACC_READ)
  if (id6 == FAIL) then
     print *,"(/,a,'Failed to open, ',a)",trim(filename6)
     stop
  endif
endif
    
! 3. end subroutine     
end subroutine open_alb_file
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

!~~~~~~~~~~~~~~~~~~~ subroutine 13:  read alb file ~~~~~~~~~~~~~~~~~~~~~
subroutine read_alb_file(id2, id3, id4, id5, id6)

!-----------------------------------------------------------------------
! !F90 read_alb_file
!
! !Description:
!    This program is to read the surface albedo file.
!
! !Input  parameters:
!    id1/2           =  hdf file id
!          
! !Output parameters:
!    emiss           =  emissivity data
!
!-----------------------------------------------------------------------

! 1. define variables
  integer(kind=4), intent(in)                   :: id2, id3, id4, id5, id6
  real(kind=4),allocatable, dimension(:,:)     :: lat, lon
!  integer(kind=1), allocatable, dimension(:,:)  :: space_mask
  real(kind=4), allocatable, dimension(:,:)    :: alb2, alb3, alb4, alb5, alb6
!  real(kind=4), dimension(num_lon_alb,num_lat_alb) :: alb_org2, alb_org3, alb_org4, alb_org5, alb_org6
  real(kind=4), allocatable, dimension(:,:) :: alb_org2, alb_org3, alb_org4, alb_org5, alb_org6
  integer(kind=4) ::  nx, ny, xr, yr, i, j
  
  integer(kind=4) :: astatus

  character(len=16) :: sds_name2
  character(len=16) :: sds_name3
  character(len=15) :: sds_name4
  character(len=15) :: sds_name5
  character(len=15) :: sds_name6
  
! 2. begin program
xr = sat%nElem
yr = sat%nLine

!=== 2.1. allocate variables
allocate(lat(xr,yr),lon(xr,yr),       &
         alb2(xr,yr), alb3(xr,yr), alb4(xr,yr), alb5(xr,yr), alb6(xr,yr),      &
         alb_org2(num_lon_alb,num_lat_alb), alb_org3(num_lon_alb,num_lat_alb), &
         alb_org4(num_lon_alb,num_lat_alb), alb_org5(num_lon_alb,num_lat_alb), &
         alb_org6(num_lon_alb,num_lat_alb))

!=== 2.2. assign data
lon        = geo%lon
lat        = geo%lat
!space_mask = sat%space_mask
alb2       = missing_value_real4
alb3       = missing_value_real4
alb4       = missing_value_real4
alb5       = missing_value_real4
alb6       = missing_value_real4
nx         = missing_value_int4
ny         = missing_value_int4

!=== 2.3. read start
! read channel 2 and 5

  sds_name2 = 'Albedo_Map_0.659'
  sds_name3 = 'Albedo_Map_0.858'
  sds_name4 = 'Albedo_Map_1.24'
  sds_name5 = 'Albedo_Map_1.64'  
  sds_name6 = 'Albedo_Map_2.13'

  if (sat%chan_flag(3) > 0) then 
  	  call read_alb_hdf(id2, trim(sds_name2), alb_org2) 
  endif
  if (sat%chan_flag(4) > 0) then 
      call read_alb_hdf(id3, trim(sds_name3), alb_org3) 
  endif
  if (sat%chan_flag(5) > 0) then 
      call read_alb_hdf(id4, trim(sds_name4), alb_org4) 
  endif
  if (sat%chan_flag(6) > 0) then 
      call read_alb_hdf(id5, trim(sds_name5), alb_org5) 
  endif
  if (sat%chan_flag(7) > 0) then 
      call read_alb_hdf(id6, trim(sds_name6), alb_org6) 
  endif
    
  ! match lon lat 
  do i = 1, xr
  do j = 1, yr
  
     !if (space_mask(i,j) /= sym%SPACE ) then
     
     nx = INT((lon(i,j) - first_lon_alb)/del_lon_alb)+1
     ny = INT((first_lat_alb - lat(i,j))/del_lat_alb)+1
     if (nx < 1) then 
     	nx = 1
     endif
     if (nx > num_lon_alb) then 
     	nx = num_lon_alb
     endif
     if (ny < 1) then 
        ny = 1
     endif
     if (ny > num_lat_alb) then 
        ny = num_lat_alb   
     endif
        
     if (sat%chan_flag(3) > 0) then 
        if (alb_org2(nx,ny) > 100. .or. alb_org2(nx,ny) <= 0.) then
           alb_org2(nx,ny) = Missing_Value_real4
        endif     
        alb2(i,j) = alb_org2(nx,ny)
     endif
     
     if (sat%chan_flag(4) > 0) then 
        if (alb_org3(nx,ny) > 100. .or. alb_org3(nx,ny) <= 0.) then
           alb_org3(nx,ny) = Missing_Value_real4
        endif     
        alb3(i,j) = alb_org3(nx,ny)
     endif

     if (sat%chan_flag(5) > 0) then 
        if (alb_org4(nx,ny) > 100. .or. alb_org4(nx,ny) <= 0.) then
           alb_org4(nx,ny) = Missing_Value_real4
        endif     
        alb4(i,j) = alb_org4(nx,ny)
     endif
    
     if (sat%chan_flag(6) > 0) then            
        if (alb_org5(nx,ny) > 100. .or. alb_org5(nx,ny) <= 0.) then
           alb_org5(nx,ny) = Missing_Value_real4
        endif
        alb5(i,j) = alb_org5(nx,ny)
     endif

     if (sat%chan_flag(7) > 0) then 
        if (alb_org6(nx,ny) > 100. .or. alb_org6(nx,ny) <= 0.) then
           alb_org6(nx,ny) = Missing_Value_real4
        endif
        alb6(i,j) = alb_org6(nx,ny)
     endif
      
    !endif
     
  enddo
  enddo

if (sat%chan_flag(3) > 0) then 
  sat%ws_albedo66   = alb2
endif
if (sat%chan_flag(4) > 0) then 
  sat%ws_albedo87   = alb3
endif
if (sat%chan_flag(5) > 0) then 
  sat%ws_albedo124  = alb4
endif
if (sat%chan_flag(6) > 0) then 
  sat%ws_albedo164  = alb5
endif
if (sat%chan_flag(7) > 0) then 
  sat%ws_albedo213  = alb6
endif
  
!=== 2.4. deallocate variables
  deallocate(lat,lon,alb2,alb3,alb4,alb5,alb6, &
             alb_org2, alb_org3, alb_org4, alb_org5, alb_org6, stat=astatus)
             
  if (astatus /= 0) then
     print *,"(a,'Error deallocating ws albedo data file.')"
     stop
  endif
  
! 3. end subroutine  
end subroutine read_alb_file
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

!~~~~~~~~~~~~~~~~~~~ subroutine 14:  read alb hdf file ~~~~~~~~~~~~~~~~~
subroutine read_alb_hdf(sd_id, sds_name, alb_org)

!-----------------------------------------------------------------------
! !F90 read_alb_hdf
!
! !Description:
!    This program is to read the surface albedo hdf file.
!
! !Input  parameters:
!    sd_id           =  hdf file id
!    sds_name        =  satellite channel number name
!          
! !Output parameters:
!    alb_org           =  albedo data
!
!-----------------------------------------------------------------------

! 1. define variables
!===== 1.1.hdf4 Function declaration.
character(*), intent(in) :: sds_name

integer(kind=4) sd_id, sds_id, attr_id
integer(kind=4), dimension(2) :: start, stride, edges

integer(kind=2),allocatable, dimension(:,:) :: alb_org1
!real(kind=4), dimension(num_lon_alb,num_lat_alb) :: alb_org
real(kind=4), dimension(:,:),  intent(out) :: alb_org

!===== 1.2. other
integer(kind=4) :: status  
character*100 message

!*******
! 2. begin program
  !===== 2.1. initialize
  allocate(alb_org1(num_lon_alb,num_lat_alb))
  start    = 0
  stride   = 1
  edges(1) = num_lon_alb
  edges(2) = num_lat_alb ! 16 column * 1 row
  status   = 0
  message = 'WS Surface albedo'

  ! read albedo
  sds_id = sfselect(sd_id, sfn2index(sd_id,sds_name))   
  status = sfrdata (sds_id, start, stride, edges, alb_org1)
  call hdf_info_message(message, sds_name, status) 
  
  status = sfendacc(sds_id)
  
  alb_org = alb_org1*0.1

  deallocate(alb_org1)

! 3. end subroutine    
end subroutine read_alb_hdf
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

!~~~~~~~~~~~~~~~~~~~ subroutine 15:  close alb hdf file ~~~~~~~~~~~~~~~~
subroutine close_alb_file(id2, id3, id4, id5, id6)

!-----------------------------------------------------------------------
! !F90 close_alb_file
!
! !Description:
!    This program is to close the surface albedo hdf file.
!
! !Input  parameters:
!    id1/2           =  hdf file id
!          
! !Output parameters:
!    none
!
!-----------------------------------------------------------------------

  integer(kind=4), intent(in) :: id2, id3, id4, id5, id6
  
  integer(kind=4) :: istatus
  integer(kind=4) :: sfend

if (sat%chan_flag(3) > 0) then 
  istatus = sfend(id2)
  if (istatus /= 0) then
     print *,"(/,a,'Error closing surface ws albedo 0.66 hdf file.')"
     stop
  endif
endif

if (sat%chan_flag(4) > 0) then 
  istatus = sfend(id3)
  if (istatus /= 0) then
     print *,"(/,a,'Error closing surface ws albedo 0.858 hdf file.')"
     stop
  endif
endif

if (sat%chan_flag(5) > 0) then   
  istatus = sfend(id4)
  if (istatus /= 0) then
     print *,"(/,a,'Error closing surface ws albedo 1.24 hdf file.')"
     stop
  endif
endif

if (sat%chan_flag(6) > 0) then    
  istatus = sfend(id5)
  if (istatus /= 0) then
     print *,"(/,a,'Error closing surface ws albedo 1.64 hdf file.')"
     stop
  endif
endif

if (sat%chan_flag(7) > 0) then 
  istatus = sfend(id6)
  if (istatus /= 0) then
     print *,"(/,a,'Error closing surface ws albedo 2.13 hdf file.')"
     stop
  endif
endif

end subroutine close_alb_file
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~ subroutine  20:  close smk hdf file ~~~~~~~~~~~~~~~~
!subroutine  read_ecosystem_file(rlon,rlat,npixel,scans_cube,cube_eco)

!-----------------------------------------------------------------------
! !F90 read_ecosystem_file
!
! !description:
!    This program is to close the surface albedo hdf file.
!
! !Input  parameters:
!    npixel     = sat%nElem
!    scans_cube = sat%nLine/10
!          
! !Output parameters:
!    cube_eco
!
!-----------------------------------------------------------------------

!  integer(kind=4), intent(in) :: npixel,scans_cube
!  integer(kind=4)  :: eco_1km_lun     &
!                      eco_10min_lun
  
!  integer(kind=4) :: istatus, io_err, iocheck.
!       local scalars
!        real b_lat,b_lon,land,prmsl,pwat,ugrd,vgrd,ozone,icec,     &
!             sst,landtmp,lat_temp,lon_temp
!        integer lat_indx,lon_indx,bytloc,ii,jj,mm,nn,h_output,debug,   &
!               i,j,newlin_eco2,iocheck,newele_eco2,nise,nise_out
!        logical init
!       scalar arguements
!        integer cube,max_pixel 
!        integer buf_anc_size(2) 
!       local arrays
!        real plat(1,1),plon(1,1)
!        byte map_eco(2),madison_eco,ecotest(1,1)
!        byte cube_eco(npixel,scans_cube),qa_bits(10)
!        real rlat(npixel,scans_cube)
!        real rlon(npixel,scans_cube)
    
!  data init / .true. /
  
! --- initialize the lun
!eco_1km_lun = -5555
!eco_10min_lun = -5555
!cube_eco = 0
 !        do 134 mm = 1 , buf_anc_size(2)
 !        do 135 nn = 1 , buf_anc_size(1)
 !           cube_eco(nn,mm) = 0
 ! 135    continue
 ! 134    continue
! ... set max buf size of ancillary data file
!buf_anc_size(1) = npixel
!buf_anc_size(2) = scans_cube

! GLOBAL 1KM ECOSYSTEM FILE     
! ------ issue the OPEN request
!open( FILE=trim(code_root_path)//'coeff/'//gogel_name,    &
!      UNIT=ECO_1KM_UNIT,      &
!      STATUS=ECO_1KM_STATUS,  &
!      FORM=ECO_1KM_FORM,      & 
!      ACCESS=ECO_1KM_ACCESS,  &
!      RECL=ECO_1KM_RECL,      &
!      IOSTAT=io_err) 
!      eco_1km_lun = ECO_1KM_UNIT
      
!----------------------------------------------------------------------
!        Okay, now read the ecosystem information

!         if (init) then
! ...      First, let's make sure we are reading from the correct
! ...       file.  Extract Madison Wisconsin eco_type from file
!           if ( eco_1km_lun .ne. -5555) then
!              plat(1,1) = 43.083
!              plon(1,1) = -89.305
!              call read_goode( 1,1,plat,plon,eco_1km_lun,ecotest )
       !      if ( ecotest(1,1) .ne. 22 ) call message( 'get_anc_data',      &
       !        'Extracted incorrect ecosystem value from 1 km file. ' //    &
       !        '[OPERATOR ACTION: Verify size and format of file. ' //      &
       !        'if error persists, contact SDST.]', 0, 2 )
!           else
! ...         Read from the 10 minute global Olson ecosystem file
!              b_lat = 43.083
!              b_lon = -89.305
!              lat_indx = int((90.0 - b_lat) * 6.0 + 1.0)
!              lon_indx = int((b_lon + 180.0) * 6.0 + 1.0)
!              if(lat_indx .gt. 1080) lat_indx = 1080
!              if(lon_indx .gt. 2160) lon_indx = 2160
!              bytloc = ((lat_indx-1) * 2160) + lon_indx
!              newlin_eco2 = bytloc / 2 + 1
!              newele_eco2 = mod(bytloc, 2)
!              if (mod( bytloc,2) .eq. 0) then
!                 newlin_eco2 = newlin_eco2 - 1
!                 newele_eco2 = 2
!              endif

! ...        read value from 10minute file
!             read (eco_10min_lun, rec=newlin_eco2, iostat=iocheck) map_eco
!             if (iocheck .ne. 0) then
     !          call message( 'get_anc_data',                       &
     !         'Error reading ecosystem value from file. ' //       &
     !         'Make sure correct 10 minute file is loaded. ' //    &
     !         ' [OPERATOR ACTION: Notify SDST.]',                  &
     !         0, 2 )
!             endif

!             madison_eco = map_eco(newele_eco2)
! ...        Compare value extracted with known correct value
!             if (madison_eco .ne. 55) then
     !          call message( 'get_anc_data',                            &
     !         'Extracted incorrect ecosystem value from file. ' //      &
     !         'Make sure correct 10 minute file is loaded. ' //         &
     !         ' [OPERATOR ACTION: Notify SDST.]',                       &
     !         0, 2 )
!             endif
!           endif

!           init = .false.
!         endif
! -------------------------------------------------------------------

!         if (eco_1km_lun .ne. -5555) then
! ...      Read all ecosystem values for this scan cube

!           call read_goode( buf_anc_size( 1 ), buf_anc_size( 2 ), rlat, rlon, eco_1km_lun, cube_eco )
           
!         else

! ...      Read out of 10 minute global file
!           do  ii = 1 , scans_cube
!             do  jj = 1 , npixel
             
!                b_lat = rlat(jj,ii)
!                b_lon = rlon(jj,ii)
!                lat_indx = int((90.0 - b_lat) * 6.0 + 1.0)
!                lon_indx = int((b_lon + 180.0) * 6.0 + 1.0)
!                if(lat_indx .gt. 1080) lat_indx = 1080
!                if(lon_indx .gt. 2160) lon_indx = 2160
!                bytloc = ((lat_indx-1) * 2160) + lon_indx
!                newlin_eco2 = bytloc / 2 + 1
!                newele_eco2 = mod(bytloc, 2)
!                if (mod( bytloc,2) .eq. 0) then
!                  newlin_eco2 = newlin_eco2 - 1
!                  newele_eco2 = 2
!                endif

! ...           read value from file
!                read (eco_10min_lun, rec=newlin_eco2, iostat=iocheck) map_eco
!                if (iocheck .ne. 0) then
     !             call message( 'get_anc_data',                    &
     !            'Error reading ecosystem value from file. ' //    &
     !            ' [OPERATOR ACTION: Notify SDST.]',               &
     !            0, 2 )                                            &
!                endif
!                cube_eco(jj,ii) = map_eco(newele_eco2)
!                
!              enddo
!           enddo

!         endif 

!         if (eco_1km_lun .ne. -5555) then
! ...      Read all ecosystem values for this scan cube

!             close(ECO_1KM_UNIT)
!         else             


!end subroutine  read_ecosystem_file
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!-------------------------- end MODULE ---------------------------------
end module get_ancil_data_module