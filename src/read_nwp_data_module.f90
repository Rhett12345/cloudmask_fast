module read_nwp_data_module

!C-----------------------------------------------------------------------
!C !F90                                                                  
!C
!C !description: 
!C    This MODULE is to read Numerical Weather Prediction [nwp] DATA arrays
!C    for fylat platform.
!C
!C !Input parameters
!C    none
!C 
!C !Output parameters
!C    none
!C
!C  
!C !end
!C----------------------------------------------------------------------

! use modules
use data_arrays_module
use message_module
use constant
use numerical
use names_module 

implicit none
!+++++++++++++++++++ step 1: define global variables +++++++++++++++++++
! |------|
! | none |
! |------|
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


CONTAINS
!+++++++++++++++++++ step 2: subroutines +++++++++++++++++++++++++++++++

!~~~~~~~~~~~~~~~~~~~ subroutine 1: read nwp DATA main program ~~~~~~~~~~
subroutine fylat_read_nwp_data()

!-----------------------------------------------------------------------
! !F90 read_nwp_data_main
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
character*300 bin_name1
character*300 bin_name2


! 2. begin program
  print*,'  ... read nwp data'
  
  !===== 2.0. extract nwp time
  call extract_nwptime(nwp_data_path, nwp_grib_data1, 1,     &
                       nwptime%year(1), nwptime%month(1),    &
                       nwptime%day(1),  nwptime%hour(1))
  call extract_nwptime(nwp_data_path, nwp_grib_data2, 2,     &
                       nwptime%year(2), nwptime%month(2),    &
                       nwptime%day(2),  nwptime%hour(2))
  
  !===== 2.1. allocate nwpo variables
  call allocate_nwpo_arrays(fylat_nwp_opt)

  !===== 2.2. convert grib to binary for nwp DATA
  call convert_grib_to_binary(nwp_grib_data1, bin_name1, 1)
  call convert_grib_to_binary(nwp_grib_data2, bin_name2, 2)

  
  !===== 2.3. read nwp arrays in binary format
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
     call read_nwp_arrays(bin_name1, 1) 
     call read_nwp_arrays(bin_name2, 2) 
  endif
  
  if (fylat_nwp_opt == 3) then   ! 3 = T639
     call read_nwp_T639_arrays(bin_name1, 1) 
     call read_nwp_T639_arrays(bin_name2, 2)   
  endif

  IF (fylat_nwp_opt == 6) THEN   ! 6=grapes gfs 0.25*0.25 (1440*720)
     print*,'read nwp grapes gfs binary data'
     call read_NWP_grapes_gfs_arrays(bin_name1,1) 
     call read_NWP_grapes_gfs_arrays(bin_name2,2)   
  ENDIF
  
  IF (fylat_nwp_opt == 7 .or. fylat_nwp_opt == 8) THEN   ! 8=gfs 0p25 with 31 layers
     print*,'read nwp binary data with 31 layers'
     CALL read_NWP_arrays_0p25(bin_name1,1) 
     CALL read_NWP_arrays_0p25(bin_name2,2) 
  ENDIF

  IF (fylat_nwp_opt == 9) THEN   ! 9=gfs 0.5*0.5 with 41 layers (720*361)
     print*,'read nwp gfs 0p50 binary data with 41 layers '
     CALL read_nwp_arrays_0p50_41Layers(bin_name1,1) 
     CALL read_nwp_arrays_0p50_41Layers(bin_name2,2)   
  ENDIF

  IF (fylat_nwp_opt == 10) THEN   ! 10=gfs 0p25 with 41 layers (1440*721)
     print*,'read nwp gfs 0p25 binary data with 41 layers '
     CALL read_NWP_arrays_0p25_41Layers(bin_name1,1) 
     CALL read_NWP_arrays_0p25_41Layers(bin_name2,2) 
  ENDIF

  !===== 2.4. find x, y nwp position
  call find_nwp_eqarea_cell(geo%lon,geo%lat,sat%x_nwp,sat%y_nwp,nwp%first_lat)
  !call find_nwp_eqarea_cell(geo%lon,geo%lat,sat%x2_nwp,sat%y2_nwp,nwp%first_lat2) 
   
! 3. end subroutine   
end subroutine fylat_read_nwp_data
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~ subroutine 2: convert grib to binary ~~~~~~~~~~~~~~
subroutine convert_grib_to_binary(grib_name, &
                                  bin_name,  &
                                  id )

!-----------------------------------------------------------------------
! !F90 read_nwp_arrays
!
! !description:
!    This program is to read Navigation arrays.
!
! !Input  parameters:
!    year            = year
!    month           = month
!    day             = day
!    hour            = hour 
!
! !Output parameters:
!    bin_name         = binary DATA's name
!
!-----------------------------------------------------------------------

! 1. define variables
!===== 1.1.input variables
character*100 filename, xname  ! nwp grib DATA's name 
character*100 filename_T639 ! nwp T639 DATA's name
character*100 filename_gfs1p00 ! nwp gfs 1p00 DATA's name
character stime*2  ! nwp gfs prediction start time
integer(kind=1) :: id

!===== 1.2. nwp grib and binary file's name
character*300 grib_name  !full name with path directory
character*300 bin_name   !full name with path directory 

!===== 1.3. other
integer(kind=4) :: status, L1, L2, i, j, k
character*300 wgrib_exe, T639_exe, wgrib2_exe, wgrib2_grapes_gfs_exe, wgrib3_exe, wgrib4_exe, wgrib5_exe 

!===== 1.4. grib DATA out name
character(len=30) :: list(113+52)
character(len=30) :: list2(124)  ! grib2
! we define the variables follow sample hdf ncep DATA [gblav.xxx]
!                                  ! (1) pressure 10~1000 mb
DATA  list /'PRES:sfc',          & ! (2) surface pressure                  1
            'PRMSL:MSL',         & ! (3) MSL pressure                      2
            'TMP:sfc',           & ! (4) surface temperature               3
            'HGT:sfc',           & ! (5) surface height                    4
            'PRES:sfc',             & ! xxxxxx(6) surface albedo  ALBdo                  5
            'TMP:sigma=0.9950',  & ! (7) temperature at sigma=0.995        6
            'RH:sigma=0.9950',   & ! (8) rh at sigma=0.995                 7
            'UGRD:sigma=0.9950', & ! (9) U-wind at sigma=0.995             8
            'VGRD:sigma=0.9950', & ! (10)V-wind at sigma=0.995             9
            'PWAT',              & ! (11)total precipitable water          10
            'WEASD',             & ! (12)water equivalent snow depth       11
            'TOZNE',             & ! (13)total ozone                       12
            'TMP:tropopause',    & ! (14)tropopause temperature            13
            'TMP:10 mb',         & ! (15)temperature profile [26 layers]   14
            'TMP:20 mb',         & !        ...                            15
            'TMP:30 mb',         & !        ...                            16
            'TMP:50 mb',         & !        ...                            17
            'TMP:70 mb',         & !        ...                            18
            'TMP:100 mb',        & !        ...                            19
            'TMP:150 mb',        & !        ...                            20 
            'TMP:200 mb',        & !        ...                            21
            'TMP:250 mb',        & !        ...                            22
            'TMP:300 mb',        & !        ...                            23
            'TMP:350 mb',        & !        ...                            24
            'TMP:400 mb',        & !        ...                            25
            'TMP:450 mb',        & !        ...                            26
            'TMP:500 mb',        & !        ...                            27
            'TMP:550 mb',        & !        ...                            28
            'TMP:600 mb',        & !        ...                            29
            'TMP:650 mb',        & !        ...                            30
            'TMP:700 mb',        & !        ...                            31
            'TMP:750 mb',        & !        ...                            32  
            'TMP:800 mb',        & !        ...                            33
            'TMP:850 mb',        & !        ...                            34
            'TMP:900 mb',        & !        ...                            35
            'TMP:925 mb',        & !        ...                            36 
            'TMP:950 mb',        & !        ...                            37
            'TMP:975 mb',        & !        ...                            38
            'TMP:1000 mb',       & !        ...                            39
            'HGT:10 mb',         & ! (16)height profile [26 layers]        40
            'HGT:20 mb',         & !        ...                            41
            'HGT:30 mb',         & !        ...                            42
            'HGT:50 mb',         & !        ...                            43
            'HGT:70 mb',         & !        ...                            44
            'HGT:100 mb',        & !        ...                            45
            'HGT:150 mb',        & !        ...                            46 
            'HGT:200 mb',        & !        ...                            47
            'HGT:250 mb',        & !        ...                            48
            'HGT:300 mb',        & !        ...                            49
            'HGT:350 mb',        & !        ...                            50
            'HGT:400 mb',        & !        ...                            51
            'HGT:450 mb',        & !        ...                            52
            'HGT:500 mb',        & !        ...                            53
            'HGT:550 mb',        & !        ...                            54
            'HGT:600 mb',        & !        ...                            55
            'HGT:650 mb',        & !        ...                            56
            'HGT:700 mb',        & !        ...                            57
            'HGT:750 mb',        & !        ...                            58  
            'HGT:800 mb',        & !        ...                            59
            'HGT:850 mb',        & !        ...                            60
            'HGT:900 mb',        & !        ...                            61
            'HGT:925 mb',        & !        ...                            62 
            'HGT:950 mb',        & !        ...                            63
            'HGT:975 mb',        & !        ...                            64
            'HGT:1000 mb',       & !        ...                            65 
            'O3MR:10 mb',        & ! (17)O3 profile [26 layers]            66
            'O3MR:20 mb',        & !  only 6 layers are available          67
            'O3MR:30 mb',        & !        ...                            68
            'O3MR:50 mb',        & !        ...                            69
            'O3MR:70 mb',        & !        ...                            70
            'O3MR:100 mb',       & !        ...                            71
            'RH:100 mb',         & ! (18)RH profile [26 layers]            72
            'RH:150 mb',         & !  only 21 layers are available         73 
            'RH:200 mb',         & !        ...                            74
            'RH:250 mb',         & !        ...                            75
            'RH:300 mb',         & !        ...                            76
            'RH:350 mb',         & !        ...                            77
            'RH:400 mb',         & !        ...                            78
            'RH:450 mb',         & !        ...                            79
            'RH:500 mb',         & !        ...                            80
            'RH:550 mb',         & !        ...                            81
            'RH:600 mb',         & !        ...                            82
            'RH:650 mb',         & !        ...                            83
            'RH:700 mb',         & !        ...                            84
            'RH:750 mb',         & !        ...                            85  
            'RH:800 mb',         & !        ...                            86
            'RH:850 mb',         & !        ...                            87
            'RH:900 mb',         & !        ...                            88
            'RH:925 mb',         & !        ...                            89 
            'RH:950 mb',         & !        ...                            90
            'RH:975 mb',         & !        ...                            91
            'RH:1000 mb',        & !        ...                            92 
            'CLWMR:100 mb',      & ! (19)Cloud water profile [26 layers]   93
            'CLWMR:150 mb',      & !  only 21 layers are available         94 
            'CLWMR:200 mb',      & !        ...                            95
            'CLWMR:250 mb',      & !        ...                            96
            'CLWMR:300 mb',      & !        ...                            97
            'CLWMR:350 mb',      & !        ...                            98
            'CLWMR:400 mb',      & !        ...                            99
            'CLWMR:450 mb',      & !        ...                            100
            'CLWMR:500 mb',      & !        ...                            101
            'CLWMR:550 mb',      & !        ...                            102
            'CLWMR:600 mb',      & !        ...                            103
            'CLWMR:650 mb',      & !        ...                            104
            'CLWMR:700 mb',      & !        ...                            105
            'CLWMR:750 mb',      & !        ...                            106 
            'CLWMR:800 mb',      & !        ...                            107
            'CLWMR:850 mb',      & !        ...                            108
            'CLWMR:900 mb',      & !        ...                            109
            'CLWMR:925 mb',      & !        ...                            110
            'CLWMR:950 mb',      & !        ...                            111
            'CLWMR:975 mb',      & !        ...                            112
            'CLWMR:1000 mb',     & !        ...                            113
            'UGRD:10 mb',         & ! (16)height profile [26 layers]        40
            'UGRD:20 mb',         & !        ...                            41
            'UGRD:30 mb',         & !        ...                            42
            'UGRD:50 mb',         & !        ...                            43
            'UGRD:70 mb',         & !        ...                            44
            'UGRD:100 mb',        & !        ...                            45
            'UGRD:150 mb',        & !        ...                            46 
            'UGRD:200 mb',        & !        ...                            47
            'UGRD:250 mb',        & !        ...                            48
            'UGRD:300 mb',        & !        ...                            49
            'UGRD:350 mb',        & !        ...                            50
            'UGRD:400 mb',        & !        ...                            51
            'UGRD:450 mb',        & !        ...                            52
            'UGRD:500 mb',        & !        ...                            53
            'UGRD:550 mb',        & !        ...                            54
            'UGRD:600 mb',        & !        ...                            55
            'UGRD:650 mb',        & !        ...                            56
            'UGRD:700 mb',        & !        ...                            57
            'UGRD:750 mb',        & !        ...                            58  
            'UGRD:800 mb',        & !        ...                            59
            'UGRD:850 mb',        & !        ...                            60
            'UGRD:900 mb',        & !        ...                            61
            'UGRD:925 mb',        & !        ...                            62 
            'UGRD:950 mb',        & !        ...                            63
            'UGRD:975 mb',        & !        ...                            64
            'UGRD:1000 mb',       & !        ...                            65 
            'VGRD:10 mb',         & ! (16)height profile [26 layers]        40
            'VGRD:20 mb',         & !        ...                            41
            'VGRD:30 mb',         & !        ...                            42
            'VGRD:50 mb',         & !        ...                            43
            'VGRD:70 mb',         & !        ...                            44
            'VGRD:100 mb',        & !        ...                            45
            'VGRD:150 mb',        & !        ...                            46 
            'VGRD:200 mb',        & !        ...                            47
            'VGRD:250 mb',        & !        ...                            48
            'VGRD:300 mb',        & !        ...                            49
            'VGRD:350 mb',        & !        ...                            50
            'VGRD:400 mb',        & !        ...                            51
            'VGRD:450 mb',        & !        ...                            52
            'VGRD:500 mb',        & !        ...                            53
            'VGRD:550 mb',        & !        ...                            54
            'VGRD:600 mb',        & !        ...                            55
            'VGRD:650 mb',        & !        ...                            56
            'VGRD:700 mb',        & !        ...                            57
            'VGRD:750 mb',        & !        ...                            58  
            'VGRD:800 mb',        & !        ...                            59
            'VGRD:850 mb',        & !        ...                            60
            'VGRD:900 mb',        & !        ...                            61
            'VGRD:925 mb',        & !        ...                            62 
            'VGRD:950 mb',        & !        ...                            63
            'VGRD:975 mb',        & !        ...                            64
            'VGRD:1000 mb'/!        ...      

DATA  list2/' "PRES:surface"',           & ! (2) surface pressure                  1
            ' "PRMSL:mean"',             & ! (3) MSL pressure                      2
            ' "TMP:surface"',            & ! (4) surface temperature               3
            ' "HGT:surface"',            & ! (5) surface height                    4
            ' "PRES:surface"',           & ! xxxxxx(6) surface albedo  ALBdo       5
            ' "TMP:0.995 sigma level"',  & ! (7) temperature at sigma=0.995        6
            ' "RH:0.995 sigma level"',   & ! (8) rh at sigma=0.995                 7
            ' "UGRD:0.995 sigma level"', & ! (9) U-wind at sigma=0.995             8
            ' "VGRD:0.995 sigma level"', & ! (10)V-wind at sigma=0.995             9
            ' "PWAT"',              & ! (11)total precipitable water          10
            ' "WEASD:surface"',             & ! (12)water equivalent snow depth       11
            ' "TOZNE"',             & ! (13)total ozone                       12
            ' "TMP:tropopause"',    & ! (14)tropopause temperature            13
            ' "TMP:10 mb"',         & ! (15)temperature profile [26 layers]   14
            ' "TMP:20 mb"',         & !        ...                            15
            ' "TMP:30 mb"',         & !        ...                            16
            ' "TMP:50 mb"',         & !        ...                            17
            ' "TMP:70 mb"',         & !        ...                            18
            ' "TMP:100 mb"',        & !        ...                            19
            ' "TMP:150 mb"',        & !        ...                            20 
            ' "TMP:200 mb"',        & !        ...                            21
            ' "TMP:250 mb"',        & !        ...                            22
            ' "TMP:300 mb"',        & !        ...                            23
            ' "TMP:350 mb"',        & !        ...                            24
            ' "TMP:400 mb"',        & !        ...                            25
            ' "TMP:450 mb"',        & !        ...                            26
            ' "TMP:500 mb"',        & !        ...                            27
            ' "TMP:550 mb"',        & !        ...                            28
            ' "TMP:600 mb"',        & !        ...                            29
            ' "TMP:650 mb"',        & !        ...                            30
            ' "TMP:700 mb"',        & !        ...                            31
            ' "TMP:750 mb"',        & !        ...                            32  
            ' "TMP:800 mb"',        & !        ...                            33
            ' "TMP:850 mb"',        & !        ...                            34
            ' "TMP:900 mb"',        & !        ...                            35
            ' "TMP:925 mb"',        & !        ...                            36 
            ' "TMP:950 mb"',        & !        ...                            37
            ' "TMP:975 mb"',        & !        ...                            38
            ' "TMP:1000 mb"',       & !        ...                            39
            ' "HGT:10 mb"',         & ! (16)height profile [26 layers]        40
            ' "HGT:20 mb"',         & !        ...                            41
            ' "HGT:30 mb"',         & !        ...                            42
            ' "HGT:50 mb"',         & !        ...                            43
            ' "HGT:70 mb"',         & !        ...                            44
            ' "HGT:100 mb"',        & !        ...                            45
            ' "HGT:150 mb"',        & !        ...                            46 
            ' "HGT:200 mb"',        & !        ...                            47
            ' "HGT:250 mb"',        & !        ...                            48
            ' "HGT:300 mb"',        & !        ...                            49
            ' "HGT:350 mb"',        & !        ...                            50
            ' "HGT:400 mb"',        & !        ...                            51
            ' "HGT:450 mb"',        & !        ...                            52
            ' "HGT:500 mb"',        & !        ...                            53
            ' "HGT:550 mb"',        & !        ...                            54
            ' "HGT:600 mb"',        & !        ...                            55
            ' "HGT:650 mb"',        & !        ...                            56
            ' "HGT:700 mb"',        & !        ...                            57
            ' "HGT:750 mb"',        & !        ...                            58  
            ' "HGT:800 mb"',        & !        ...                            59
            ' "HGT:850 mb"',        & !        ...                            60
            ' "HGT:900 mb"',        & !        ...                            61
            ' "HGT:925 mb"',        & !        ...                            62 
            ' "HGT:950 mb"',        & !        ...                            63
            ' "HGT:975 mb"',        & !        ...                            64
            ' "HGT:1000 mb"',       & !        ...                            65 
            ' "O3MR:10 mb"',        & ! (17)O3 profile [26 layers]            66
            ' "O3MR:20 mb"',        & !  only 6 layers are available          67
            ' "O3MR:30 mb"',        & !        ...                            68
            ' "O3MR:50 mb"',        & !        ...                            69
            ' "O3MR:70 mb"',        & !        ...                            70
            ' "O3MR:100 mb"',       & !        ...                            71
            ' "RH:100 mb"',         & ! (18)RH profile [26 layers]            72
            ' "RH:150 mb"',         & !  only 26 layers are available         73 
            ' "RH:200 mb"',         & !        ...                            74
            ' "RH:250 mb"',         & !        ...                            75
            ' "RH:300 mb"',         & !        ...                            76
            ' "RH:350 mb"',         & !        ...                            77
            ' "RH:400 mb"',         & !        ...                            78
            ' "RH:450 mb"',         & !        ...                            79
            ' "RH:500 mb"',         & !        ...                            80
            ' "RH:550 mb"',         & !        ...                            81
            ' "RH:600 mb"',         & !        ...                            82
            ' "RH:650 mb"',         & !        ...                            83
            ' "RH:700 mb"',         & !        ...                            84
            ' "RH:750 mb"',         & !        ...                            85  
            ' "RH:800 mb"',         & !        ...                            86
            ' "RH:850 mb"',         & !        ...                            87
            ' "RH:900 mb"',         & !        ...                            88
            ' "RH:925 mb"',         & !        ...                            89 
            ' "RH:950 mb"',         & !        ...                            90
            ' "RH:975 mb"',         & !        ...                            91
            ' "RH:1000 mb"',        & !        ...                            92 
            ' "CLWMR:100 mb"',      & ! (19)Cloud water profile [26 layers]   93
            ' "CLWMR:150 mb"',      & !  only 21 layers are available         94 
            ' "CLWMR:200 mb"',      & !        ...                            95
            ' "CLWMR:250 mb"',      & !        ...                            96
            ' "CLWMR:300 mb"',      & !        ...                            97
            ' "CLWMR:350 mb"',      & !        ...                            98
            ' "CLWMR:400 mb"',      & !        ...                            99
            ' "CLWMR:450 mb"',      & !        ...                            100
            ' "CLWMR:500 mb"',      & !        ...                            101
            ' "CLWMR:550 mb"',      & !        ...                            102
            ' "CLWMR:600 mb"',      & !        ...                            103
            ' "CLWMR:650 mb"',      & !        ...                            104
            ' "CLWMR:700 mb"',      & !        ...                            105
            ' "CLWMR:750 mb"',      & !        ...                            106 
            ' "CLWMR:800 mb"',      & !        ...                            107
            ' "CLWMR:850 mb"',      & !        ...                            108
            ' "CLWMR:900 mb"',      & !        ...                            109
            ' "CLWMR:925 mb"',      & !        ...                            110
            ' "CLWMR:950 mb"',      & !        ...                            111
            ' "CLWMR:975 mb"',      & !        ...                            112
            ' "CLWMR:1000 mb"',     & !        ...                            113 !+++++
            !+++++++++++
            ' "O3MR:150 mb"',       & !        ...                            114
            ' "O3MR:200 mb"',       & !        ...                            115
            ' "O3MR:250 mb"',       & !        ...                            116
            ' "O3MR:300 mb"',       & !        ...                            117
            ' "O3MR:350 mb"',       & !        ...                            118
            ' "O3MR:400 mb"',       & !        ...                            119
            ' "RH:10 mb"',          & !        ...                            120
            ' "RH:20 mb"',          & !        ...                            121
            ' "RH:30 mb"',          & !        ...                            122
            ' "RH:50 mb"',          & !        ...                            123
            ' "RH:70 mb"'/            !        ...                            124
                      
                                                                                
!*******
! 2. begin program
  
!===== 2.1. initialize [set names]
  L1 = LEN(trim(nwp_data_path))
  L2 = LEN(trim(grib_name)) 
  !print*,id,grib_name
  if (fylat_nwp_opt == 1 .or. fylat_nwp_opt == 4) then
      bin_name = 'fnl_'//nwptime%year(id)//nwptime%month(id)//nwptime%day(id)//'_'//nwptime%hour(id)//'_00'
  endif 

  if (fylat_nwp_opt == 2) then
      bin_name = 'gfs1p00_'//nwptime%year(id)//nwptime%month(id)//nwptime%day(id)//'_'//nwptime%hour(id)//'_00'
  endif 
    
  if (fylat_nwp_opt == 3) then
      bin_name = 'T639_'//nwptime%year(id)//nwptime%month(id)//nwptime%day(id)//'_'//nwptime%hour(id)//'_00'
  endif 

  if (fylat_nwp_opt == 5) then
      bin_name = 'gfs0p50_'//nwptime%year(id)//nwptime%month(id)//nwptime%day(id)//'_'//nwptime%hour(id)//'_00'
  endif 

  if (fylat_nwp_opt == 6) then
      bin_name = 'GRAPES_GFS_'//nwptime%year(id)//nwptime%month(id)//nwptime%day(id)//'_'//nwptime%hour(id)//'_00'
  endif 
  
  if (fylat_nwp_opt == 7) then
      bin_name = 'gdas1_0p25_'//nwptime%year(id)//nwptime%month(id)//nwptime%day(id)//'_'//nwptime%hour(id)//'_00'
  endif 

  IF (fylat_nwp_opt == 8) THEN   !9=gfs 0p25   added by minmin 20190516 
     !CALL find_gfs0p25_name (year, month, day, hour, filename, stime) 
     !filename_gfs0p25  = 'gfs0p25_'//year//month//day//'_'//hour//'_00'   !gfs 1p00
     bin_name  = 'gfs0p25_'//nwptime%year(id)//nwptime%month(id)//nwptime%day(id)//'_'//nwptime%hour(id)//'_00'
  ENDIF

  IF (fylat_nwp_opt == 9) THEN    !9=gfs 0.5*0.5 @41-layers (grib2)
     bin_name  = 'gfs0p50_41L_'//nwptime%year(id)//nwptime%month(id)//nwptime%day(id)//'_'//nwptime%hour(id)//'_00'
  ENDIF

  IF (fylat_nwp_opt == 10) THEN   !9=gfs 0.25*0.25 @41-layers (grib2)
     bin_name  = 'gfs0p25_41L_'//nwptime%year(id)//nwptime%month(id)//nwptime%day(id)//'_'//nwptime%hour(id)//'_00'
  ENDIF
  
  bin_name = TRIM(nwp_data_path)//trim(bin_name)  
  ! set grib and binary name [full name with path directory]
  
  ! check grib file 
  !status = checkfile(TRIM(grib_name))
  !if (status /= 0) then
  !   call file_message('Numerical Weather Prediction grib data', status)
  !endif  

!===== 2.2. convert grib DATA to binary DATA
  status = checkfile(bin_name)

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


if (status /= 0) then

   ! set path of wgrib.exe 
   wgrib_exe  = TRIM(code_root_path)//'/wgrib/wgrib' 
   wgrib2_exe = TRIM(TRIM(code_root_path)//'/wgrib/NCEP_grib2_to_bin_uv.sh')
   wgrib3_exe = TRIM(TRIM(code_root_path)//'/wgrib/NCEP_grib2_to_bin_uv_0p25.sh')
   wgrib4_exe = TRIM(TRIM(code_root_path)//'/wgrib/NCEP_grib2_to_bin_uv_41Layers.sh')
   wgrib5_exe = TRIM(TRIM(code_root_path)//'/wgrib/NCEP_grib2_to_bin_uv_0p25_41Layers.sh')
   wgrib2_grapes_gfs_exe = TRIM(TRIM(code_root_path)//'/wgrib/GRAPES_GFS_grib2_to_bin_uv.sh')
   T639_exe   = TRIM(TRIM(code_root_path)//'/wgrib/T639_to_bin.sh')  ! script for converting grib2 to binary for T639 data

   ! convert nwp grib to binary format
   print*,'*************************************'
   print*,'* convert nwp grib to binary format *'
   print*,'* in processing ......              *' 
   print*,'* ORG nwp data: ',trim(grib_name)
   print*,'* Converted data: ',trim(bin_name)
   if (fylat_nwp_opt == 1) then   ! 1=ncep and 2=gfs
      do i=1,113+52
         call system(wgrib_exe//'-s '//TRIM(grib_name)//' | grep -E '''//TRIM(list(i))//''' | '&
                    //wgrib_exe//'-i -nh '//TRIM(grib_name)//' -nh -append -o '//TRIM(bin_name)//'>/dev/null') 
         !call system(wgrib_exe//' -s '//TRIM(grib_name)//' | find "ALBdo" | '&
         !            //wgrib_exe//' -i -nh '//TRIM(grib_name)//' -nh -append -o '//TRIM(bin_name)) 
      enddo    
   endif
   
   if (fylat_nwp_opt == 3) then   ! T639
      call system(T639_exe//' '//trim(grib_name)//' '//TRIM(nwp_data_path)//'/')
   endif

   if (fylat_nwp_opt == 2 .or. fylat_nwp_opt == 4 .or. fylat_nwp_opt == 5) then   ! ncep/gfs grib2
      !do k=1,113
         !call system('wgrib2 -set local_table 1 -s '//TRIM(grib_name)//' | grep '''//TRIM(list(i))//''' |'&
         !             //'wgrib2 -i '//TRIM(grib_name)//' -no_header -append -bin '//TRIM(bin_name)//'>/dev/null')
         !call system('wgrib2 -s '//TRIM(grib_name)//' | grep'//TRIM(list2(k))//' | wgrib2 -i '//TRIM(grib_name)//' -no_header -append -bin '//TRIM(bin_name)//'>/dev/null')
      !enddo
      call system(wgrib2_exe//' '//TRIM(grib_name)//' '//TRIM(bin_name))
   endif

   IF (fylat_nwp_opt == 6) THEN   ! grapes gfs 0p25
   !print*,TRIM(wgrib2_grapes_gfs_exe)//' '//TRIM(grib_name)//' '//TRIM(bin_name)
      CALL system(TRIM(wgrib2_grapes_gfs_exe)//' '//TRIM(grib_name)//' '//TRIM(bin_name))
   ENDIF

   IF (fylat_nwp_opt == 7 .or. fylat_nwp_opt == 8) THEN   ! gfs 0.p25
      CALL system(wgrib3_exe//' '//TRIM(grib_name)//' '//TRIM(bin_name))
   ENDIF

   IF (fylat_nwp_opt == 9) THEN   ! gfs 0.5  @41 layers
      CALL system(wgrib4_exe//' '//TRIM(grib_name)//' '//TRIM(bin_name))
   ENDIF

   IF (fylat_nwp_opt == 10) THEN   ! gfs 0.p25 @41 layers
      CALL system(wgrib5_exe//' '//TRIM(grib_name)//' '//TRIM(bin_name))
   ENDIF

   print*,'* convert over                      *'
   print*,'*************************************'   
  
else

   print*,'*************************************'
   print*,'* The binary nwp file has existed   *'
   print*,'* ORG nwp data: ',trim(grib_name)
   print*,'* Converted data: ',trim(bin_name),' is OK!'
   print*,'*************************************' 
  
endif
              
! 3. end subroutine   
end subroutine convert_grib_to_binary
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~ subroutine 3: read nwp DATA with binary format ~~~~
subroutine read_nwp_arrays(bin_name, num)

!-----------------------------------------------------------------------
! !F90 read_nwp_arrays
!
! !description:
!    This program is to read nwp arrays in binary format and make the 
!    necessary nwp%DATA for next calculations.
!
! !Input  parameters:
!    bin_name        = nwp DATA's name
!    num             = number of nwp DATA [1 or 2]
!
! !Output parameters:
!    none
!
!-----------------------------------------------------------------------

! use modules

! 1. define variables
!===== 1.1.input information binary DATA's name
character(*), intent(in)          :: bin_name
integer(kind=4), intent(in)       :: num

!===== 1.2. other
integer(kind=4) :: status, io_err, SK, CK, i, j, k

!===== 1.3. variables 
!real(kind=4), dimension(1:360,1:181,1:124) :: bin_data
real(kind=4), allocatable, dimension(:,:,:) :: bin_data
real(kind=4), dimension(1:26) :: pre

data pre /10.0,  20.0,   30.0,  50.0,  70.0, 100.0, 150.0, 200.0, 250.0, &
          300.0, 350.0, 400.0, 450.0, 500.0, 550.0, 600.0, 650.0, 700.0, &
          750.0, 800.0, 850.0, 900.0, 925.0, 950.0, 975.0, 1000.0/

!*******
! 2. begin program

!===== 2.1. initialize
  ! check and open file 
  status = checkfile(TRIM(bin_name))
  if (status /= 0) then
     call file_message('Numerical Weather Prediction binary data', status)
  endif 
          
  ! open the ncep binary DATA     
  if (fylat_nwp_opt == 1 ) then  ! grib1
     open(11,file=TRIM(bin_name),status='old',access='direct',form="binary",&
          recl=181*360*(113+52)*4,IOSTAT=io_err)
  endif
  
  if (fylat_nwp_opt == 2 .or. fylat_nwp_opt == 4) then  ! grib2
     open(11,file=TRIM(bin_name),status='old',access='direct',form="binary",&
          recl=181*360*(124+52)*4,IOSTAT=io_err)     
  endif

  if (fylat_nwp_opt == 5) then  ! grib2
     open(11,file=TRIM(bin_name),status='old',access='direct',form="binary",&
          recl=361*720*(124+52)*4,IOSTAT=io_err)     
  endif
       
  ! check the opened binary DATA 
  if (io_err /= 0) then    
     call bin_info_message('NCEP binary DATA')
  endif
  
!===== 2.2. read binary DATA
  if (fylat_nwp_opt /= 5) then  ! grib2
  
     allocate(bin_data(360,181,124+52))
  
     if (fylat_nwp_opt == 1) then  ! grib1
        read(11,rec=1) (((bin_data(i,j,k),i=1,360),j=1,181),k=1,113+52)
        close(11) 
     endif
  
     if (fylat_nwp_opt == 2 .or. fylat_nwp_opt == 4) then  ! grib2
        read(11,rec=1) (((bin_data(i,j,k),i=1,360),j=1,181),k=1,124+52)
        close(11) 
     endif
  
  endif

  if (fylat_nwp_opt == 5) then  ! grib2
  
     allocate(bin_data(720,361,124+52))
     read(11,rec=1) (((bin_data(i,j,k),i=1,720),j=1,361),k=1,124+52)
     close(11) 
     
  endif
     
!===== 2.3. delete binary DATA
!  call system('rm -rf '//TRIM(bin_name))

if (fylat_nwp_opt /= 5) then  ! grib2  1p00
!===== 2.4. assign variables
!xxx lon lat
do i=1,360
do j=1,181
   nwpo%lon(i,j) = -179.5+(i-1)
   nwpo%lat(i,j) = 90.-(j-1)
enddo
enddo

!xxx (1) plev_nointerp
nwpo%plev_nointerp(1:26,num) = pre

do j=1,181
!xxx (2) psfc  ! the matrix transform is due to agree with longtitude distributions
  if (fylat_nwp_opt == 1) then  ! grib1
      SK=j
      !nwp%first_lat2 = 90.0
  endif
  if (fylat_nwp_opt == 2 .or. fylat_nwp_opt == 4) then  ! grib2
      SK=182-j
      !nwp%first_lat2  = -90.0
  endif

nwpo%psfc(1:180,j,num)   = bin_data(181:360,SK,1)
nwpo%psfc(181:360,j,num) = bin_data(1:180,SK,1)

!xxx (3) pmsl 
nwpo%pmsl(1:180,j,num)   = bin_data(181:360,SK,2)
nwpo%pmsl(181:360,j,num) = bin_data(1:180,SK,2)

!xxx (4) tsfc
nwpo%tsfc(1:180,j,num)   = bin_data(181:360,SK,3)
nwpo%tsfc(181:360,j,num) = bin_data(1:180,SK,3)

!xxx (5) zsfc
nwpo%zsfc(1:180,j,num)   = bin_data(181:360,SK,4)
nwpo%zsfc(181:360,j,num) = bin_data(1:180,SK,4)

!xxx (6) albedo
nwpo%albedo(1:180,j,num)   = bin_data(181:360,SK,5)
nwpo%albedo(181:360,j,num) = bin_data(1:180,SK,5)

!xxx (7) t_sigma
nwpo%t_sigma(1:180,j,num)   = bin_data(181:360,SK,6)
nwpo%t_sigma(181:360,j,num) = bin_data(1:180,SK,6)

!xxx (8) rh_sigma
nwpo%rh_sigma(1:180,j,num)  = bin_data(181:360,SK,7)
nwpo%rh_sigma(181:360,j,num)= bin_data(1:180,SK,7)

!xxx (9) u_sigma
nwpo%u_sigma(1:180,j,num)   = bin_data(181:360,SK,8)
nwpo%u_sigma(181:360,j,num) = bin_data(1:180,SK,8)

!xxx (10) v_sigma
nwpo%v_sigma(1:180,j,num)   = bin_data(181:360,SK,9)
nwpo%v_sigma(181:360,j,num) = bin_data(1:180,SK,9)

!xxx (11) tpw
nwpo%tpw(1:180,j,num)   = bin_data(181:360,SK,10)
nwpo%tpw(181:360,j,num) = bin_data(1:180,SK,10)

!xxx (12) weasd
nwpo%weasd(1:180,j,num)   = bin_data(181:360,SK,11)
nwpo%weasd(181:360,j,num) = bin_data(1:180,SK,11)

!xxx (13) o3col
nwpo%o3col(1:180,j,num)   = bin_data(181:360,SK,12)
nwpo%o3col(181:360,j,num) = bin_data(1:180,SK,12)

!xxx (14) ttropo
nwpo%ttropo(1:180,j,num)    = bin_data(181:360,SK,13)
nwpo%ttropo(181:360,j,num)  = bin_data(1:180,SK,13)

!xxx (15) tlev
do i=1,26
   nwpo%tlev(1:180,j,i,num)    = bin_data(181:360,SK,13+i)
   nwpo%tlev(181:360,j,i,num)  = bin_data(1:180,SK,13+i)  !last = 39
enddo

!xxx (16) zlev
do i=1,26
   nwpo%zlev(1:180,j,i,num)    = bin_data(181:360,SK,39+i)
   nwpo%zlev(181:360,j,i,num)  = bin_data(1:180,SK,39+i)  !last = 65
enddo

!xxx (17) o3lev
nwpo%o3lev(:,j,:,num)   = -999.0
do i=1,6
   nwpo%o3lev(1:180,j,i,num)   = bin_data(181:360,SK,65+i)
   nwpo%o3lev(181:360,j,i,num) = bin_data(1:180,SK,65+i)  !last = 71
enddo

!xxx (18) rhlev
nwpo%rhlev(:,j,:,num)   = -999.0
do i=1,21
   nwpo%rhlev(1:180,j,i+5,num)   = bin_data(181:360,SK,71+i)
   nwpo%rhlev(181:360,j,i+5,num) = bin_data(1:180,SK,71+i)  !last = 92
enddo

!xxx (19) clwlev
nwpo%clwlev(:,j,:,num)   = -999.0
do i=1,21
   nwpo%clwlev(1:180,j,i+5,num)  = bin_data(181:360,SK,92+i)
   nwpo%clwlev(181:360,j,i+5,num)= bin_data(1:180,SK,92+i)  !last = 113
enddo

!xxx (20) add O3MR + RH
if (fylat_nwp_opt == 2 .or. fylat_nwp_opt == 4) then  ! grib2

! O3MR level 7,8,9,10,11,12
do i=1,6
   nwpo%o3lev(1:180,j,i+6,num)   = bin_data(181:360,SK,113+i)
   nwpo%o3lev(181:360,j,i+6,num) = bin_data(1:180,SK,113+i)  !last = 119
enddo
! RH  level 1,2,3,4,5
do i=1,5
   nwpo%rhlev(1:180,j,i,num)   = bin_data(181:360,SK,119+i)
   nwpo%rhlev(181:360,j,i,num) = bin_data(1:180,SK,119+i)  !last = 124
enddo

!xxx ulev
do i=1,26
   nwpo%ulev(1:180,j,i,num)    = bin_data(181:360,SK,124+i)
   nwpo%ulev(181:360,j,i,num)  = bin_data(1:180,SK,124+i)  !
enddo

!xxx vlev
do i=1,26
   nwpo%vlev(1:180,j,i,num)    = bin_data(181:360,SK,124+26+i)
   nwpo%vlev(181:360,j,i,num)  = bin_data(1:180,SK,124+26+i)  !
enddo


endif

enddo

endif !if (fylat_nwp_opt /= 5) then  ! grib2  1p00


if (fylat_nwp_opt == 5) then  ! grib2 0p50

!xxx lon lat
do i=1,720
do j=1,361
   nwpo%lon(i,j) = -179.75+(i-1)*0.5
   nwpo%lat(i,j) = 90.-(j-1)*0.5
enddo
enddo

!xxx (1) plev_nointerp
nwpo%plev_nointerp(1:26,num) = pre

do j=1,361
!xxx (2) psfc  ! the matrix transform is due to agree with longtitude distributions

! gfs0p50 grib2
SK=362-j

nwpo%psfc(1:360,j,num)   = bin_data(361:720,SK,1)
nwpo%psfc(361:720,j,num) = bin_data(1:360,SK,1)

!xxx (3) pmsl 
nwpo%pmsl(1:360,j,num)   = bin_data(361:720,SK,2)
nwpo%pmsl(361:720,j,num) = bin_data(1:360,SK,2)

!xxx (4) tsfc
nwpo%tsfc(1:360,j,num)   = bin_data(361:720,SK,3)
nwpo%tsfc(361:720,j,num) = bin_data(1:360,SK,3)

!xxx (5) zsfc
nwpo%zsfc(1:360,j,num)   = bin_data(361:720,SK,4)
nwpo%zsfc(361:720,j,num) = bin_data(1:360,SK,4)

!xxx (6) albedo
nwpo%albedo(1:360,j,num)   = bin_data(361:720,SK,5)
nwpo%albedo(361:720,j,num) = bin_data(1:360,SK,5)

!xxx (7) t_sigma
nwpo%t_sigma(1:360,j,num)   = bin_data(361:720,SK,6)
nwpo%t_sigma(361:720,j,num) = bin_data(1:360,SK,6)

!xxx (8) rh_sigma
nwpo%rh_sigma(1:360,j,num)  = bin_data(361:720,SK,7)
nwpo%rh_sigma(361:720,j,num)= bin_data(1:360,SK,7)

!xxx (9) u_sigma
nwpo%u_sigma(1:360,j,num)   = bin_data(361:720,SK,8)
nwpo%u_sigma(361:720,j,num) = bin_data(1:360,SK,8)

!xxx (10) v_sigma
nwpo%v_sigma(1:360,j,num)   = bin_data(361:720,SK,9)
nwpo%v_sigma(361:720,j,num) = bin_data(1:360,SK,9)

!xxx (11) tpw
nwpo%tpw(1:360,j,num)   = bin_data(361:720,SK,10)
nwpo%tpw(361:720,j,num) = bin_data(1:360,SK,10)

!xxx (12) weasd
nwpo%weasd(1:360,j,num)   = bin_data(361:720,SK,11)
nwpo%weasd(361:720,j,num) = bin_data(1:360,SK,11)

!xxx (13) o3col
nwpo%o3col(1:360,j,num)   = bin_data(361:720,SK,12)
nwpo%o3col(361:720,j,num) = bin_data(1:360,SK,12)

!xxx (14) ttropo
nwpo%ttropo(1:360,j,num)    = bin_data(361:720,SK,13)
nwpo%ttropo(361:720,j,num)  = bin_data(1:360,SK,13)

!xxx (15) tlev
do i=1,26
   nwpo%tlev(1:360,j,i,num)    = bin_data(361:720,SK,13+i)
   nwpo%tlev(361:720,j,i,num)  = bin_data(1:360,SK,13+i)  !last = 39
enddo

!xxx (16) zlev
do i=1,26
   nwpo%zlev(1:360,j,i,num)    = bin_data(361:720,SK,39+i)
   nwpo%zlev(361:720,j,i,num)  = bin_data(1:360,SK,39+i)  !last = 65
enddo

!xxx (17) o3lev
nwpo%o3lev(:,j,:,num)   = -999.0
do i=1,6
   nwpo%o3lev(1:360,j,i,num)   = bin_data(361:720,SK,65+i)
   nwpo%o3lev(361:720,j,i,num) = bin_data(1:360,SK,65+i)  !last = 71
enddo

!xxx (18) rhlev
nwpo%rhlev(:,j,:,num)   = -999.0
do i=1,21
   nwpo%rhlev(1:360,j,i+5,num)   = bin_data(361:720,SK,71+i)
   nwpo%rhlev(361:720,j,i+5,num) = bin_data(1:360,SK,71+i)  !last = 92
enddo

!xxx (19) clwlev
nwpo%clwlev(:,j,:,num)   = -999.0
do i=1,21
   nwpo%clwlev(1:360,j,i+5,num)  = bin_data(361:720,SK,92+i)
   nwpo%clwlev(361:720,j,i+5,num)= bin_data(1:360,SK,92+i)  !last = 113
enddo

!xxx (20) add O3MR + RH

! O3MR level 7,8,9,10,11,12
do i=1,6
   nwpo%o3lev(1:360,j,i+6,num)   = bin_data(361:720,SK,113+i)
   nwpo%o3lev(361:720,j,i+6,num) = bin_data(1:360,SK,113+i)  !last = 119
enddo
! RH  level 1,2,3,4,5
do i=1,5
   nwpo%rhlev(1:360,j,i,num)   = bin_data(361:720,SK,119+i)
   nwpo%rhlev(361:720,j,i,num) = bin_data(1:360,SK,119+i)  !last = 124
enddo

!xxx ulev
do i=1,26
   nwpo%ulev(1:360,j,i,num)    = bin_data(361:720,SK,124+i)
   nwpo%ulev(361:720,j,i,num)  = bin_data(1:360,SK,124+i)  !
enddo

!xxx vlev
do i=1,26
   nwpo%vlev(1:360,j,i,num)    = bin_data(361:720,SK,124+26+i)
   nwpo%vlev(361:720,j,i,num)  = bin_data(1:360,SK,124+26+i)  !
enddo

enddo

endif !if (fylat_nwp_opt == 5) then  ! grib2 0p50

deallocate(bin_data)

! 3. end subroutine   
end subroutine read_nwp_arrays
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~ subroutine 4: allocate_nwpo_arrays ~~~~~~~~~~~~~~~~
subroutine allocate_nwpo_arrays(nwp_id)
 
!-----------------------------------------------------------------------
! !F90 allocate_nwpo_arrays
!
! !description:
!    This program is to allocate original nwp data.
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
integer(kind=1) :: nwp_id
integer(kind=4) :: x, y, z, c
  
! 2. begin program 
  if (nwp_id == 1 .or. nwp_id == 2 .or. nwp_id == 4) then   ! 1/4=ncep grib1/2 and 2/5=gfs grib1/2
      x = 360
      y = 181
      z = 26
      c = 2
  endif
  
  if (nwp_id == 3) then     ! 3=T639
      x = 1280
      y = 641
      z = 36
      c = 2
  endif

  if (nwp_id == 5 ) then   ! 1/4=ncep grib1/2 and 2/5=gfs grib1/2
      x = 720
      y = 361
      z = 26
      c = 2
  endif

  IF (nwp_id == 6) THEN     ! 6=grapes gfs 0.25*0.25 grib2
      x = 1440
      y = 720
      z = 40
      c = 2
  ENDIF

  IF (nwp_id == 7 .or. nwp_id == 8) THEN     ! 7=gdas1 0p25 grib2  ;8=gfs0p25 grib2    
      x = 1440 !721
      y = 721 !360
      z = 31
      c = 2
  ENDIF

  IF (nwp_id == 9) THEN     ! 9=gfs0p5 grib2  @ 41 layers
      x = 720
      y = 361
      z = 41
      c = 2
  ENDIF

  IF (nwp_id == 10) THEN     ! 10=gfs0p25 grib2  @ 41 layers
      x = 1440 !721
      y = 721 !360
      z = 41
      c = 2
  ENDIF
  
  allocate(nwpo%lon(x,y), nwpo%lat(x,y),              & 
           nwpo%plev_nointerp(z,c),                   &
           nwpo%psfc(x,y,c),nwpo%pmsl(x,y,c),         &
           nwpo%tsfc(x,y,c),nwpo%zsfc(x,y,c),         &
           nwpo%albedo(x,y,c),nwpo%t_sigma(x,y,c),    &
           nwpo%rh_sigma(x,y,c),nwpo%u_sigma(x,y,c),  &
           nwpo%v_sigma(x,y,c),nwpo%tpw(x,y,c),       &
           nwpo%weasd(x,y,c),nwpo%o3col(x,y,c),       &
           nwpo%ttropo(x,y,c),                        &
           nwpo%tlev(x,y,z,c),                        &
           nwpo%zlev(x,y,z,c),                        &  
           nwpo%o3lev(x,y,z,c),                       &    
           nwpo%rhlev(x,y,z,c),                       &  
           nwpo%clwlev(x,y,z,c),                      &
           nwpo%ulev(x,y,z,c),                        &  
           nwpo%vlev(x,y,z,c),                        &
           stat=astatus)
           
  if (astatus /= 0) then
     print *,"(a,'Not enough memory to allocate L1b data structure.')"
     stop
  endif
 
! 3. end subroutine
end subroutine allocate_nwpo_arrays
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~ subroutine 5: deallocate nwpo arrays ~~~~~~~~~~~~~~ 
subroutine deallocate_nwpo_arrays
 
!-----------------------------------------------------------------------
! !F90 deallocate_nwpo_arrays
!
! !description:
!    This program is to deallocate original nwp data.
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
integer(kind=4)  :: astatus
 
! 2. begin program 
  deallocate(nwpo%lon, nwpo%lat,          &
             nwpo%plev_nointerp,          &
             nwpo%psfc,nwpo%pmsl,         &
             nwpo%tsfc,nwpo%zsfc,         &
             nwpo%albedo,nwpo%t_sigma,    &
             nwpo%rh_sigma,nwpo%u_sigma,  &
             nwpo%v_sigma,nwpo%tpw,       &
             nwpo%weasd,nwpo%o3col,       &
             nwpo%ttropo,                 &
             nwpo%tlev,                   &
             nwpo%zlev,                   &
             nwpo%o3lev,                  &
             nwpo%rhlev,                  &
             nwpo%clwlev,                 &
             nwpo%ulev,                   &
             nwpo%vlev,                   &
             stat=astatus)
           
  if (astatus /= 0) then
     print *,"(a,'Error deallocating nwpo DATA structure.')"
     stop
  endif
 
! 3. end subroutine
end subroutine deallocate_nwpo_arrays
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~ subroutine 6: FIND_nwp_EQAREA_CELL ~~~~~~~~~~~~~~~~ 
subroutine find_nwp_eqarea_cell(lon, lat, ilon, ilat, first_lat)
 
!-----------------------------------------------------------------------
! !F90 FIND_nwp_EQAREA_CELL
!
! !description:
!    This program is to convert lat, lon into GFS grid-cell.
!
! !Input  parameters:
!    none
!
! !Output parameters:
!    none
!
!-----------------------------------------------------------------------

!use names_module ,ONLY: fylat_nwp_opt

real(kind=4), dimension(:,:), target, intent(in) :: lon, lat
real(kind=4)                        , intent(in) :: first_lat
!integer(kind=1), dimension(:,:), intent(in)      :: space_mask
integer(kind=4), dimension(:,:), intent(out)     :: ilon, ilat

real(kind=4) :: rlon, rlat
integer(kind=4) :: factor, nx, ny, i, j, k

! 2. begin program 
    
    nx = size(lon,dim=1)
    ny = size(lat,dim=2)
    
    ilon = missing_value_int4
    ilat = missing_value_int4
    
    factor = 1
    if (nwp%first_lat > 0.0) factor = -1
    !if (first_lat > 0.0) factor = -1
    
    do j = 1, ny
    do i = 1, nx
    !print*,j,i,lon(i,j),lat(i,j),space_mask(i,j)
      !if (space_mask(i,j) == sym%NO_SPACE) then    
          rlon = lon(i,j)
          rlat = lat(i,j)
          !if (rlon < 0.0) rlon = rlon + 360.0
          if (fylat_nwp_opt == 1 .or. fylat_nwp_opt == 2 .or. fylat_nwp_opt == 4) then   ! 1=ncep and 2=gfs
             ilat(i,j) = max(1, min(nwp%nlat, int ((rlat - nwp%first_lat + 0.5*nwp%dlat*factor) / (nwp%dlat*factor) + 1)))
             !ilat(i,j) = max(1, min(nwp%nlat, int ((rlat - first_lat + 0.5*nwp%dlat*factor) / (nwp%dlat*factor) + 1)))
             ilon(i,j) = max(1, min(nwp%nlon, int (((rlon - nwp%first_lon + 0.5*nwp%dlon) / nwp%dlon) + 1)))  
             !print*,ilat(i,j),ilon(i,j),rlat,rlon
          endif

          if (fylat_nwp_opt == 5 .or. fylat_nwp_opt == 9) then   ! 1=ncep and 2=gfs
             ilat(i,j) = max(1, min(nwp%nlat05, int ((rlat - nwp%first_lat + 0.5*nwp%dlat05*factor) / (nwp%dlat05*factor) + 1)))
             ilon(i,j) = max(1, min(nwp%nlon05, int (((rlon - nwp%first_lon + 0.5*nwp%dlon05) / nwp%dlon05) + 1)))  
             !print*,ilat(i,j),ilon(i,j),rlat,rlon
          endif
                    
          if (fylat_nwp_opt == 3) then   ! 3=T639  resolution = 0.2815/0.2823
             ilat(i,j) = max(1, min(nwp%nlat_T639, int ((rlat - nwp%first_lat + 0.5*nwp%dlat_T639*factor) / (nwp%dlat_T639*factor) + 1)))
             !ilat(i,j) = max(1, min(nwp%nlat_T639, int ((rlat - first_lat + 0.5*nwp%dlat_T639*factor) / (nwp%dlat_T639*factor) + 1)))
             ilon(i,j) = max(1, min(nwp%nlon_T639, int (((rlon - nwp%first_lon + 0.5*nwp%dlon_T639) / nwp%dlon_T639) + 1))) 
             
             ! small resolution = 0.563 1.0 degree
             !ilat(i,j) = max(2, 2*min(321, int ((rlat - nwp%first_lat + 0.5*0.563*factor) / (0.563*factor) + 1)))
             !ilon(i,j) = max(2, 2*min(640, int (((rlon - nwp%first_lon + 0.5*0.563) / 0.563) + 1))) 
             !print*,ilat(i,j),ilon(i,j),rlat,rlon
             
             ! small resolution = 1.0 degree
             !ilat(i,j) = max(1, min(nwp%nlat, int ((rlat - nwp%first_lat + 0.5*nwp%dlat*factor) / (nwp%dlat*factor) + 1)))
             !ilon(i,j) = max(1, min(nwp%nlon, int (((rlon - nwp%first_lon + 0.5*nwp%dlon) / nwp%dlon) + 1)))  
          endif 

          IF (fylat_nwp_opt == 6) THEN   ! 6=grapes gfs resolution = 0.25
             ilat(i,j) = max(1, min(nwp%nlat25, int ((rlat - nwp%first_lat + 0.5*0.25*factor) / (0.25*factor) + 1)))
             ilon(i,j) = max(1, min(nwp%nlon25, int (((rlon - nwp%first_lon + 0.5*0.25) / 0.25) + 1))) 
          ENDIF            

          IF (fylat_nwp_opt == 7 .or. fylat_nwp_opt == 8 .or. fylat_nwp_opt == 10) THEN   ! 8=gfs 0.25
             ilat(i,j) = max(1, min(nwp%nlat25, int ((rlat - nwp%first_lat + 0.5*nwp%dlat0p25*factor) / (nwp%dlat0p25*factor) + 1)))
             ilon(i,j) = max(1, min(nwp%nlon25, int (((rlon - nwp%first_lon + 0.5*nwp%dlon0p25) / nwp%dlon0p25) + 1)))  
             !PRINT*,ilat(i,j),ilon(i,j),rlat,rlon
          ENDIF     

      ! endif
        
    end do
    end do

! 3. end subroutine
end subroutine find_nwp_eqarea_cell
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~ subroutine 3: read nwp DATA with binary format ~~~~
subroutine read_nwp_T639_arrays(bin_name, num)

!-----------------------------------------------------------------------
! !F90 read_nwp_T639_arrays
!
! !description:
!    This program is to read nwp arrays in binary format and make the 
!    necessary nwp%DATA for next calculations.
!
! !Input  parameters:
!    bin_name        = nwp T639 data's name
!    num             = number of nwp DATA [1 or 2]
!
! !Output parameters:
!    none
!
!-----------------------------------------------------------------------

! use modules

! 1. define variables
!===== 1.1.input information binary DATA's name
character(*), intent(in)             :: bin_name
integer(kind=4), intent(in)       :: num

!===== 1.2. other
integer(kind=4) :: status, io_err, SK, i, j, k

!===== 1.3. variables 
!real(kind=4), dimension(1:1280,1:641,1:118) :: bin_data
real(kind=4), allocatable, dimension(:,:,:) :: bin_data
real(kind=4), dimension(1:36) :: pre
real(kind=4)                  :: reso

!*******
! 2. begin program

!===== 2.1. initialize

pre(1:36) = (/  0.1,   0.2,   0.5,   1.0,   1.5,   2.0,   3.0,   4.0,    &
                5.0,   7.0,   10.,   20.,   30.,   50.,   70.,  100.,    &
               150.,  200.,  250.,  300.,  350.,  400.,  450.,  500.,    &
               550.,  600.,  650.,  700.,  750.,  800.,  850.,  900.,    &
               925.,  950.,  975., 1000./)
               
  ! check and open file 
  status = checkfile(TRIM(bin_name))
  if (status /= 0) then
     call file_message('Numerical Weather Prediction binary data', status)
  endif 
          
  ! open the T639 binary data   
  open(11,file=TRIM(bin_name),status='old',access='direct',form="binary",&
         recl=1280*641*118*4,IOSTAT=io_err)

  ! check the opened binary DATA 
  if (io_err /=0) then    
     call bin_info_message('T639 binary data')
  endif
  
!===== 2.2. read binary DATA
  allocate(bin_data(1280,641,118))
  read(11,rec=1) (((bin_data(i,j,k),i=1,1280),j=1,641),k=1,118)
  close(11)
  
!===== 2.3. delete binary DATA
 ! call system('rm -rf '//TRIM(bin_name))
  
!===== 2.4. assign variables
!xxx lon lat
reso = 0.2815
do i=1,1280
do j=1,641
   nwpo%lon(i,j) = -179.85925+(i-1)*reso
   nwpo%lat(i,j) = 90.-(j-1)*reso
   !nwpo%lat(i,j) = -90.+(j-1)*reso
enddo
enddo

!xxx (1) plev_nointerp
nwpo%plev_nointerp(1:36,num) = pre

!xxx (2) psfc  ! the matrix transform is due to agree with longtitude distributions
do j=1,641

SK=642-j
!nwp%first_lat2  = -90.0
!SK=j

nwpo%psfc(1:640,j,num)   = bin_data(641:1280,SK,1)
nwpo%psfc(641:1280,j,num) = bin_data(1:640,SK,1)

!xxx (3) pmsl 
nwpo%pmsl(1:640,j,num)   = bin_data(641:1280,SK,2)
nwpo%pmsl(641:1280,j,num) = bin_data(1:640,SK,2)

!xxx (4) tsfc
nwpo%tsfc(1:640,j,num)   = bin_data(641:1280,SK,3)
nwpo%tsfc(641:1280,j,num) = bin_data(1:640,SK,3)

!xxx (5) zsfc
nwpo%zsfc(1:640,j,num)   = bin_data(641:1280,SK,4)
nwpo%zsfc(641:1280,j,num) = bin_data(1:640,SK,4)

!xxx (6) albedo [sfc pre]
nwpo%albedo(1:640,j,num)   = bin_data(641:1280,SK,5)
nwpo%albedo(641:1280,j,num) = bin_data(1:640,SK,5)

!xxx (7) t_sigma
nwpo%t_sigma(1:640,j,num)   = bin_data(641:1280,SK,6)
nwpo%t_sigma(641:1280,j,num) = bin_data(1:640,SK,6)

!xxx (8) rh_sigma
nwpo%rh_sigma(1:640,j,num)  = bin_data(641:1280,SK,7)
nwpo%rh_sigma(641:1280,j,num)= bin_data(1:640,SK,7)

!xxx (9) u_sigma
nwpo%u_sigma(1:640,j,num)   = bin_data(641:1280,SK,8)
nwpo%u_sigma(641:1280,j,num) = bin_data(1:640,SK,8)

!xxx (10) v_sigma
nwpo%v_sigma(1:640,j,num)   = bin_data(641:1280,SK,9)
nwpo%v_sigma(641:1280,j,num) = bin_data(1:640,SK,9)

!xxx (11) tpw
nwpo%tpw(1:640,j,num)   = -999.0
nwpo%tpw(641:1280,j,num) = -999.0

!xxx (12) weasd
nwpo%weasd(1:640,j,num)   = -999.0
nwpo%weasd(641:1280,j,num) = -999.0

!xxx (13) o3col
nwpo%o3col(1:640,j,num)   = -999.0
nwpo%o3col(641:1280,j,num) = -999.0

!xxx (14) ttropo
nwpo%ttropo(1:640,j,num)    = bin_data(641:1280,SK,10)
nwpo%ttropo(641:1280,j,num)  = bin_data(1:640,SK,10)

!xxx (15) tlev
do i=1,36
   nwpo%tlev(1:640,j,i,num)    = bin_data(641:1280,SK,10+i)
   nwpo%tlev(641:1280,j,i,num)  = bin_data(1:640,SK,10+i)  !last = 46
enddo

!xxx (16) zlev
do i=1,36
   nwpo%zlev(1:640,j,i,num)    = bin_data(641:1280,SK,46+i)
   nwpo%zlev(641:1280,j,i,num)  = bin_data(1:640,SK,46+i)  !last = 82
enddo

!xxx (17) o3lev
nwpo%o3lev(:,j,:,num)   = -999.0
!do i=1,36
!   nwpo%o3lev(1:640,:,i,num)   =   -999.0 ! bin_data(181:360,:,65+i)
!   nwpo%o3lev(641:1280,:,i,num) =  -999.0 ! bin_data(1:180,:,65+i)  !last = 71
!enddo

!xxx (18) rhlev
nwpo%rhlev(:,j,:,num)   = -999.0
do i=1,36
   nwpo%rhlev(1:640,j,i,num)   = bin_data(641:1280,SK,82+i)
   nwpo%rhlev(641:1280,j,i,num) = bin_data(1:640,SK,82+i)  !last = 92
enddo

!xxx (19) clwlev
nwpo%clwlev(:,j,:,num)   = -999.0
!do i=1,36
!   nwpo%clwlev(1:640,:,i,num)  =   -999.0 ! bin_data(181:360,:,92+i)
!   nwpo%clwlev(641:1280,:,i,num)=  -999.0 ! bin_data(1:180,:,92+i)  !last = 113
!enddo

enddo

deallocate(bin_data)
! 3. end subroutine   
end subroutine read_nwp_T639_arrays
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~ SUBROUTINE 1: extract satellite observation time ~~
subroutine extract_nwptime (nwp_path,xname,id,year,month,day,hour)

!-----------------------------------------------------------------------
! !F90 extract_sattime
!
! !Description:
!    This program is to extract nwp data time.
!
! !Input  parameters:
!    L1b_path        = index number of nwp option
!    xname           = nwp data name
!
! !Output parameters:
!    year   sattime%year     = integer year
!    month  sattime%month    = integer month
!    day    sattime%day      = integer day
!    hour   sattime%hour     = integer hour 
!
!-----------------------------------------------------------------------

!USE module

implicit none

! 1. define variables
!===== 1.1.input and out variables
!integer(kind=int1), intent(in):: index   ! index
character*(*), intent(in)     :: nwp_path        ! path
character*(*), intent(in)     :: xname     ! satname
character(len=4), intent(out) :: year
character(len=2), intent(out) :: month, day, hour

!===== 1.2.middle variables  
integer(kind=4), dimension(1:4)  :: fp, ep      ! the first [fp] and END [ep] position of time string  
integer(kind=4), dimension(1:9)  :: tt
character(len=2), dimension(1:9) :: ttn
integer(kind=1)                  :: id

!===== 1.3.other variables
integer length, ierr, leap_flg, day1, L1, L2, fname, i, j
character(len=200) :: nwpname

!******* these parameters are set for Metesat-8 seviri *********
!data fp /21, 25, 29, 31/  ! year, day, hour, mint
!data ep /24, 27, 30, 32/
!***************************************************************
! 2. begin program
  PRINT*,'  ... fylat extract nwp data time'
  
!xxxxxxxxxxxx
!Note: here, the satellite is FY3D/MERSI_II
!      IF we USE dIFferent satellite, we should change character length.
!      Variable fp(x) and ep(x) of position of string should be changed. 
!xxxxxxxxxxxx
!  fylat_nwp_opt	          
  !   1 = ncep reanalysis 1*1 (grib1)
  !   2 = gfs1p00 1*1 (grib2)
  !   3 = T639 0.125*0.125 (grib2)
  !   4 = ncep reanalysis 1*1 (grib2) 
  !   5 = gfs0p50 0.5*0.5 (grib2) 
  !   6 = grapes 0.25*0.25 (grib2) 
  !   7 = gdas1 1*1 (grib1) from UW ftp.ssec.wisc.edu

  L1 = LEN(trim(nwp_path))
  L2 = LEN(trim(xname)) 
  nwpname = xname(L1+5:L2)   
  write(*,*)'           nwp data full name = '//trim(nwpname)
  
  if (fylat_nwp_opt == 1 .or. fylat_nwp_opt == 4) then
     fp(1:4) = (/5,  9, 11, 14/) ! year, month, day, hour, mint
     ep(1:4) = (/8, 10, 12, 15/)
     year  = nwpname(fp(1):ep(1))
     month = nwpname(fp(2):ep(2))
     day   = nwpname(fp(3):ep(3))
     hour  = nwpname(fp(4):ep(4))
  endif
  !print*,nwpname
  !print*,year,month,day,hour

  if (fylat_nwp_opt == 77) then
     fp(1:4) = (/15, 17, 19, 22/) ! year, month, day, hour, mint
     ep(1:4) = (/16, 18, 20, 23/)
     year  = '20'//nwpname(fp(1):ep(1))
     month = nwpname(fp(2):ep(2))
     day   = nwpname(fp(3):ep(3))
     hour  = nwpname(fp(4):ep(4))
  endif
  !print*,nwpname
  !print*,year,month,day,hour  

  if (fylat_nwp_opt == 2 .or. fylat_nwp_opt == 3 .or. fylat_nwp_opt == 5 .or. fylat_nwp_opt == 6 &
      .or. fylat_nwp_opt == 7 .or. fylat_nwp_opt == 8  &
      .or. fylat_nwp_opt == 9 .or. fylat_nwp_opt == 10) then
     fp(1:4) = (/5,  9, 11, 14/) ! year, month, day, hour, mint
     ep(1:4) = (/8, 10, 12, 15/)
     
     tt(1:9)  = (/0, 3, 6, 9, 12, 15, 18, 21, 24/)
     ttn(1:9) = (/'00', '03', '06', '09', '12', '15', '18', '21', '24'/)     
     
     !ICNVRT(WAY,NUM,STRING,LENGTH,IERR)
     fname = sat%year
     CALL ICNVRT(0,fname,year,length,ierr)
     
     fname = sat%month
     CALL ICNVRT(0,fname,month,length,ierr)
     if (sat%month < 10) month = '0'//month

     fname = sat%day
     CALL ICNVRT(0,fname,day,length,ierr)   
     if (sat%day < 10) day = '0'//day
       
     if (id == 1) then
        do i = 1, 8
           if (sat%hour >= tt(i) .and. sat%hour < tt(i+1)) then
               hour = ttn(i)
               exit
           endif
        enddo
     endif

     if (id == 2) then
        do i = 1, 8
           if (sat%hour >= tt(i) .and. sat%hour < tt(i+1)) then
               hour = ttn(i+1)
               exit
           endif
        enddo
     endif
          
  endif


! 3. END SUBROUTINE   
end subroutine extract_nwptime
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~ SUBROUTINE 3: READ NWP DATA with binary format ~~~~
SUBROUTINE read_NWP_grapes_gfs_arrays(bin_name, num)

!-----------------------------------------------------------------------
! !F90 read_NWP_grapes_gfs_arrays
!
! !Description:
!    This program is to READ NWP arrays in binary format and make the 
!    necessary nwp%DATA for next calculations.
!
! !Input  parameters:
!    bin_name        = NWP T639 data's name
!    num             = number of NWP DATA [1 or 2]
!
! !Output parameters:
!    none
!
!-----------------------------------------------------------------------

! USE modules

! 1. define variables
!===== 1.1.input information binary DATA's name
CHARACTER(*), INTENT(in)             :: bin_name
INTEGER(kind=4), INTENT(in)       :: num

!===== 1.2. other
INTEGER(kind=4) :: status, io_err, S, nvar

!===== 1.3. variables 
!REAL(kind=4), DIMENSION(1:1280,1:641,1:118) :: bin_data
REAL(kind=4), ALLOCATABLE, DIMENSION(:,:,:) :: bin_data
REAL(kind=4), DIMENSION(1:40) :: pre
REAL(kind=4)                  :: reso
INTEGER(kind=4)                :: tperr, i,j ,k, SK
REAL(kind=4), ALLOCATABLE, DIMENSION(:,:,:) :: tt_tmp

!*******
! 2. begin program

!===== 2.1. initialize

pre(1:40) = (/  0.1,   0.2,   0.5,   1.0,   1.5,   2.0,   3.0,   4.0,                         &
                5.0,   7.0,   10.,   20.,   30.,   50.,   70.,  100.,   125.,                 &
               150.,  175.,  200.,  225.,  250.,  275.,  300.,  350.,  400.,  450.,  500.,    &
               550.,  600.,  650.,  700.,  750.,  800.,  850.,  900.,                         &
               925.,  950.,  975., 1000./)
               
  ! check and open file 
  status = checkfile(TRIM(bin_name))
  IF (status /= 0) THEN
     CALL file_message('Numerical Weather Prediction binary data', status)
  ENDIF 
          
  ! OPEN the T639 binary data   
  nvar = 153+82
  OPEN(11,file=TRIM(bin_name),status='old',access='direct',form="binary",&
         recl=1440*720*nvar*4,IOSTAT=io_err)

  ! check the opened binary DATA 
  IF (io_err /=0) THEN    
     CALL bin_info_message('grapes gfs binary data')
  ENDIF
  
!===== 2.2. READ binary DATA
  ALLOCATE(bin_data(1440,720,nvar))
  READ(11,rec=1) (((bin_data(i,j,k),i=1,1440),j=1,720),k=1,nvar)
  CLOSE(11)
  
!===== 2.3. delete binary DATA
 ! CALL system('rm -rf '//TRIM(bin_name))
  
!===== 2.4. assign variables
!xxx lon lat
reso = 0.25
DO i=1,1440
DO j=1,720
   nwpo%lon(i,j) = -179.875+(i-1)*reso
   nwpo%lat(i,j) = 89.875 -(j-1)*reso
   !nwpo%lat(i,j) = -90.+(j-1)*reso
ENDDO
ENDDO

!xxx (1) plev_nointerp
nwpo%plev_nointerp(1:40,num) = pre

!xxx (2) psfc  ! the matrix transform is due to agree with longtitude distributions
DO j=1,720

SK=721-j   !revised by minmin 20200620   frist lat = -89.875 for grapes
!SK=j

nwpo%psfc(1:720,j,num)   = bin_data(721:1440,SK,1)
nwpo%psfc(721:1440,j,num) = bin_data(1:720,SK,1)

!xxx (3) pmsl 
nwpo%pmsl(1:720,j,num)   = bin_data(721:1440,SK,2)
nwpo%pmsl(721:1440,j,num) = bin_data(1:720,SK,2)

!xxx (4) tsfc
nwpo%tsfc(1:720,j,num)   = bin_data(721:1440,SK,3)
nwpo%tsfc(721:1440,j,num) = bin_data(1:720,SK,3)

!xxx (5) zsfc
nwpo%zsfc(1:720,j,num)   = bin_data(721:1440,SK,4)
nwpo%zsfc(721:1440,j,num) = bin_data(1:720,SK,4)

!xxx (6) albedo [sfc pre]
nwpo%albedo(1:720,j,num)   = bin_data(721:1440,SK,5)
nwpo%albedo(721:1440,j,num) = bin_data(1:720,SK,5)

!xxx (7) t_sigma
nwpo%t_sigma(1:720,j,num)   = bin_data(721:1440,SK,6)
nwpo%t_sigma(721:1440,j,num) = bin_data(1:720,SK,6)

!xxx (8) rh_sigma
nwpo%rh_sigma(1:720,j,num)  = bin_data(721:1440,SK,7)
nwpo%rh_sigma(721:1440,j,num)= bin_data(1:720,SK,7)

!xxx (9) u_sigma
nwpo%u_sigma(1:720,j,num)   = bin_data(721:1440,SK,8)
nwpo%u_sigma(721:1440,j,num) = bin_data(1:720,SK,8)

!xxx (10) v_sigma
nwpo%v_sigma(1:720,j,num)   = bin_data(721:1440,SK,9)
nwpo%v_sigma(721:1440,j,num) = bin_data(1:720,SK,9)

!xxx (11) tpw
nwpo%tpw(1:720,j,num)   = bin_data(721:1440,SK,10)
nwpo%tpw(721:1440,j,num) = bin_data(1:720,SK,10)

!xxx (12) weasd
nwpo%weasd(1:720,j,num)   = bin_data(721:1440,SK,11)
nwpo%weasd(721:1440,j,num) = bin_data(1:720,SK,11)

!xxx (13) o3col
nwpo%o3col(1:720,j,num)   = -999.0
nwpo%o3col(721:1440,j,num) = -999.0

!xxx (14) ttropo
nwpo%ttropo(1:720,j,num)    = bin_data(721:1440,SK,12)
nwpo%ttropo(721:1440,j,num)  = bin_data(1:720,SK,12)

!xxx (15) tlev
DO i=1,40
   nwpo%tlev(1:720,j,i,num)    = bin_data(721:1440,SK,12+i)
   nwpo%tlev(721:1440,j,i,num)  = bin_data(1:720,SK,12+i)  !last = 46
ENDDO

!xxx (16) zlev
DO i=1,40
   nwpo%zlev(1:720,j,i,num)    = bin_data(721:1440,SK,52+i)
   nwpo%zlev(721:1440,j,i,num)  = bin_data(1:720,SK,52+i)  !last = 82
ENDDO



!xxx (17) o3lev
nwpo%o3lev(:,j,:,num)   = -999.0
!DO i=1,36
!   nwpo%o3lev(1:640,:,i,num)   =   -999.0 ! bin_data(181:360,:,65+i)
!   nwpo%o3lev(641:1280,:,i,num) =  -999.0 ! bin_data(1:180,:,65+i)  !last = 71
!ENDDO

!xxx (18) rhlev
nwpo%rhlev(:,j,:,num)   = -999.0
DO i=1,30
   nwpo%rhlev(1:720,j,10+i,num)   = bin_data(721:1440,SK,92+i)
   nwpo%rhlev(721:1440,j,10+i,num) = bin_data(1:720,SK,92+i)  !last = 92
ENDDO

!if (j==10) then
!print*,nwpo%psfc(10,j,1), nwpo%pmsl(10,j,1),nwpo%tsfc(10,j,1),nwpo%zsfc(10,j,1)
!print*,nwpo%albedo(10,j,1),nwpo%t_sigma(10,j,1),nwpo%rh_sigma(10,j,1)
!print*,nwpo%u_sigma(10,j,1),nwpo%v_sigma(10,j,1)
!print*,nwpo%tpw(10,j,1),nwpo%weasd(10,j,1), nwpo%ttropo(10,j,1)
!print*,'test'
!print*,'t = ',nwpo%tlev(10,j,:,1)
!print*,'z = ',nwpo%zlev(10,j,:,1)
!print*,'rh = ',nwpo%rhlev(10,j,:,1)
!stop
!endif

!xxx (19) clwlev
nwpo%clwlev(:,j,:,num)   = -999.0
DO i=1,30
   nwpo%clwlev(1:720,j,10+i,num)  =    bin_data(721:1440,SK,122+i)
   nwpo%clwlev(721:1440,j,10+i,num)=   bin_data(1:720,SK,122+i)  !last = 113
ENDDO

!xxx ulev
DO i=1,40
   nwpo%ulev(1:720,j,i,num)    = bin_data(721:1440,SK,152+i)
   nwpo%ulev(721:1440,j,i,num)  = bin_data(1:720,SK,152+i)  !
ENDDO

!xxx vlev
DO i=1,40
   nwpo%vlev(1:720,j,i,num)    = bin_data(721:1440,SK,192+i)
   nwpo%vlev(721:1440,j,i,num)  = bin_data(1:720,SK,192+i)  !
ENDDO

!xxx u10m
!nwpo%u10m(1:720,j,num)         = bin_data(721:1440,SK,232+1)
!nwpo%u10m(721:1440,j,num)      = bin_data(1:720,SK,232+1)

!xxx v10m
!nwpo%v10m(1:720,j,num)         = bin_data(721:1440,SK,232+2)
!nwpo%v10m(721:1440,j,num)      = bin_data(1:720,SK,232+2)

ENDDO

! re-calculate tropopause temperature 2003 GRL Thomas
   !nwpo%ttropo(i,j,num)
allocate(tt_tmp(1440,720,28))
tt_tmp(:,:,:) = nwpo%tlev(:,:,13:40,num)
CALL tropo( tt_tmp(:,:,:) ,  &
            nwp40%nlon,                &
            nwp40%nlat,                &
            28,                        &
            pre(13:40),                &
            45000., 7500., 7500.,      &
            .true.,                    &
            nwpo%ttropo(:,:,num),      &  
            tperr)
deallocate(tt_tmp) 
            
DEALLOCATE(bin_data)
! 3. END SUBROUTINE   
END SUBROUTINE read_NWP_grapes_gfs_arrays
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~ SUBROUTINE 3: READ NWP DATA with binary format ~~~~
SUBROUTINE read_NWP_arrays_0p25(bin_name, num)

!-----------------------------------------------------------------------
! !F90 read_NWP_arrays
!
! !Description:
!    This program is to READ NWP arrays in binary format and make the 
!    necessary nwp%DATA for next calculations.
!
! !Input  parameters:
!    bin_name        = NWP DATA's name
!    num             = number of NWP DATA [1 or 2]
!
! !Output parameters:
!    none
!
!-----------------------------------------------------------------------

! USE modules

! 1. define variables
!===== 1.1.input information binary DATA's name
CHARACTER(*), INTENT(in)  :: bin_name
INTEGER(KIND=4), INTENT(in)       :: num

!===== 1.2. other
INTEGER(KIND=4) :: status, io_err, SK,CK

!===== 1.3. variables 
!REAL(KIND=real4), DIMENSION(1:360,1:181,1:124) :: bin_data
REAL(KIND=4), ALLOCATABLE, DIMENSION(:,:,:) :: bin_data
REAL(KIND=4), DIMENSION(1:31) :: pre

DATA pre / 1.0,   2.0,    3.0,   5.0,   7.0,                             &
          10.0,  20.0,   30.0,  50.0,  70.0, 100.0, 150.0, 200.0, 250.0, &
          300.0, 350.0, 400.0, 450.0, 500.0, 550.0, 600.0, 650.0, 700.0, &
          750.0, 800.0, 850.0, 900.0, 925.0, 950.0, 975.0, 1000.0/

!*******
! 2. begin program

!===== 2.1. initialize
  ! check and open file 
  status = checkfile(TRIM(bin_name))
  IF (status /= 0) THEN
     CALL file_message('Numerical Weather Prediction binary data', status)
  ENDIF 
          
  ! OPEN the ncep binary DATA       
  IF (fylat_nwp_opt == 7 .or. fylat_nwp_opt == 8) THEN  ! grib2
     !OPEN(11,file=TRIM(bin_name),status='old',access='direct',form="binary",&
     !     recl=361*720*124*4,IOSTAT=io_err)     
     OPEN(11,file=TRIM(bin_name),status='old',access='direct',form="binary",&
          recl=721*1440*(124+52+30)*4,IOSTAT=io_err)     
  ENDIF
  print*,'read nwp bin data ok!'
         
  ! check the opened binary DATA 
  IF (io_err /= 0) THEN    
     CALL bin_info_message('NCEP binary DATA')
  ENDIF
  
  IF (fylat_nwp_opt == 7 .or. fylat_nwp_opt == 8) THEN   !gfs0p50

     !ALLOCATE(bin_data(720,361,124))
     ALLOCATE(bin_data(1440,721,124+52+30))
     !READ(11,rec=1) (((bin_data(i,j,k),i=1,720),j=1,361),k=1,124)
     READ(11,rec=1) (((bin_data(i,j,k),i=1,1440),j=1,721),k=1,124+52+30)
     CLOSE(11) 
   
  ENDIF
    
!===== 2.2. READ binary DATA
IF (fylat_nwp_opt == 7 .or. fylat_nwp_opt == 8) THEN  !gfs0p25

!xxx lon lat
DO i=1,1440
DO j=1,721
   nwpo%lon(i,j) = -179.875+(i-1)*0.25
   nwpo%lat(i,j) = 90.-(j-1)*0.25
ENDDO
ENDDO

!xxx (1) plev_nointerp
nwpo%plev_nointerp(1:31,num) = pre

DO j=1,721
!xxx (2) psfc  ! the matrix transform is due to agree with longtitude distributions

! gfs0p50 grib2
SK=722-j

nwpo%psfc(1:720,j,num)   = bin_data(721:1440,SK,1)
nwpo%psfc(721:1440,j,num) = bin_data(1:720,SK,1)

!xxx (2) pmsl 
nwpo%pmsl(1:720,j,num)   = bin_data(721:1440,SK,2)
nwpo%pmsl(721:1440,j,num) = bin_data(1:720,SK,2)

!xxx (3) tsfc
nwpo%tsfc(1:720,j,num)   = bin_data(721:1440,SK,3)
nwpo%tsfc(721:1440,j,num) = bin_data(1:720,SK,3)

!xxx (4) zsfc
nwpo%zsfc(1:720,j,num)   = bin_data(721:1440,SK,4)
nwpo%zsfc(721:1440,j,num) = bin_data(1:720,SK,4)

!xxx (5) albedo
nwpo%albedo(1:720,j,num)   = bin_data(721:1440,SK,5)
nwpo%albedo(721:1440,j,num) = bin_data(1:720,SK,5)

!xxx (6) t_sigma
nwpo%t_sigma(1:720,j,num)   = bin_data(721:1440,SK,6)
nwpo%t_sigma(721:1440,j,num) = bin_data(1:720,SK,6)

!xxx (7) rh_sigma
nwpo%rh_sigma(1:720,j,num)  = bin_data(721:1440,SK,7)
nwpo%rh_sigma(721:1440,j,num)= bin_data(1:720,SK,7)

!xxx (8) u_sigma
nwpo%u_sigma(1:720,j,num)   = bin_data(721:1440,SK,8)
nwpo%u_sigma(721:1440,j,num) = bin_data(1:720,SK,8)

!xxx (9) v_sigma
nwpo%v_sigma(1:720,j,num)   = bin_data(721:1440,SK,9)
nwpo%v_sigma(721:1440,j,num) = bin_data(1:720,SK,9)

!xxx (10) tpw
nwpo%tpw(1:720,j,num)   = bin_data(721:1440,SK,10)
nwpo%tpw(721:1440,j,num) = bin_data(1:720,SK,10)

!xxx (11) weasd
nwpo%weasd(1:720,j,num)   = bin_data(721:1440,SK,11)
nwpo%weasd(721:1440,j,num) = bin_data(1:720,SK,11)

!xxx (12) o3col
nwpo%o3col(1:720,j,num)   = bin_data(721:1440,SK,12)
nwpo%o3col(721:1440,j,num) = bin_data(1:720,SK,12)

!xxx (13) ttropo
nwpo%ttropo(1:720,j,num)    = bin_data(721:1440,SK,13)
nwpo%ttropo(721:1440,j,num)  = bin_data(1:720,SK,13)

!xxx (14) tlev
DO i=1,31
   nwpo%tlev(1:720,j,i,num)    = bin_data(721:1440,SK,13+i)
   nwpo%tlev(721:1440,j,i,num)  = bin_data(1:720,SK,13+i)  !last = 44
ENDDO

!xxx (15) zlev
DO i=1,31
   nwpo%zlev(1:720,j,i,num)    = bin_data(721:1440,SK,44+i)
   nwpo%zlev(721:1440,j,i,num)  = bin_data(1:720,SK,44+i)  !last = 75
ENDDO

!xxx (16) o3lev
nwpo%o3lev(:,j,:,num)   = -999.0
DO i=1,17
   nwpo%o3lev(1:720,j,i,num)   = bin_data(721:1440,SK,75+i)
   nwpo%o3lev(721:1440,j,i,num) = bin_data(1:720,SK,75+i)  !last = 92
ENDDO

!xxx (18) rhlev
nwpo%rhlev(:,j,:,num)   = -999.0
DO i=1,31
   nwpo%rhlev(1:720,j,i,num)   = bin_data(721:1440,SK,92+i)
   nwpo%rhlev(721:1440,j,i,num) = bin_data(1:720,SK,92+i)  !last = 123
ENDDO


!xxx (19) clwlev
nwpo%clwlev(:,j,:,num)   = -999.0
DO i=1,21
   nwpo%clwlev(1:720,j,i+10,num)  = bin_data(721:1440,SK,123+i)
   nwpo%clwlev(721:1440,j,i+10,num)= bin_data(1:720,SK,123+i)  !last = 144
ENDDO

!xxx ulev
DO i=1,31
   nwpo%ulev(1:720,j,i,num)    = bin_data(721:1440,SK,144+i)
   nwpo%ulev(721:1440,j,i,num)  = bin_data(1:720,SK,144+i)  ! !last = 175
ENDDO

!xxx vlev
DO i=1,31
   nwpo%vlev(1:720,j,i,num)    = bin_data(721:1440,SK,175+i)
   nwpo%vlev(721:1440,j,i,num)  = bin_data(1:720,SK,175+i)  !
ENDDO

ENDDO

ENDIF ! IF (fylat_nwp_opt == 8) THEN

DEALLOCATE(bin_data)

! 3. END SUBROUTINE   
END SUBROUTINE read_NWP_arrays_0p25
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~ SUBROUTINE 3: READ NWP DATA with binary format ~~~~
SUBROUTINE read_NWP_arrays_0p25_41Layers(bin_name, num)

!-----------------------------------------------------------------------
! !F90 read_NWP_arrays
!
! !Description:
!    This program is to READ NWP arrays in binary format and make the 
!    necessary nwp%DATA for next calculations.
!
! !Input  parameters:
!    bin_name        = NWP DATA's name
!    num             = number of NWP DATA [1 or 2]
!
! !Output parameters:
!    none
!
!-----------------------------------------------------------------------

! USE modules

! 1. define variables
!===== 1.1.input information binary DATA's name
CHARACTER(*), INTENT(in)  :: bin_name
INTEGER(KIND=4), INTENT(in)       :: num

!===== 1.2. other
INTEGER(KIND=4) :: status, io_err, SK,CK

!===== 1.3. variables 
!REAL(KIND=real4), DIMENSION(1:360,1:181,1:124) :: bin_data
REAL(KIND=4), ALLOCATABLE, DIMENSION(:,:,:) :: bin_data
REAL(KIND=4), DIMENSION(1:41) :: pre

DATA pre / 0.01,  0.02,   0.04, 0.07,   0.1,   0.2,   0.4,   0.7,        &
           1.0,   2.0,    3.0,   5.0,   7.0,                             &
          10.0,   15.0,  20.0,  30.0,  40.0,  50.0,  70.0,               &
          100.0, 150.0, 200.0, 250.0, &
          300.0, 350.0, 400.0, 450.0, 500.0, 550.0, 600.0, 650.0, 700.0, &
          750.0, 800.0, 850.0, 900.0, 925.0, 950.0, 975.0, 1000.0/
         
!*******
! 2. begin program

!===== 2.1. initialize
  ! check and open file 
  status = checkfile(TRIM(bin_name))
  IF (status /= 0) THEN
     CALL file_message('Numerical Weather Prediction binary data', status)
  ENDIF 
          
  ! OPEN the ncep binary DATA       
  IF (fylat_nwp_opt == 10) THEN  ! grib2
     !OPEN(11,file=TRIM(bin_name),status='old',access='direct',form="binary",&
     !     recl=361*720*124*4,IOSTAT=io_err)     
     OPEN(11,file=TRIM(bin_name),status='old',access='direct',form="binary",&
          recl=721*1440*(177+22+82)*4,IOSTAT=io_err)     
  ENDIF
  print*,'read nwp bin data ok!'
         
  ! check the opened binary DATA 
  IF (io_err /= 0) THEN    
     CALL bin_info_message('NCEP binary DATA')
  ENDIF
  
  IF (fylat_nwp_opt == 10) THEN   !gfs0p50

     !ALLOCATE(bin_data(720,361,124))
     ALLOCATE(bin_data(1440,721,177+22+82))
     !READ(11,rec=1) (((bin_data(i,j,k),i=1,720),j=1,361),k=1,124)
     READ(11,rec=1) (((bin_data(i,j,k),i=1,1440),j=1,721),k=1,177+22+82)
     CLOSE(11) 
   
  ENDIF
    
!===== 2.2. READ binary DATA
IF (fylat_nwp_opt == 10) THEN  !gfs0p25

!xxx lon lat
DO i=1,1440
DO j=1,721
   nwpo%lon(i,j) = -179.875+(i-1)*0.25
   nwpo%lat(i,j) = 90.-(j-1)*0.25
ENDDO
ENDDO

!xxx (1) plev_nointerp
nwpo%plev_nointerp(1:41,num) = pre

DO j=1,721
!xxx (2) psfc  ! the matrix transform is due to agree with longtitude distributions

! gfs0p50 grib2
SK=722-j

nwpo%psfc(1:720,j,num)   = bin_data(721:1440,SK,1)
nwpo%psfc(721:1440,j,num) = bin_data(1:720,SK,1)

!xxx (2) pmsl 
nwpo%pmsl(1:720,j,num)   = bin_data(721:1440,SK,2)
nwpo%pmsl(721:1440,j,num) = bin_data(1:720,SK,2)

!xxx (3) tsfc
nwpo%tsfc(1:720,j,num)   = bin_data(721:1440,SK,3)
nwpo%tsfc(721:1440,j,num) = bin_data(1:720,SK,3)

!xxx (4) zsfc
nwpo%zsfc(1:720,j,num)   = bin_data(721:1440,SK,4)
nwpo%zsfc(721:1440,j,num) = bin_data(1:720,SK,4)

!xxx (5) albedo
nwpo%albedo(1:720,j,num)   = bin_data(721:1440,SK,5)
nwpo%albedo(721:1440,j,num) = bin_data(1:720,SK,5)

!xxx (6) t_sigma
nwpo%t_sigma(1:720,j,num)   = bin_data(721:1440,SK,6)
nwpo%t_sigma(721:1440,j,num) = bin_data(1:720,SK,6)

!xxx (7) rh_sigma
nwpo%rh_sigma(1:720,j,num)  = bin_data(721:1440,SK,7)
nwpo%rh_sigma(721:1440,j,num)= bin_data(1:720,SK,7)

!xxx (8) u_sigma
nwpo%u_sigma(1:720,j,num)   = bin_data(721:1440,SK,8)
nwpo%u_sigma(721:1440,j,num) = bin_data(1:720,SK,8)

!xxx (9) v_sigma
nwpo%v_sigma(1:720,j,num)   = bin_data(721:1440,SK,9)
nwpo%v_sigma(721:1440,j,num) = bin_data(1:720,SK,9)

!xxx (10) tpw
nwpo%tpw(1:720,j,num)   = bin_data(721:1440,SK,10)
nwpo%tpw(721:1440,j,num) = bin_data(1:720,SK,10)

!xxx (11) weasd
nwpo%weasd(1:720,j,num)   = bin_data(721:1440,SK,11)
nwpo%weasd(721:1440,j,num) = bin_data(1:720,SK,11)

!xxx (12) o3col
nwpo%o3col(1:720,j,num)   = bin_data(721:1440,SK,12)
nwpo%o3col(721:1440,j,num) = bin_data(1:720,SK,12)

!xxx (13) ttropo
nwpo%ttropo(1:720,j,num)    = bin_data(721:1440,SK,13)
nwpo%ttropo(721:1440,j,num)  = bin_data(1:720,SK,13)

!xxx (14) tlev
DO i=1,41
   nwpo%tlev(1:720,j,i,num)    = bin_data(721:1440,SK,13+i)
   nwpo%tlev(721:1440,j,i,num)  = bin_data(1:720,SK,13+i)  !last = 54
ENDDO

!xxx (15) zlev
DO i=1,41
   nwpo%zlev(1:720,j,i,num)    = bin_data(721:1440,SK,54+i)
   nwpo%zlev(721:1440,j,i,num)  = bin_data(1:720,SK,54+i)  !last = 95
ENDDO

!xxx (16) o3lev
nwpo%o3lev(:,j,:,num)   = -999.0
DO i=1,41
   nwpo%o3lev(1:720,j,i,num)   = bin_data(721:1440,SK,95+i)
   nwpo%o3lev(721:1440,j,i,num) = bin_data(1:720,SK,95+i)  !last = 136
ENDDO

!xxx (18) rhlev
nwpo%rhlev(:,j,:,num)   = -999.0
DO i=1,41
   nwpo%rhlev(1:720,j,i,num)   = bin_data(721:1440,SK,136+i)
   nwpo%rhlev(721:1440,j,i,num) = bin_data(1:720,SK,136+i)  !last = 177
ENDDO


!xxx (19) clwlev
nwpo%clwlev(:,j,:,num)   = -999.0
DO i=1,22
   nwpo%clwlev(1:720,j,i+18,num)  = bin_data(721:1440,SK,177+i)
   nwpo%clwlev(721:1440,j,i+18,num)= bin_data(1:720,SK,177+i)  !last = 199
ENDDO

!xxx ulev
DO i=1,41
   nwpo%ulev(1:720,j,i,num)    = bin_data(721:1440,SK,199+i)
   nwpo%ulev(721:1440,j,i,num)  = bin_data(1:720,SK,199+i)  ! !last = 240
ENDDO

!xxx vlev
DO i=1,41
   nwpo%vlev(1:720,j,i,num)    = bin_data(721:1440,SK,240+i)
   nwpo%vlev(721:1440,j,i,num)  = bin_data(1:720,SK,240+i)  !
ENDDO

ENDDO

ENDIF ! IF (fylat_nwp_opt == 8) THEN

DEALLOCATE(bin_data)

! 3. END SUBROUTINE   
END SUBROUTINE read_NWP_arrays_0p25_41Layers
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~ SUBROUTINE 3: READ NWP DATA with binary format ~~~~
SUBROUTINE read_NWP_arrays_0p50_41Layers(bin_name, num)

!-----------------------------------------------------------------------
! !F90 read_NWP_arrays
!
! !Description:
!    This program is to READ NWP arrays in binary format and make the 
!    necessary nwp%DATA for next calculations.
!
! !Input  parameters:
!    bin_name        = NWP DATA's name
!    num             = number of NWP DATA [1 or 2]
!
! !Output parameters:
!    none
!
!-----------------------------------------------------------------------

! USE modules

! 1. define variables
!===== 1.1.input information binary DATA's name
CHARACTER(*), INTENT(in)  :: bin_name
INTEGER(KIND=4), INTENT(in)       :: num

!===== 1.2. other
INTEGER(KIND=4) :: status, io_err, SK,CK

!===== 1.3. variables 
!REAL(KIND=real4), DIMENSION(1:360,1:181,1:124) :: bin_data
REAL(KIND=4), ALLOCATABLE, DIMENSION(:,:,:) :: bin_data
REAL(KIND=4), DIMENSION(1:41) :: pre

DATA pre / 0.01,  0.02,   0.04, 0.07,   0.1,   0.2,   0.4,   0.7,        &
           1.0,   2.0,    3.0,   5.0,   7.0,                             &
          10.0,   15.0,  20.0,  30.0,  40.0,  50.0,  70.0,               &
          100.0, 150.0, 200.0, 250.0, &
          300.0, 350.0, 400.0, 450.0, 500.0, 550.0, 600.0, 650.0, 700.0, &
          750.0, 800.0, 850.0, 900.0, 925.0, 950.0, 975.0, 1000.0/
         
!*******
! 2. begin program

!===== 2.1. initialize
  ! check and open file 
  status = checkfile(TRIM(bin_name))
  IF (status /= 0) THEN
     CALL file_message('Numerical Weather Prediction binary data', status)
  ENDIF 
          
  ! OPEN the ncep binary DATA       
  IF (fylat_nwp_opt == 9) THEN  ! grib2
     !OPEN(11,file=TRIM(bin_name),status='old',access='direct',form="binary",&
     !     recl=361*720*124*4,IOSTAT=io_err)     
     OPEN(11,file=TRIM(bin_name),status='old',access='direct',form="binary",&
          recl=361*720*(177+22+82)*4,IOSTAT=io_err)     
  ENDIF
  print*,'read nwp bin data ok!'
         
  ! check the opened binary DATA 
  IF (io_err /= 0) THEN    
     CALL bin_info_message('NCEP binary DATA')
  ENDIF
  
  IF (fylat_nwp_opt == 9) THEN   !gfs0p50

     !ALLOCATE(bin_data(720,361,124))
     ALLOCATE(bin_data(720,361,177+22+82))
     !READ(11,rec=1) (((bin_data(i,j,k),i=1,720),j=1,361),k=1,124)
     READ(11,rec=1) (((bin_data(i,j,k),i=1,720),j=1,361),k=1,177+22+82)
     CLOSE(11) 
   
  ENDIF
    
!===== 2.2. READ binary DATA
IF (fylat_nwp_opt == 9) THEN  !gfs0p25

!xxx lon lat
do i=1,720
do j=1,361
   nwpo%lon(i,j) = -179.75+(i-1)*0.5
   nwpo%lat(i,j) = 90.-(j-1)*0.5
enddo
enddo

!xxx (1) plev_nointerp
nwpo%plev_nointerp(1:41,num) = pre

DO j=1,361
!xxx (2) psfc  ! the matrix transform is due to agree with longtitude distributions

! gfs0p50 grib2
SK=362-j

nwpo%psfc(1:360,j,num)   = bin_data(361:720,SK,1)
nwpo%psfc(361:720,j,num) = bin_data(1:360,SK,1)

!xxx (2) pmsl 
nwpo%pmsl(1:360,j,num)   = bin_data(361:720,SK,2)
nwpo%pmsl(361:720,j,num) = bin_data(1:360,SK,2)

!xxx (3) tsfc
nwpo%tsfc(1:360,j,num)   = bin_data(361:720,SK,3)
nwpo%tsfc(361:720,j,num) = bin_data(1:360,SK,3)

!xxx (4) zsfc
nwpo%zsfc(1:360,j,num)   = bin_data(361:720,SK,4)
nwpo%zsfc(361:720,j,num) = bin_data(1:360,SK,4)

!xxx (5) albedo
nwpo%albedo(1:360,j,num)   = bin_data(361:720,SK,5)
nwpo%albedo(361:720,j,num) = bin_data(1:360,SK,5)

!xxx (6) t_sigma
nwpo%t_sigma(1:360,j,num)   = bin_data(361:720,SK,6)
nwpo%t_sigma(361:720,j,num) = bin_data(1:360,SK,6)

!xxx (7) rh_sigma
nwpo%rh_sigma(1:360,j,num)  = bin_data(361:720,SK,7)
nwpo%rh_sigma(361:720,j,num)= bin_data(1:360,SK,7)

!xxx (8) u_sigma
nwpo%u_sigma(1:360,j,num)   = bin_data(361:720,SK,8)
nwpo%u_sigma(361:720,j,num) = bin_data(1:360,SK,8)

!xxx (9) v_sigma
nwpo%v_sigma(1:360,j,num)   = bin_data(361:720,SK,9)
nwpo%v_sigma(361:720,j,num) = bin_data(1:360,SK,9)

!xxx (10) tpw
nwpo%tpw(1:360,j,num)   = bin_data(361:720,SK,10)
nwpo%tpw(361:720,j,num) = bin_data(1:360,SK,10)

!xxx (11) weasd
nwpo%weasd(1:360,j,num)   = bin_data(361:720,SK,11)
nwpo%weasd(361:720,j,num) = bin_data(1:360,SK,11)

!xxx (12) o3col
nwpo%o3col(1:360,j,num)   = bin_data(361:720,SK,12)
nwpo%o3col(361:720,j,num) = bin_data(1:360,SK,12)

!xxx (13) ttropo
nwpo%ttropo(1:360,j,num)    = bin_data(361:720,SK,13)
nwpo%ttropo(361:720,j,num)  = bin_data(1:360,SK,13)

!xxx (14) tlev
DO i=1,41
   nwpo%tlev(1:360,j,i,num)    = bin_data(361:720,SK,13+i)
   nwpo%tlev(361:720,j,i,num)  = bin_data(1:360,SK,13+i)  !last = 54
ENDDO

!xxx (15) zlev
DO i=1,41
   nwpo%zlev(1:360,j,i,num)    = bin_data(361:720,SK,54+i)
   nwpo%zlev(361:720,j,i,num)  = bin_data(1:360,SK,54+i)  !last = 95
ENDDO

!xxx (16) o3lev
nwpo%o3lev(:,j,:,num)   = -999.0
DO i=1,41
   nwpo%o3lev(1:360,j,i,num)   = bin_data(361:720,SK,95+i)
   nwpo%o3lev(361:720,j,i,num) = bin_data(1:360,SK,95+i)  !last = 136
ENDDO

!xxx (18) rhlev
nwpo%rhlev(:,j,:,num)   = -999.0
DO i=1,41
   nwpo%rhlev(1:360,j,i,num)   = bin_data(361:720,SK,136+i)
   nwpo%rhlev(361:720,j,i,num) = bin_data(1:360,SK,136+i)  !last = 177
ENDDO


!xxx (19) clwlev
nwpo%clwlev(:,j,:,num)   = -999.0
DO i=1,22
   nwpo%clwlev(1:360,j,i+18,num)  = bin_data(361:720,SK,177+i)
   nwpo%clwlev(361:720,j,i+18,num)= bin_data(1:360,SK,177+i)  !last = 199
ENDDO

!xxx ulev
DO i=1,41
   nwpo%ulev(1:360,j,i,num)    = bin_data(361:720,SK,199+i)
   nwpo%ulev(361:720,j,i,num)  = bin_data(1:360,SK,199+i)  ! !last = 240
ENDDO

!xxx vlev
DO i=1,41
   nwpo%vlev(1:360,j,i,num)    = bin_data(361:720,SK,240+i)
   nwpo%vlev(361:720,j,i,num)  = bin_data(1:360,SK,240+i)  !
ENDDO

ENDDO

ENDIF ! IF (fylat_nwp_opt == 8) THEN

DEALLOCATE(bin_data)

! 3. END SUBROUTINE   
END SUBROUTINE read_NWP_arrays_0p50_41Layers
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!subroutine tropo(temp, nlon, nlat, nlev, pres, plimu, pliml, plimlex, dofill, tp, tperr)
subroutine tropo(temp, nlon, nlat, nlev, pres, plimu, pliml, plimlex, dofill, tt, tperr)

!-----------------------------------------------------------------------
!
! determination of tropopause height from gridded temperature data
!
! reference:  Reichler, T., M. Dameris, R. Sausen (2003): 
!             Determining the tropopause height from gridded data, 
!             Geophys. Res. L., 30, No. 20, 2042
!
! input:    temp(nlon,nlat,nlev)    3D-temperature field
!           nlon                    # of grid points in x
!           nlat                    # of grid points in y
!           nlev                    # of vertical pressure levels
!           pres(nlev)              array of pressure levels in hPa, length = nlev
!           plimu                   upper limit for tropopause pressure in Pa, usually 45000
!           pliml                   lower limit for tropopause pressure in Pa, usually 7500
!           plimlex                 lower limit in extratropics, usually same as pliml, i.e., 7500
!           dofill                  fill undefined values with neighboring points if .true.
!
! output:   tp(nlon, nlat)          tropopause pressure in Pa, same horizontal dimension as temp
!           tperr                   # of undetermined values
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

implicit none

integer(kind=4),intent(in)                         :: nlon, nlat, nlev
real(kind=4),intent(in),dimension(nlon,nlat,nlev) :: temp
real(kind=4),intent(in),dimension(nlev)           :: pres
real(kind=4), intent(in)                          :: plimu, pliml, plimlex
logical, intent(in)                                   :: dofill
real(kind=4),dimension(nlon,nlat)     :: tp
integer(kind=4),intent(out)                        :: tperr
real(kind=4),intent(out),dimension(nlon,nlat)     :: tt

integer(kind=4)                                    :: i, j, invert, ifil
integer(kind=4)                                    :: lon, lat
real(kind=4),dimension(nlev)                      :: t
real(kind=4),dimension(nlev)                      :: p
real(kind=4)                                      :: trp
real(kind=4)                                      :: ptemp, t1, t2, w1, w2, p1, p2, dp , dt

real(kind=4), parameter                           :: gamma=-0.002 ! K/m

! check vertical orientation of data
if (pres(1) .gt. pres(2)) then
   invert=1
   do i=1,nlev
      p(i)=pres(nlev+1-i)*100.  ! hPa > Pa
   enddo
else
   invert=0
   do i=1,nlev
      p(i)=pres(i)*100.         ! hPa > Pa
   enddo
endif

tperr = 0
do lon=1,nlon
do lat=1,nlat
   if (invert.eq.1) then
      do i=1,nlev
         t(i)=temp(lon,lat,nlev+1-i)
      enddo
   else
      do i=1,nlev
         t(i)=temp(lon,lat,i)
      enddo
   endif
   call twmo(nlev, t, p, plimu, pliml, gamma, trp)
   if (lat.lt..15*nlat.and.trp.lt.plimlex) trp=-99.
   if (lat.gt..85*nlat.and.trp.lt.plimlex) trp=-99.
   tp(lon,lat)=trp
   if (trp.lt..0) then
      tperr = tperr+1    
   endif
end do
end do

! fill holes
if (dofill) then
   call fill(tp, nlon, nlat, ifil)
   if (ifil.ne.tperr) then
      print*, 'Inconsistent'
      stop
   endif
endif

! get tropopause temperature
t1 = 0.
t2 = 0.
p1 = 0.
p2 = 0.
do lon=1,nlon
do lat=1,nlat
   ptemp = tp(lon,lat)
   do i=1,nlev
      t(i)=temp(lon,lat,nlev+1-i)
      p(i)=pres(i)*100.         ! hPa > Pa
   enddo

   do i=1,nlev-1
      if (ptemp >= p(i) .and. ptemp < p(i+1)) then
         p1 = p(i)
         p2 = p(i+1)
         t1 = t(i)
         t2 = t(i+1)
         dp = log(p2/p1)
         dt = (t2-t1)/dp
         exit
      endif
   enddo
   tt(lon,lat)= t1 + dt * log(ptemp/p1)
end do
end do

return
END SUBROUTINE tropo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! twmo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine twmo(level, t, p, plimu, pliml, gamma, trp)

implicit none
integer(kind=4),intent(in)                   :: level
real(kind=4),intent(in),dimension(level)    :: t, p
real(kind=4),intent(in)                     :: plimu, pliml, gamma
real(kind=4),intent(out)                    :: trp

real(kind=4),parameter                      :: kap=0.286
real(kind=4),parameter                      :: faktor = -9.81/287.0
real(kind=4),parameter                      :: deltaz = 2000.0
real(kind=4),parameter                      :: ka1=kap-1.

real(kind=4)                                :: pmk, pm, a, b, tm, dtdp, dtdz
real(kind=4)                                :: ag, bg, ptph
real(kind=4)                                :: pm0, pmk0, dtdz0
real(kind=4)                                :: p2km, asum, aquer
real(kind=4)                                :: pmk2, pm2, a2, b2, tm2, dtdp2, dtdz2
integer(kind=4)                              :: icount, jj
integer(kind=4)                              :: j

trp=-99.0                           ! negative means not valid
do j=level,2,-1

   ! dt/dz
   pmk= .5 * (p(j-1)**kap+p(j)**kap)
   pm = pmk**(1/kap)              
   a = (t(j-1)-t(j))/(p(j-1)**kap-p(j)**kap)
   b = t(j)-(a*p(j)**kap)
   tm = a * pmk + b              
   dtdp = a * kap * (pm**ka1)
   dtdz = faktor*dtdp*pm/tm

   ! dt/dz valid?
   if (j.eq.level)    go to 999     ! no, start level, initialize first
   if (dtdz.le.gamma) go to 999     ! no, dt/dz < -2 K/km
   if (pm.gt.plimu)   go to 999     ! no, too low

   ! dtdz is valid, calculate tropopause pressure
   if (dtdz0.lt.gamma) then
      ag = (dtdz-dtdz0) / (pmk-pmk0)    
      bg = dtdz0 - (ag * pmk0)         
      ptph = exp(log((gamma-bg)/ag)/kap)
   else
      ptph = pm
   endif

   if (ptph.lt.pliml) go to 999    
   if (ptph.gt.plimu) go to 999          

   ! 2nd test: dtdz above 2 km must not exceed gamma
   p2km = ptph + deltaz*(pm/tm)*faktor          ! p at ptph + 2km
   asum = 0.0                                   ! dtdz above
   icount = 0                                   ! number of levels above

   ! test until apm < p2km
   do jj=j,2,-1

       pmk2 = .5 * (p(jj-1)**kap+p(jj)**kap)    ! p mean ^kappa
       pm2 = pmk2**(1/kap)                      ! p mean
       if(pm2.gt.ptph) go to 110                ! doesn't happen
       if(pm2.lt.p2km) go to 888                ! ptropo is valid

       a2 = (t(jj-1)-t(jj))                     ! a
       a2 = a2/(p(jj-1)**kap-p(jj)**kap)
       b2 = t(jj)-(a2*p(jj)**kap)               ! b
       tm2 = a2 * pmk2 + b2                     ! T mean
       dtdp2 = a2 * kap * (pm2**(kap-1))        ! dt/dp
       dtdz2 = faktor*dtdp2*pm2/tm2
       asum = asum+dtdz2
       icount = icount+1
       aquer = asum/float(icount)               ! dt/dz mean
  
       ! discard ptropo ?
        if (aquer.le.gamma) go to 999           ! dt/dz above < gamma

110 continue
    enddo                           ! test next level

888 continue                        ! ptph is valid
    trp = ptph
    return

999 continue                        ! continue search at next higher level
    pm0 = pm
    pmk0 = pmk
    dtdz0  = dtdz

enddo

! no tropopouse found
return
end subroutine twmo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fill(dat, ix, iy, ir)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer(kind=4), intent(in)                 :: ix, iy
integer(kind=4), intent(out)                :: ir
real(kind=4), dimension(ix,iy)              :: dat
real(kind=4), dimension(4)                  :: help
integer(kind=4)                             :: jx, jy, icount, ipk, ic, ii, jj
real(kind=4)                                :: drop, sum

icount = 0
do jx=1,ix
do jy=1,iy
   if (loch(dat(jx,jy))) icount = icount+1
enddo
enddo
if (icount.gt.(ix*iy)/2) stop 'ERROR: Too many holes (>50%)'
ir = icount
if (icount.eq.0) return

ipk = 0
10   continue
do jx=1,ix
do jy=1,iy
     if(loch(dat(jx,jy))) then
          drop = dat(jx,jy)

        ! left edge
          if (jx.eq.1) then
          if (jy.eq.1) then
          help(1) = dat(jx,jy+1)
          help(2) = dat(jx+1,jy)
          help(3) = drop
          help(4) = drop
          go to 200
          endif
          if (jy.eq.iy) then
          help(1) = drop
          help(2) = dat(jx+1,jy)
          help(3) = dat(jx,jy-1)
          help(4) = drop
          go to 200
          endif
          help(1) = dat(jx,jy+1)
          help(2) = dat(jx+1,jy)
          help(3) = dat(jx,jy-1)
          help(4) = drop
          go to 200
          endif

          ! right edge
          if (jx.eq.ix) then
          if (jy.eq.1) then
          help(1) = dat(jx,jy+1)
          help(2) = drop
          help(3) = drop
          help(4) = dat(jx-1,jy)
          go to 200
          endif
          if (jy.eq.iy) then
          help(1) = drop
          help(2) = drop
          help(3) = dat(jx,jy-1)
          help(4) = dat(jx-1,jy)
          go to 200
          endif
          help(1) = dat(jx,jy+1)
          help(2) = drop
          help(3) = dat(jx,jy-1)
          help(4) = dat(jx-1,jy)
          go to 200
          endif

        ! bottom edge
          if (jy.eq.1) then
          help(1) = dat(jx,jy+1)
          help(2) = dat(jx+1,jy)
          help(3) = drop
          help(4) = dat(jx-1,jy)
          go to 200
          endif

          ! upper edge
          if(jy.eq.iy) then
          help(1) = drop
          help(2) = dat(jx+1,jy)
          help(3) = dat(jx,jy-1)
          help(4) = dat(jx-1,jy)
          go to 200
          endif

        ! no edge
          help(1) = dat(jx,jy+1)
          help(2) = dat(jx+1,jy)
          help(3) = dat(jx,jy-1)
          help(4) = dat(jx-1,jy)

 200      continue

          ic = 0
          sum = 0.0
          do jj=1,4
            if(.not.loch(help(jj))) then
              sum = sum+help(jj)
              ic = ic+1
            endif
          enddo

          if (ic.gt.0) then
            dat(jx,jy) = sum/float(ic)  ! fill with mean of valid
            ipk = ipk+1                 ! neighbourpoints
          endif

     endif
     if (ipk .ge. icount) return    ! until all filled
enddo
enddo
go to 10

end subroutine fill

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
logical function loch(x)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real(kind=4), intent(in)        :: x
real(kind=4) :: edge

edge = -98.0
if (x.lt.edge) then
   loch = .true.
else
   loch = .false.
endif
return
end function loch



!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!+++++++++++++++++++++ step 3: end MODULE ++++++++++++++++++++++++++++++
end module read_nwp_data_module
