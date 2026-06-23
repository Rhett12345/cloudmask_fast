module io_module

!-----------------------------------------------------------------------
! f90 io_module
!
!Description: 
!   io module for nc data2 reading and writing to HDF5.
!
! Author: Min Min
! E-mail: minmin@cma.gov.cn
!
!-----------------------------------------------------------------------


use HDF5 
!use netcdf

implicit none

real(kind=4),allocatable,dimension(:,:,:) :: sst0
real(kind=4),allocatable,dimension(:,:)   :: sst


contains
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine read_nc_file(fname,dname,nd)
!use netcdf
character(len=300) :: path,fname
character(len=3)   :: dname,cname
integer :: nd

!===== 1.1.hdf4 Function declaration
INTEGER sd_id, sds_id, attr_id, status,i,j
INTEGER, DIMENSION(2) :: start1, stride1
INTEGER (HID_T)       :: sd_id5, sds_id5, attr_id5
INTEGER (HSIZE_T), DIMENSION(3) ::  edges25

allocate(sst0(1440,720,366) , &
         sst(1440,720))

  start1   = 0
  stride1  = 1
  edges25(1)= 1440
  edges25(2)= 720
  edges25(3)= 366
  
  CALL h5open_f(status)
  CALL h5fopen_f (trim(fname), H5F_ACC_RDONLY_F, sd_id5, status)
  
  CALL h5dopen_f(sd_id5, trim(dname), sds_id5, status)
  CALL H5DREAD_F(sds_id5, H5T_NATIVE_REAL, sst0, edges25, status)
  CALL h5dclose_f(sds_id5, status)
  
  ! Terminate access to the SD interface and close the file. 
  CALL h5fclose_f(sd_id5, status)  ! close 
  CALL h5close_f(status)
  
  !print*,sst0(736,507,nd)
  !print*,sst0(736,508,nd)
  !print*,sst0(736,509,nd)
  !print*,sst0(736,510,nd)
  !print*,sst0(736,511,nd)
  
  do i = 1, 720
     j = 721-i
     sst(1:720,j) = sst0(721:1440,i,nd)
     sst(721:1440,j) = sst0(1:720,i,nd)
  enddo
  deallocate (sst0)

!integer:: ierr,ncid,varid,dimid,len_file,xtype,ndims,dimids(5),natts,len

!include 'netcdf.inc'

!allocate(sst(1440,720,366))
!print*,' ... read nc file name', trim(fname)
!ierr = nf90_open(trim(fname),nf90_nowrite,ncid)
!ierr=nf90_inq_dimid (ncid,trim(dname),dimid) 
!ierr=nf90_Inquire_Dimension(ncid,dimid,cname,len)	
!print*, trim(nf90_strerror(ierr))	
!print*,'r',ierr,'ff ',ncid,' 22',trim(dname)
!ierr=nf90_inq_varid(ncid,trim(dname), varid) 

!print*,'33'
!ierr=nf90_Inquire_Variable(ncid,varid,vname,xtype,ndims,dimids,natts)   
!ierr=nf90_get_var(ncid,varid,sst)
!print*,'44'
!ierr=nf90_close(ncid)   


end subroutine read_nc_file
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine write_hdf5_file(fname,dname)
!use netcdf
character(len=300) :: path,fname
character(len=3)   :: dname,cname

!===== 1.1.hdf4 Function declaration
integer :: error
INTEGER (HID_T) :: file_id1               
INTEGER (HID_T) :: sds_id            ! sds id for longitude   [unit:degree]

INTEGER (HID_T) :: dsp_id              ! dsp id for longitude   [unit:degree]

INTEGER (HID_T) :: attr_id           ! sds id for longitude   [unit:degree]
INTEGER (HID_T) :: asp_id           ! sds id for longitude   [unit:degree]
INTEGER (HID_T) :: type_id           ! sds id for longitude   [unit:degree]

!xxxxxxxxxxxxx other]
INTEGER (HSIZE_T), DIMENSION(2)         :: dims_sp

INTEGER, PARAMETER :: RANK1 = 1 ! Dataset rank
INTEGER, PARAMETER :: RANK2 = 2 ! Dataset rank
INTEGER, PARAMETER :: RANK3 = 3 ! Dataset rank

  CALL h5open_f(error)
  CALL h5fcreate_f(trim(fname), H5F_ACC_TRUNC_F, file_id1, error)
  
  dims_sp = (/1440,720/)
  CALL h5screate_simple_f(RANK2, dims_sp, dsp_id, error)
  CALL h5dcreate_f(file_id1, "sst", H5T_NATIVE_REAL, dsp_id, sds_id, error)
  CALL h5dwrite_f (sds_id, H5T_NATIVE_REAL, sst, dims_sp, error)
  CALL h5sclose_f (dsp_id, error)
  CALL h5dclose_f (sds_id, error)
  
    ! Terminate access to the SD interface and close the file. 
  CALL h5fclose_f(file_id1, error)

  CALL h5close_f(error)
   
  deallocate(sst)


end subroutine write_hdf5_file
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine julian (IY, IM, ID, IH, MIT, JD)

!-----------------------------------------------------------------------
! !F90 julian
!
! !Description:
!    The main program computes julian day (1-365/366)
!
! !Input  parameters:
!    IY               = INTEGER year      
!    IM               = INTEGER month
!    ID               = INTEGER day
!    IH               = INTEGER hour
!    MIT              = INTEGER minute
!
!
! !Output parameters:
!    JD               = julian day
!
!-----------------------------------------------------------------------

! 1. define variables
INTEGER(KIND=4), INTENT(in) :: IY, IM, ID, IH, MIT
REAL(KIND=8), INTENT(out)   :: JD ! julian day
REAL(KIND=4)    :: XI, XJ
INTEGER(KIND=4) :: IY1, IM1

! 2. begin program
IF (IM <= 2) THEN   ! january & february
  IY1 = int(IY-1)
  IM1 = int(IM+12)
  JD = dble(int( 365.25*(IY1 + 4716.0)) + int( 30.6001*( IM1 + 1.0)) + 2.0 - &
       int( IY1/100.0 ) + int( int( IY1/100.0 )/4.0 ) + ID - 1524.5) + &
       dble((IH + MIT/60.+0./3600.)/24.)
     
ELSE

  JD = dble( int( 365.25*(IY + 4716.0)) + int( 30.6001*( IM + 1.0)) + 2.0 - &
       int( IY/100.0 ) + int( int( IY/100.0 )/4.0 ) + ID - 1524.5) + &
       dble((IH + MIT/60.+0./3600.)/24.)
     
ENDIF

end subroutine julian
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE  julian_to_date (JD, year, month, day, hour, mint)

!-----------------------------------------------------------------------
! !F90 julian converter
!
! This function converts the Julian dates to Gregorian dates.
!
! Syntax:
! [day,month,year,hour,mint] = julian_to_date(JD)
!
! !Input  parameters:
!    JD               = julian day
!
! !Output parameters:
!    year             = INTEGER year      
!    month            = INTEGER month
!    day              = INTEGER day
!    hour             = INTEGER hour
!    mint             = INTEGER minute
!
!-----------------------------------------------------------------------

! 1. define variables
INTEGER(KIND=4), INTENT(out) :: year, month, day, hour, mint
REAL(KIND=8), INTENT(in)    :: JD ! julian day
REAL(KIND=4)     :: I, D, E, G   
REAL(KIND=4)    :: B, C, Fr
INTEGER(KIND=4)  :: A, a4

I = int(JD + 0.5)
Fr = abs( I - ( JD + 0.5) )

IF (I >= 2299160. ) THEN
     A = int( ( I- 1867216.25 ) / 36524.25 )
     a4 = int( A / 4 )
     B = I + 1. + float(A - a4)
ELSE
     B = I
ENDIF

C = B + 1524.
D = int( ( C - 122.1 ) / 365.25 )
E = int( 365.25 * D )
G = int( ( C - E ) / 30.6001 )
day = int( C - E + Fr - int( 30.6001 * G ) )

IF (G <= 13.5 ) THEN
    month = int(G - 1)
ELSE
    month = int(G - 13)
ENDIF

IF (month > 2.5) THEN
    year = int(D - 4716)
ELSE
    year = int(D - 4715)
ENDIF

hour = int( Fr * 24. )
mint = int( abs( hour -( Fr * 24. ) ) * 60. )

! 3. END
END SUBROUTINE julian_to_date
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~ SUBROUTINE 2: convert integers to characters and vice versa ~~~
SUBROUTINE ICNVRT(WAY,NUM,STRING,LENGTH,IERR)

!-----------------------------------------------------------------------
! !F90 ICNVRT
!
! !Description:
!        This SUBROUTINE does an INTEGER-to-CHARACTER conversion
!        or a characater-to-INTEGER conversion depENDing on the
!        INTEGER WAY:
!                IF WAY = 0 THEN an INTEGER-to-CHARACTER conversion
!                is done. IF WAY .NE. 0 THEN a CHARACTER-to-INTEGER
!                conversion is done.
!
! !USAGE:
!U
!U        CALL ICNVRT(WAY,NUM,STRING)
!U             where WAY, NUM, STRING, and LENGTH are defined below.
!U
!U        Example: CALL ICNVRT(0,1000,STRING,LENGTH)
!U                 on RETURN STRING = '1000' and
!U                 LENGTH = 4.
!         
! !Input  parameters:
!    WAY - INTEGER::; Determines which way the conversion goes:
!              IF WAY = 0 THEN an INTEGER-to-CHARACTER conversion
!                         is performed;
!              IF WAY.NE.0 THEN a CHARACTER-to-INTEGER conversion
!                         is performed.
!
!    NUM - INTEGER::; an input only IF WAY = 0. NUM is the INTEGER
!                number to be converted to a CHARACTER expression.
!
!    STRING - CHARACTER; an input only IF WAY .NE. 0. STRING
!                is the CHARACTER expression to be converted to an
!                INTEGER value. It contain no decimal points or 
!                non-numeric characters other than possibly a
!                sign. IF STRING contains  a '+' sign, it will be
!                stripped of it on RETURN.
!
! !Output parameters:
!    NUM - INTEGER::; contains the INTEGER:: representation of 
!                STRING.
!
!    STRING - CHARACTER; contains the CHARACTER representation of NUM.
!
!    LENGTH - INTEGER::; The length of STRING to the first blank.
!                  The signIFicant part of STRING can be accessed with
!                  the declaration STRING(1:LENGTH).
!
!    IERR - INTEGER:: variable giving RETURN condition:
!                IERR = 0 for normal RETURN;
!                IERR = 1 IF NUM cannot be converted to STRING because
!                       STRING is too short or STRING cannot be
!                       converted to NUM because STRING is too long.
!                IERR = 2 IF STRING contained a non-numeric CHARACTER
!                       other than a leading sign or something went
!                       wrong with an INTEGER-to-CHARACTER conversion.
!
! !Other
!       ALGORITHM:
!A
!A         Nothing noteworthy, except that this SUBROUTINE will work
!A          for strange CHARACTER sets where the CHARACTER '1' doesn't
!A          follow '0', etc.
!A
!       MACHINE DEPENDENCIES: CM
!M          The parameter MAXINT (below) should be set to the
!M          number of digits that an INTEGER:: data type can have
!M          not including leading signs. For VAX FORTRAN V4.4-177
!M          MAXINT = 10.
!M
!M          NOTE: Under VAX FORTRAN V4.4-177, the
!M          error condition IERR = 1 will never occur for an
!M          INTEGER-to-CHARACTER conversion IF STRING
!M          is ALLOCATEd at least 11 bytes (CHARACTER*11).
!M
!       HISTORY:
!H
!H      written by:             bobby bodenheimer
!H      date:                   september 1986
!H      current version:        1.0
!H      modIFications:          none
!H
!       ROUTINES CALLED:
!C
!C          NONE.
!C
!----------------------------------------------------------------------
!       written for:    The CASCADE Project
!                       Oak Ridge National Laboratory
!                       U.S. Department of Energy
!                       contract number DE-AC05-840R21400
!                       subcontract number 37B-7685 S13
!                       organization:  The University of Tennessee
!----------------------------------------------------------------------
!       THIS SOFTWARE IS IN THE PUBLIC DOMAIN
!       NO RESTRICTIONS ON ITS USE ARE IMPLIED
!----------------------------------------------------------------------


! Global Variables.
!
 INTEGER(KIND=4), INTENT(in):: WAY
 INTEGER(KIND=4), INTENT(out)::  LENGTH, IERR
 INTEGER(KIND=4), INTENT(inout):: NUM
 CHARACTER(LEN=*), INTENT(inout):: STRING
!
!
! Local Variables
!
      INTEGER(KIND=4)::       I
      INTEGER(KIND=4)::       MNUM
      INTEGER(KIND=4)::       M
      logical::       NEG
!
      INTEGER, parameter::MAXINT=10
!
      NEG = .FALSE.
      IERR = 0
!
!  INTEGER-to-CHARACTER conversion.
!
      IF (WAY == 0) THEN
         STRING = " "
         IF (NUM < 0) THEN
            NEG = .TRUE.
            MNUM = -NUM
            LENGTH = INT(LOG10(REAL(MNUM))) + 1
         ELSE IF (NUM == 0) THEN
            MNUM = NUM
            LENGTH = 1
         ELSE
            MNUM = NUM
            LENGTH = INT(LOG10(REAL(MNUM))) + 1
         END IF
         IF (LENGTH > LEN(STRING)) THEN
            IERR = 1
            RETURN
         END IF
ten:     DO I=LENGTH,1,-1    
            M=INT(REAL(MNUM)/10**(I-1))
            IF (M == 0) THEN
               STRING(LENGTH-I+1:LENGTH-I+1) = "0"
            ELSE IF (M == 1) THEN
               STRING(LENGTH-I+1:LENGTH-I+1) = "1"
            ELSE IF (M == 2) THEN
               STRING(LENGTH-I+1:LENGTH-I+1) = "2"
            ELSE IF (M == 3) THEN
               STRING(LENGTH-I+1:LENGTH-I+1) = "3"
            ELSE IF (M == 4) THEN
               STRING(LENGTH-I+1:LENGTH-I+1) = "4"
            ELSE IF (M == 5) THEN
               STRING(LENGTH-I+1:LENGTH-I+1) = "5"
            ELSE IF (M == 6) THEN
               STRING(LENGTH-I+1:LENGTH-I+1) = "6"
            ELSE IF (M == 7) THEN
               STRING(LENGTH-I+1:LENGTH-I+1) = "7"
            ELSE IF (M == 8) THEN
               STRING(LENGTH-I+1:LENGTH-I+1) = "8"
            ELSE IF (M == 9) THEN
               STRING(LENGTH-I+1:LENGTH-I+1) = "9"
            ELSE
               IERR = 2
               RETURN
            END IF
            MNUM = MNUM - M*10**(I-1)
         END DO ten

         IF (NEG .eqv. .true.) THEN
            STRING = "-"//STRING
            LENGTH = LENGTH + 1
         END IF
!
!  CHARACTER-to-INTEGER conversion.
!
      ELSE
         IF (STRING(1:1) == "-") THEN
            NEG = .TRUE.
            STRING = STRING(2:LEN(STRING))
         END IF
         IF (STRING(1:1) == "+") STRING = STRING(2:LEN(STRING))
         NUM = 0
         LENGTH = INDEX(STRING," ") - 1
         IF (LENGTH > MAXINT) THEN
            IERR = 1
            RETURN
         END IF
twenty:  DO I=LENGTH,1,-1
            IF (STRING(LENGTH-I+1:LENGTH-I+1) == "0") THEN
               M = 0
            ELSE IF (STRING(LENGTH-I+1:LENGTH-I+1) == "1") THEN
               M = 1
            ELSE IF (STRING(LENGTH-I+1:LENGTH-I+1) == "2") THEN
               M = 2
            ELSE IF (STRING(LENGTH-I+1:LENGTH-I+1) == "3") THEN
               M = 3
            ELSE IF (STRING(LENGTH-I+1:LENGTH-I+1) == "4") THEN
               M = 4
            ELSE IF (STRING(LENGTH-I+1:LENGTH-I+1) == "5") THEN
               M = 5
            ELSE IF (STRING(LENGTH-I+1:LENGTH-I+1) == "6") THEN
               M = 6
            ELSE IF (STRING(LENGTH-I+1:LENGTH-I+1) == "7") THEN
               M = 7
            ELSE IF (STRING(LENGTH-I+1:LENGTH-I+1) == "8") THEN
               M = 8
            ELSE IF (STRING(LENGTH-I+1:LENGTH-I+1) == "9") THEN
               M = 9
            ELSE
               IERR = 2
               RETURN
            END IF
            NUM = NUM + INT(10**(I-1))*M
         END DO twenty

         IF (NEG .eqv. .true.) THEN
            NUM = -NUM
            STRING = '-'//STRING
            LENGTH = LENGTH + 1
         END IF
      END IF
!
!  Last lines of ICNVRT
!
   RETURN
   
END SUBROUTINE ICNVRT
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
end module io_module
