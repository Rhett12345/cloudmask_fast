module message_module

!C-----------------------------------------------------------------------
!C !F90                                                                  
!C
!C !Description: 
!C    This module is a message code for FY3/MERSI-II product code/
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
!C    National Satellite Meteorological Center, CMA 
!C  
!C !END
!C----------------------------------------------------------------------

!use names_module
!use platform_module

use constant

implicit none


CONTAINS
!+++++++++++++++++++ step 2: subroutines +++++++++++++++++++++++++++++++
!~~~~~~~~~~~~~~~~~~~ function 1: check file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
integer function checkfile(name)

!-----------------------------------------------------------------------
! !F90 checkfile
!
! !Description:
!    This program is to check the file status.
!
! !Input  parameters:
!    name            =  file path + file name
!          
! !Output parameters:
!    checkfile       =  id of file's information
!
!-----------------------------------------------------------------------

implicit none

! 1. define variables  
character(*), intent(in)  :: name
logical                   :: exist
character(len=11)         :: readability

! 2. begin program  
checkfile = 0

inquire( file = name,   &
         exist= exist,  &
         read = readability)
         
if (.not. exist) then
   checkfile = 1
endif
if (readability == 'NO') then
   checkfile = 2
endif

! 3. end function
end function checkfile
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

!~~~~~~~~~~~~~~~~~~~ subroutine 1: file information message output ~~~~~
subroutine file_message(message, status)

!-----------------------------------------------------------------------
! !F90 file_message
!
! !Description:
!    This program is to output the message of reading file.
!
! !Input  parameters:
!    message         =  message to point out the mistake file
!    severity        =  id of file have been checked
!          
! !Output parameters:
!    none
!
!-----------------------------------------------------------------------
   
implicit none

! 1. define variables
integer(kind=4), intent(in) :: status
character*(*), intent(in) :: message

! 2. begin program  
select case(status)
   case(0)  
!         do nothing for operational algorithm
   case(1)      
	      print*, "File is not exist: ", message
	      stop
   case(2)
		  print*, "File can not be opened: ", message
		  stop
end select
  
! 3. end subroutine     
end subroutine file_message
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

!~~~~~~~~~~~~~~~~~~~ subroutine 2: hdf message output ~~~~~~~~~~~~~~~~~~
subroutine hdf_info_message (message1, message2, status)

!-----------------------------------------------------------------------
! !F90 hdf_info_message
!
! !Description:
!    This program is to output the message of reading hdf data.
!
! !Input  parameters:
!    message1        =  message to point out the mistake hdf data
!    message2        =  message to point out the mistake dataset 
!    status          =  status of opening dataset
!          
! !Output parameters:
!    none
!
!-----------------------------------------------------------------------
   
implicit none

! 1. define variables
integer(kind=4), intent(in) :: status
character*(*), intent(in) :: message1, message2

! 2. begin program  
  if (status /= 0) then
    print*,'ERROR: Can not open ',trim(message1),' dataset:', message2
    stop
  endif

! 3. end subroutine     
end subroutine hdf_info_message
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

!~~~~~~~~~~~~~~~~~~~ subroutine 3: ncep binary data message output ~~~~~
subroutine bin_info_message (message)

!-----------------------------------------------------------------------
! !F90 bin_info_message
!
! !Description:
!    This program is to output the message of reading binary[ncep] data.
!
! !Input  parameters:
!    message        =  message to point out the mistake hdf data
!          
! !Output parameters:
!    none
!
!-----------------------------------------------------------------------
   
implicit none

! 1. define variables
character*(*), intent(in) :: message

! 2. begin program  
print*,'ERROR: Failed opening ',trim(message)
stop

! 3. end subroutine     
end subroutine bin_info_message
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
subroutine local_message_handler(message, severity, functionname)
   
   character*(*), intent(in) :: message, functionname
   integer,       intent(in) :: severity
   
   select case(severity)
      case(success)  
!        do nothing for operational algorithm
      case(warning)      
                  print*, "WARNING: ", message
      case(error)
                  print*, "ERROR: ", message
                  stop
      case(failure)
                  print*, "FAILURE, CRITICAL ERROR: ", message
                  stop
   end select
   
end subroutine local_message_handler
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!+++++++++++++++++++++ step 3: end module ++++++++++++++++++++++++++++++  
end module message_module
