module fylat_fy3mersi_cloud_amount

!C-----------------------------------------------------------------------
!C !F90                                                                  
!C
!C !Description: 
!C    cloud mask algorithm for fy3/mersi sensor
!C
!C !Input parameters
!C    none
!C 
!C !Output parameters
!C    none
!C
!C  
!C !END
!C----------------------------------------------------------------------

use names_module
use data_arrays_module
use constant


implicit none

!+++++++++++++++++++ step 1: Variables +++++++++++++++++++++++++++++++++


contains
!+++++++++++++++++++ step 2: Subroutines +++++++++++++++++++++++++++++++
!~~~~~~~~~~~~~~~~~~~ subroutine 1: fy3mersi_cloud_amount ~~~~~~~~~~~~~~~~~
subroutine fy3mersi_cloud_amount()

!-----------------------------------------------------------------------
!                          PROCESSING SECTION
!
!        The Cloud Mask is processed pixel by pixel using a sliding
!        box approach. Regional context data is stored for use
!        in uniformity tests for uncertain pixels.  The edge of 
!        the area in which you are processing, (outline of region)
!        will be processed but will not include uniformity tests.
!----------------------------------------------------------------------



integer(kind=4) :: iline,ielem, i, j, k, ix_nwp, iy_nwp
integer(kind=4) :: cc 
integer(kind=1) :: qflag
!real(kind=4)    :: u_wind
!real(kind=4)    :: v_wind
!real(kind=4)    :: tpw
integer(kind=4) :: mpos_x, mpos_y
integer(kind=1), dimension(5,5)  :: cm_5x5, qa_5x5
integer(kind=1), dimension(:,:), pointer         :: ca
integer(kind=1), dimension(:,:), pointer         :: qa
integer(kind=4), dimension(:,:), pointer         :: lat5km
integer(kind=4), dimension(:,:), pointer         :: lon5km

print*,'    ... fylat retrieve fy3/MERSI_II Cloud Amount start !!! '

!--- set local POINTERs to output structures
ca => cloud_amount
qa => cloud_amount_qa
lat5km => lat_5km
lon5km => lon_5km


!======================================================================
! Loop over pixels in this segment
!======================================================================
!b0=1
!b2b1
!00=云              0 
!01=可能晴空         1
!10=可信的晴空       2
!11=可信度高的晴空    3

! ++ calculate cloud amount
line_loop_1: do j= 1, iy_5km

element_loop_1: do i= 1, ix_5km

     !print*,'i,j',i,j
     cm_5x5 = cm_tmp(1+5*(i-1):5+5*(i-1),1+5*(j-1):5+5*(j-1),1)
     qa_5x5 = cm_tmp(1+5*(i-1):5+5*(i-1),1+5*(j-1):5+5*(j-1),2)
     
     call calculate_cloud_cover_5x5(cm_5x5, qa_5x5, cc, qflag)
     
     ca(i,j) = cc
     qa(i,j) = qflag
     
     mpos_x = 3+5*(i-1)
     mpos_y = 3+5*(j-1)
     lon5km(i,j) = int(100.0*geo%lon(mpos_x,mpos_y))
     lat5km(i,j) = int(100.0*geo%lat(mpos_x,mpos_y))
     
     cm_5x5 = 0
     qa_5x5 = 0
     
end do element_loop_1
    
end do line_loop_1


ca => null()
qa => null()
lat5km => null()
lon5km => null()  

end subroutine fy3mersi_cloud_amount
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine calculate_cloud_cover_5x5(cm_5x5_matrix,   &
                                     qa_5x5_matrix,   &
                                     cc,qflag)

integer(kind=1), intent(in), dimension(5,5) :: cm_5x5_matrix, qa_5x5_matrix
integer(kind=4), intent(out)  :: cc 
integer(kind=1), intent(out)  :: qflag
real(kind=4) :: num_valid_pixel, num_cloudy_pixel
integer :: ii,jj

num_valid_pixel  = 0.0
num_cloudy_pixel = 0.0 
qflag = 0  ! bad
cc = -999

do ii = 1, 5
do jj = 1, 5

   if (qa_5x5_matrix(ii,jj) == 1) then
       num_valid_pixel = num_valid_pixel + 1.0
       
       if (cm_5x5_matrix(ii,jj) < 2) then  ! cloudy + pro cloudy
           num_cloudy_pixel = num_cloudy_pixel + 1.0
       endif
       
   endif
   
enddo
enddo

if (num_valid_pixel > 15 .and. num_valid_pixel < 25) then ! low quality
    qflag = 1
    cc = int((num_cloudy_pixel/num_valid_pixel)*100.0)
endif 

if (num_valid_pixel == 25) then   ! high quality
    qflag = 2
    cc = int((num_cloudy_pixel/num_valid_pixel)*100.0)
endif

end subroutine calculate_cloud_cover_5x5
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!-------------------------- END MODULE ---------------------------------
end module fylat_fy3mersi_cloud_amount
