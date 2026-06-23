module io_module_intermediate

! USE modules
use data_arrays_module
use names_module
use planck_module
use constant
use platform_module
use HDF5
use io_module                   !lyj
use get_ancil_data_module       !lyj

!---------------------------------------------
!use cloudmask_data_arrays       !lyj
use fylat_fy3mersi_cloud_mask   !lyj

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

CONTAINS
!+++++++++++++++++++ step 2: subroutines +++++++++++++++++++++++++++++++
!~~~~~~~~~~~~~~~~ Subroutine  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine fylat_write_out_intermediate()

!-----------------------------------------------------------------------
! !F90 write_out_arrays
!
! !Description:
!    The main program for write intermediate array out.
!
! !Input  parameters:
!    filename        = 
!
! !Output parameters:
!    none
!-----------------------------------------------------------------------

IMPLICIT NONE

! 1. define variables
integer (HID_T)                 :: file_id         ! File identifier
integer (HID_T)                 :: sds_id          ! Dataset id
integer (HID_T)                 :: dsp_id          ! Dataspace id
character (LEN=20)              :: dset_name       ! Dataset name
integer                         :: RANK            ! Dataset rank
integer (HSIZE_T), dimension(3) :: dims_3D         ! Dataset dimensions
integer (HSIZE_T), dimension(2) :: dims_2D         ! Dataset dimensions

INTEGER (HID_T)                 :: prp_id          ! Property list identifier
INTEGER (HSIZE_T), DIMENSION(3) :: chunck_dims_3D  ! chunk dimensions
INTEGER (HSIZE_T), DIMENSION(2) :: chunck_dims_2D  ! chunk dimensions

integer :: error                                   ! Error flag 

integer, dimension(:,:), pointer  :: H5TEMP
!integer(kind=4), allocatable, dimension(:,:)  :: H5TEMP

!-------------add by jincheng -----------
integer(kind=4) :: iline,ielem
real :: ptv11_12(1354,2030)
real :: ntv11_12(1354,2030)   

!allocate (ptv11_12(sat%nElem, sat%nLine),  &
!          ntv11_12(sat%nElem, sat%nLine))  
          
!----------------------------------------
do iline =1,sat%nLine
   do ielem=1,sat%nElem
      if ((sat%tbb_ir(ielem,iline,5) > 270.0) .and. (sat%bt_clr11(ielem,iline) > 270.0)) then
         ptv11_12(ielem,iline) = (sat%tbb_ir(ielem,iline,5) - sat%tbb_ir(ielem,iline,6)) -  &
                    (sat%bt_clr11(ielem,iline) - sat%bt_clr12(ielem,iline)) *(sat%tbb_ir(ielem,iline,5) - 260.0) / &
                    (sat%bt_clr11(ielem,iline) - 260.0)
      else
         ptv11_12(ielem,iline) = (sat%tbb_ir(ielem,iline,5) - sat%tbb_ir(ielem,iline,6))
      endif
   enddo
enddo

ntv11_12 = (sat%bt_clr11 - sat%bt_clr12)-(sat%tbb_ir(:,:,5) - sat%tbb_ir(:,:,6))    

!----------------------------------------


!======================================================================
! 2. begin program
print*,'    ... fylat write out fy3/MERSI_II Intermediate HDF5 file !!! '
 
call h5open_f(error)

call h5fcreate_f(trim(fy3_intermediate), H5F_ACC_TRUNC_F, file_id, error)

call h5pcreate_f(H5P_DATASET_CREATE_F, prp_id, error) 
call h5pset_deflate_f(prp_id, 5, error)

!====================================
! --- Write 3D dataset 
RANK = 3

!------------------------------------
! --- Write sat%ref_vis(sat%nElem, sat%nLine, sat%nvis) in 3D dataset
print*, "--write ref_vis in 3D dataset"
chunck_dims_3D(1) = 100
chunck_dims_3D(2) = 100
chunck_dims_3D(3) = 10
call h5pset_chunk_f(prp_id, RANK, chunck_dims_3D, error)

dims_3D(1) = sat%nElem
dims_3D(2) = sat%nLine
dims_3D(3) = sat%nvis
call h5screate_simple_f(RANK, dims_3D, dsp_id, error)

dset_name = "ref_vis"
call h5dcreate_f(file_id, dset_name, H5T_NATIVE_REAL, dsp_id, sds_id, error, &
                 prp_id)
call h5dwrite_f (sds_id, H5T_NATIVE_REAL, sat%ref_vis, dims_3D, error)
call h5dclose_f (sds_id, error)      
call h5sclose_f (dsp_id, error)

!------------------------------------
! --- Write sat%tbb_ir(sat%nElem, sat%nLine,  sat%nir ) in 3D dataset
print*, "--write tbb_ir in 3D dataset"
chunck_dims_3D(1) = 100
chunck_dims_3D(2) = 100
chunck_dims_3D(3) = 3
call h5pset_chunk_f(prp_id, RANK, chunck_dims_3D, error)

dims_3D(1) = sat%nElem
dims_3D(2) = sat%nLine
dims_3D(3) = sat%nir
call h5screate_simple_f(RANK, dims_3D, dsp_id, error)

dset_name = "tbb_ir"
call h5dcreate_f(file_id, dset_name, H5T_NATIVE_REAL, dsp_id, sds_id, error, &
                 prp_id)
call h5dwrite_f (sds_id, H5T_NATIVE_REAL, sat%tbb_ir, dims_3D, error)
call h5dclose_f (sds_id, error)      
call h5sclose_f (dsp_id, error)

    
!====================================
! --- Write 2D dataset 
RANK = 2

chunck_dims_2D(1) = 100
chunck_dims_2D(2) = 100
call h5pset_chunk_f(prp_id, RANK, chunck_dims_2D, error)

dims_2D(1) = sat%nElem
dims_2D(2) = sat%nLine
call h5screate_simple_f(RANK, dims_2D, dsp_id, error)

!------------------------------------
! --- Write sat%snow_mask(sat%nElem,sat%nLine)
print*, "--write snow_mask"
allocate ( H5TEMP(sat%nElem,sat%nLine), stat=error )
H5TEMP = -99
H5TEMP = sat%snow_mask

dset_name = "snow_mask"
call h5dcreate_f(file_id, dset_name, H5T_STD_I8LE, dsp_id, sds_id, error, &
                 prp_id)
call h5dwrite_f (sds_id, H5T_NATIVE_INTEGER, H5TEMP, dims_2D, error)
call h5dclose_f (sds_id, error)      

!deallocate (H5TEMP, stat=error)

!------------------------------------
! --- Write sat%eco(sat%nElem,sat%nLine)
print*, "--write eco"
!allocate ( H5TEMP(sat%nElem,sat%nLine), stat=error )
H5TEMP = -99
H5TEMP = sat%eco

dset_name = "eco"
call h5dcreate_f(file_id, dset_name, H5T_STD_I8LE, dsp_id, sds_id, error, &
                 prp_id)
call h5dwrite_f (sds_id, H5T_NATIVE_INTEGER, H5TEMP, dims_2D, error)
call h5dclose_f (sds_id, error)      

deallocate (H5TEMP, stat=error)

!------------------------------------
! --- Write out_pwater(sat%nElem,sat%nLine)
print*, "--write percip_water"
dset_name = "precip_water"
call h5dcreate_f(file_id, dset_name, H5T_NATIVE_REAL, dsp_id, sds_id, error, &
                 prp_id)
call h5dwrite_f (sds_id, H5T_NATIVE_REAL, out_pwater, dims_2D, error)
call h5dclose_f (sds_id, error)      

deallocate ( out_pwater )

!------------------------------------
! --- Write out_sfctmp(sat%nElem,sat%nLine)
print*, "--write sfctmp"
dset_name = "sfctmp"
call h5dcreate_f(file_id, dset_name, H5T_NATIVE_REAL, dsp_id, sds_id, error, &
                 prp_id)
call h5dwrite_f (sds_id, H5T_NATIVE_REAL, out_sfctmp, dims_2D, error)
call h5dclose_f (sds_id, error)      

deallocate ( out_sfctmp )

!------------------------------------
print*, "--write polar   "
dset_name = "polar"
call h5dcreate_f(file_id, dset_name, H5T_STD_I8LE, dsp_id, sds_id, error, &
                 prp_id)
call h5dwrite_f (sds_id, H5T_NATIVE_INTEGER, out_polar, dims_2D, error)
call h5dclose_f (sds_id, error)

deallocate (  out_polar   )

!------------------------------------
print*, "--write day     "
dset_name = "day"
call h5dcreate_f(file_id, dset_name, H5T_STD_I8LE, dsp_id, sds_id, error, &
                 prp_id)
call h5dwrite_f (sds_id, H5T_NATIVE_INTEGER, out_day, dims_2D, error)
call h5dclose_f (sds_id, error)

deallocate (  out_day     )

!------------------------------------
print*, "--write night   "
dset_name = "night"
call h5dcreate_f(file_id, dset_name, H5T_STD_I8LE, dsp_id, sds_id, error, &
                 prp_id)
call h5dwrite_f (sds_id, H5T_NATIVE_INTEGER, out_night, dims_2D, error)
call h5dclose_f (sds_id, error)

deallocate (  out_night   )

!------------------------------------
print*, "--write land    "
dset_name = "land"
call h5dcreate_f(file_id, dset_name, H5T_STD_I8LE, dsp_id, sds_id, error, &
                 prp_id)
call h5dwrite_f (sds_id, H5T_NATIVE_INTEGER, out_land, dims_2D, error)
call h5dclose_f (sds_id, error)

deallocate (  out_land    )

!------------------------------------
print*, "--write water   "
dset_name = "water"
call h5dcreate_f(file_id, dset_name, H5T_STD_I8LE, dsp_id, sds_id, error, &
                 prp_id)
call h5dwrite_f (sds_id, H5T_NATIVE_INTEGER, out_water, dims_2D, error)
call h5dclose_f (sds_id, error)

deallocate (  out_water   )

!------------------------------------
print*, "--write coast   "
dset_name = "coast"
call h5dcreate_f(file_id, dset_name, H5T_STD_I8LE, dsp_id, sds_id, error, &
                 prp_id)
call h5dwrite_f (sds_id, H5T_NATIVE_INTEGER, out_coast, dims_2D, error)
call h5dclose_f (sds_id, error)

deallocate (  out_coast   )

!------------------------------------
print*, "--write snglnt  "
dset_name = "snglnt"
call h5dcreate_f(file_id, dset_name, H5T_STD_I8LE, dsp_id, sds_id, error, &
                 prp_id)
call h5dwrite_f (sds_id, H5T_NATIVE_INTEGER, out_snglnt, dims_2D, error)
call h5dclose_f (sds_id, error)

deallocate (  out_snglnt  )

!------------------------------------
print*, "--write snow    "
dset_name = "snow"
call h5dcreate_f(file_id, dset_name, H5T_STD_I8LE, dsp_id, sds_id, error, &
                 prp_id)
call h5dwrite_f (sds_id, H5T_NATIVE_INTEGER, out_snow, dims_2D, error)
call h5dclose_f (sds_id, error)

deallocate (  out_snow    )

!------------------------------------
print*, "--write ice     "
dset_name = "ice"
call h5dcreate_f(file_id, dset_name, H5T_STD_I8LE, dsp_id, sds_id, error, &
                 prp_id)
call h5dwrite_f (sds_id, H5T_NATIVE_INTEGER, out_ice, dims_2D, error)
call h5dclose_f (sds_id, error)

deallocate (  out_ice     )

!------------------------------------
print*, "--write desert  "
dset_name = "desert"
call h5dcreate_f(file_id, dset_name, H5T_STD_I8LE, dsp_id, sds_id, error, &
                 prp_id)
call h5dwrite_f (sds_id, H5T_NATIVE_INTEGER, out_desert, dims_2D, error)
call h5dclose_f (sds_id, error)

deallocate (  out_desert  )

!------------------------------------
print*, "--write uniform "
dset_name = "uniform"
call h5dcreate_f(file_id, dset_name, H5T_STD_I8LE, dsp_id, sds_id, error, &
                 prp_id)
call h5dwrite_f (sds_id, H5T_NATIVE_INTEGER, out_uniform, dims_2D, error)
call h5dclose_f (sds_id, error)

deallocate (  out_uniform )

!------------------------------------
print*, "--write shadow  "
dset_name = "shadow"
call h5dcreate_f(file_id, dset_name, H5T_STD_I8LE, dsp_id, sds_id, error, &
                 prp_id)
call h5dwrite_f (sds_id, H5T_NATIVE_INTEGER, out_shadow, dims_2D, error)
call h5dclose_f (sds_id, error)

deallocate (  out_shadow  )

!-------------------------------------     add by jincheng   17.5.13
print*,"--write clr11"
dset_name = "bt_clr11"
call h5dcreate_f(file_id, dset_name, H5T_NATIVE_REAL, dsp_id, sds_id, error, &
                 prp_id)
call h5dwrite_f (sds_id, H5T_NATIVE_REAL, sat%bt_clr11, dims_2D, error)
call h5dclose_f (sds_id, error)      

!deallocate ( sat%bt_clr11 )

!-------------------------------------     add by jincheng   17.5.13
print*,"--write clr12"
dset_name = "bt_clr12"
call h5dcreate_f(file_id, dset_name, H5T_NATIVE_REAL, dsp_id, sds_id, error, &
                 prp_id)
call h5dwrite_f (sds_id, H5T_NATIVE_REAL, sat%bt_clr12, dims_2D, error)
call h5dclose_f (sds_id, error)   

!-------------------------------------     add by jincheng   17.5.13
print*,"--write pfmft"
dset_name = "pfmft"
call h5dcreate_f(file_id, dset_name, H5T_NATIVE_REAL, dsp_id, sds_id, error, &
                 prp_id)
call h5dwrite_f (sds_id, H5T_NATIVE_REAL, ptv11_12, dims_2D, error)
call h5dclose_f (sds_id, error) 

!-------------------------------------     add by jincheng   17.5.13
print*,"--write nfmft"
dset_name = "nfmft"
call h5dcreate_f(file_id, dset_name, H5T_NATIVE_REAL, dsp_id, sds_id, error, &
                 prp_id)
call h5dwrite_f (sds_id, H5T_NATIVE_REAL, ntv11_12, dims_2D, error)
call h5dclose_f (sds_id, error) 




!====================================
!close space for writing 2D dataset 
call h5sclose_f (dsp_id, error)
  
!====================================================
! Close property and file
call h5pclose_f(prp_id, error)

call h5fclose_f(file_id, error)

call h5close_f(error)

!deallocate (ptv11_12,  &                   !jincheng
!            ntv11_12)
			
end subroutine fylat_write_out_intermediate
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!-------------------------- END MODULE ---------------------------------
end module io_module_intermediate


