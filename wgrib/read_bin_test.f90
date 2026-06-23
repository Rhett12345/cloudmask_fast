program test 

REAL(KIND=4), DIMENSION(1:360,1:181,1:124) :: bin_data
INTEGER(KIND=4) :: status, io_err

OPEN(11,file=TRIM("/mnt/storage1/fy4test/JMA_Himawari08/nwp/fnl_20150125_00_00.grib2.bin"),status='old',access='direct',form="binary",&
recl=181*360*2*4,IOSTAT=io_err)

READ(11,rec=1) (((bin_data(i,j,k),i=1,360),j=1,181),k=1,2)
CLOSE(11)
print *, bin_data(:,:,1)
end
