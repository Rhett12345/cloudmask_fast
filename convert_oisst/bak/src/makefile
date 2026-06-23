# Makefile
#
FC       = ifort -c -O3
F90      = ifort
hdf5     = h5fc

INC      = -I /wind1/software/mathlib/hdf5/1.8.12/intel/include
INC2	 = -I /wind1/home/minm/software/netcdf4/include

#LIB1
#LIB     = -I/home/huxq/opt/intel_lib_new/include -L/home/huxq/opt/intel_lib_new/lib /home/huxq/opt/intel_lib_new/lib/libhdf5hl_fortran.a /home/huxq/opt/intel_lib_new/lib/libhdf5_hl.a /home/huxq/opt/intel_lib_new/lib/libhdf5_fortran.a /home/huxq/opt/intel_lib_new/lib/libhdf5.a -L/home/huxq/opt/intel_lib_new/lib -L/home/huxq/opt/intel_lib_new/lib -lsz -lz -lrt -ldl -lm -Wl,-rpath -Wl,/home/huxq/opt/intel_lib_new/lib
#LIB = -I/state/partition1/apps/LIB/include -L/state/partition1/apps/LIB/lib /state/partition1/apps/LIB/lib/libhdf5hl_fortran.a /state/partition1/apps/LIB/lib/libhdf5_hl.a /state/partition1/apps/LIB/lib/libhdf5_fortran.a /state/partition1/apps/LIB/lib/libhdf5.a -L/state/partition1/apps/LIB/lib -L/state/partition1/apps/LIB/lib -lsz -lz -lrt -lm -Wl,-rpath -Wl,/state/partition1/apps/LIB/lib
LIB = -I/wind1/software/mathlib/hdf5/1.8.12/intel/include -L/wind1/software/mathlib/hdf5/1.8.12/intel/lib /wind1/software/mathlib/hdf5/1.8.12/intel/lib/libhdf5hl_fortran.a /wind1/software/mathlib/hdf5/1.8.12/intel/lib/libhdf5_hl.a /wind1/software/mathlib/hdf5/1.8.12/intel/lib/libhdf5_fortran.a /wind1/software/mathlib/hdf5/1.8.12/intel/lib/libhdf5.a -lz -lrt -ldl -lm -Wl,-rpath -Wl,/wind1/software/mathlib/hdf5/1.8.12/intel/lib
LIB2 = -L /wind1/home/minm/software/netcdf4/lib -lnetcdf -lnetcdff -lcurl


F90_SRCS = io_module.f90 oisst_daily_nc2hdf5.f90


OBJS     =  $(F90_SRCS:.f90=.o)

EXE      =  ../oisst_daily_nc2hdf5


# general rules


# Target file (executable)
$(EXE): $(OBJS)
	$(F90) $(OBJS) $(LIB) -o $@ 
	
#-------------------------	
# make clean
clean:
	-rm -f $(OBJS) *.mod
	-rm -f $(OBJS) *.o
	-rm -f ../oisst_daily_nc2hdf5

io_module.o: io_module.f90 
	$(F90) -c -O3 io_module.f90 $(INC)  #$(INC2)
	
oisst_daily_nc2hdf5.o: oisst_daily_nc2hdf5.f90
	$(FC) oisst_daily_nc2hdf5.f90 $(INC) #$(INC2)
