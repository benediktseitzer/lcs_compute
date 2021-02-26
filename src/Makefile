##################################
#                                #
#      Makefile for LCS          #
#                                #
##################################

STUFF_NEEDED=stuff_needed
DEBUG = -check all -warn all

# ZIVHPC
F90 = ifort
#o = -O2 -openmp -fpp -convert big_endian -traceback -xHost -I${MKLROOT}/include/intel64/lp64 -mkl=sequential $(DEBUG)
#o = -O2 -fpp -convert big_endian -traceback -xHost -I${MKLROOT}/include/intel64/lp64 -mkl=sequential 
o = -O2 -openmp -fpp -convert big_endian -traceback -xHost -I${MKLROOT}/include/intel64/lp64 -mkl=sequential
omod = -O2 -fpp -convert big_endian -traceback -xHost -I${MKLROOT}/include/intel64/lp64 -mkl=sequential

LIBS =  ${MKLROOT}/lib/intel64/libmkl_blas95_lp64.a ${MKLROOT}/lib/intel64/libmkl_lapack95_lp64.a -lpthread -lm  -L../$(STUFF_NEEDED)/lib libjc.a libjpeg_bene.a
.SUFFIXES:
.SUFFIXES: .f90 .o

OBJECTS = rk4.o gridadvect.o calculateFTLE.o read_data.o write_output.o input_std.o variational.o
INCLUDE = parameter_module.o interp_module.o variational_module.o

all :     lcs

lcs   : lcs.o $(OBJECTS)
	$(F90) $(o) lcs.o $(OBJECTS) $(INCLUDE) -o lcs.exe $(LIBS)

lcs.o : lcs.f90 $(INCLUDE) Makefile
	$(F90) -c $(o) lcs.f90

rk4.o : rk4.f90 Makefile
	$(F90) -c $(o) rk4.f90
interp_module.o :interp_module.f90
	$(F90) -c $(omod)  interp_module.f90
parameter_module.o : parameter_module.f90
	$(F90) -c $(omod) parameter_module.f90
variational_module.o : variational_module.f90
	$(F90) -c $(omod) variational_module.f90
gridadvect.o : gridadvect.f90
	$(F90) -c $(o) gridadvect.f90
calculateFTLE.o :calculateFTLE.f90
	$(F90) -c $(o) calculateFTLE.f90
read_data.o :read_data.f90
	$(F90) -c $(o) read_data.f90 -I../$(STUFF_NEEDED)/include
write_output.o :write_output.f90
	$(F90) -c $(o) write_output.f90 -I../$(STUFF_NEEDED)/include
input_std.o :input_std.f90
	$(F90) -c $(o) input_std.f90  
variational.o : variational.f90
	$(F90) -c $(o) variational.f90

clean :
	rm -f *.o *.mod *_genmod* *~ *# *.exe
