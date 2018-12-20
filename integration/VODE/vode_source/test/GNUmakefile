fsources = tvode.f

LINPACK_HOME = $(MICROPHYSICS_HOME)/util/LINPACK
BLAS_HOME = $(MICROPHYSICS_HOME)/util/BLAS

include ../GPackage.mak
include $(LINPACK_HOME)/GPackage.mak
include $(BLAS_HOME)/GPackage.mak

F90 = gfortran

tvode: $(fsources)
	$(F90) -o tvode $^ 

vpath %.f . ../ $(LINPACK_HOME) $(BLAS_HOME)


clean:
	$(RM) *.o tvode

