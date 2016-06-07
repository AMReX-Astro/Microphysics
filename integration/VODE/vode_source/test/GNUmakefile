fsources = tvode.f

include ../GPackage.mak
include ../../BLAS/GPackage.mak

F90 = gfortran

tvode: $(fsources)
	$(F90) -o tvode $^

vpath %.f . ../ ../../BLAS


clean:
	$(RM) *.o tvode

