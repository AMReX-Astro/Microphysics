f90SOURCES := bl_error.f90 bl_types.f90 test_enum.f90
F90SOURCES := microphysics_math.F90
OBJECTS := $(f90SOURCES:.f90=.o) $(F90SOURCES:.F90=.o)

vpath %.f90 . ..
vpath %.F90 . ..

ALL: test_enum

bl_error.o : bl_types.o parallel_stubs.o
parallel_stubs.o : bl_types.o

%.o: %.f90
	gfortran -c $<

%.o: %.F90
	gfortran -c $<

test_enum: $(OBJECTS)
	gfortran -o $@ $(OBJECTS)


clean:
	rm -f *.o *.mod test_enum
