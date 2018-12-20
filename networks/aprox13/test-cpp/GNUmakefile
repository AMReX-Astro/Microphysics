PRECISION = DOUBLE
PROFILE = FALSE

DEBUG = FALSE

DIM = 1

COMP = g++
FCOMP = gfortran

USE_MPI = FALSE
USE_OMP = FALSE

USE_REACT = TRUE

# programs to be compiled
ALL: testburn.ex table

EOS_dir := helmholtz

Network_dir := aprox13

INTEGRATOR_DIR := VBDF

f90EXE_sources += testburn.f90
f90EXE_sources += testjacobian.f90

BLOCS += .

include $(CASTRO_DIR)/Exec/Make.Castro


#my_objs = $(filter-out $(objEXETempDir)/main.o, $(objForExecs))

testburn.ex: $(objForExecs)
	@echo Linking $@ ...                                                    
	$(SILENT) $(PRELINK) $(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(libraries)



