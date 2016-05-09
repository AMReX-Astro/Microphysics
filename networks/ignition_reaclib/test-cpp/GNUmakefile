PRECISION = DOUBLE
PROFILE = FALSE

DEBUG = FALSE

DIM = 1

COMP = PGI
FCOMP = PGI

USE_MPI = FALSE
USE_OMP = FALSE

USE_REACT = TRUE

USE_ACC = FALSE

# programs to be compiled
ALL: testburn.ex table

EOS_dir := helmholtz

Network_dir := ignition_reaclib

INTEGRATOR_DIR := BS

f90EXE_sources += testburn.f90

BLOCS += .

include $(CASTRO_DIR)/Exec/Make.Castro


#my_objs = $(filter-out $(objEXETempDir)/main.o, $(objForExecs))

testburn.ex: $(objForExecs)
	@echo Linking $@ ...                                                    
	$(SILENT) $(PRELINK) $(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(libraries)



