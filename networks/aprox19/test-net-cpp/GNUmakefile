PRECISION = DOUBLE
PROFILE = FALSE

DEBUG = FALSE

DIM = 1

COMP = g++
FCOMP = gfortran

USE_MPI = FALSE
USE_OMP = FALSE

USE_REACT = TRUE

USE_CXX_EOS = TRUE
USE_CXX_REACTIONS := TRUE
USE_NSE := TRUE

ifeq ($(USE_CXX_REACTIONS),TRUE)
  DEFINES += -DCXX_REACTIONS
endif

# programs to be compiled
ALL: testburn.ex table

EOS_DIR := helmholtz

NETWORK_DIR := aprox19

Bpack := ./Make.package
Blocs := .

include $(MICROPHYSICS_HOME)/Make.Microphysics

#my_objs = $(filter-out $(objEXETempDir)/main.o, $(objForExecs))

testburn.ex: $(objForExecs)
	@echo Linking $@ ...                                                    
	$(SILENT) $(PRELINK) $(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(libraries)
