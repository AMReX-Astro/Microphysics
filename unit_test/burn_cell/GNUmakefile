PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = FALSE

DIM        = 3

COMP	   = gnu

USE_MPI    = FALSE
USE_OMP    = FALSE

USE_REACT = TRUE

EBASE = main

BL_NO_FORT = TRUE

# For DVODE parity checks on the Robertson system, evolve species only.
DEFINES += -DMICROPHYSICS_STRANG_SPECIES_ONLY
DEFINES += -DUNIT_TEST_ENABLE_VODE_CONTINUATION

# define the location of the Microphysics top directory
MICROPHYSICS_HOME  := ../..

# This sets the EOS directory
EOS_DIR     := helmholtz

# This sets the network directory
NETWORK_DIR := aprox13
SCREEN_METHOD := screen5
CONDUCTIVITY_DIR := stellar

INTEGRATOR_DIR =  VODE

EXTERN_SEARCH += . ..

Bpack   := ./Make.package
Blocs   := .

include $(MICROPHYSICS_HOME)/unit_test/Make.unit_test
