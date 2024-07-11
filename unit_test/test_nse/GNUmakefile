PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = FALSE

DIM        = 3

COMP	   = gnu

USE_MPI    = FALSE
USE_OMP    = FALSE

USE_REACT = TRUE

USE_NSE_NET = TRUE
EBASE = main

BL_NO_FORT = TRUE

# define the location of the Microphysics top directory
MICROPHYSICS_HOME  := ../..

# This sets the EOS directory
EOS_DIR     := helmholtz

# This sets the network directory
NETWORK_DIR := aprox21

SCREEN_METHOD := chabrier1998

USE_SIMPLIFIED_SDC = TRUE

CONDUCTIVITY_DIR := stellar

INTEGRATOR_DIR =  VODE

EXTERN_SEARCH += . ..

Bpack   := ./Make.package
Blocs   := .

include $(MICROPHYSICS_HOME)/unit_test/Make.unit_test
