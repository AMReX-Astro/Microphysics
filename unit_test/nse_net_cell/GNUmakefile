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


# define the location of the Microphysics top directory
MICROPHYSICS_HOME  := ../..

# This sets the EOS directory
EOS_DIR     := helmholtz

# This sets the network directory
NETWORK_DIR := ase

CONDUCTIVITY_DIR := stellar

SCREEN_METHOD := chabrier1998
#SCREEN_METHOD := null
#SCREEN_METHOD := chugunov2007

USE_SIMPLIFIED_SDC = TRUE

INTEGRATOR_DIR =  VODE

EXTERN_SEARCH += .

Bpack   := ./Make.package
Blocs   := .

include $(MICROPHYSICS_HOME)/unit_test/Make.unit_test
