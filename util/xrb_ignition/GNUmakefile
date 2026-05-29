PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = FALSE

DIM        = 3

COMP       = gnu

USE_MPI    = FALSE
USE_OMP    = FALSE

USE_REACT = TRUE

EBASE = main

# define the location of the Microphysics top directory
MICROPHYSICS_HOME  := ../..

# This sets the EOS directory
EOS_DIR     := helmholtz

# This sets the network directory
NETWORK_DIR := general_null
NETWORK_INPUTS := xrb_ignition.net

USE_CONDUCTIVITY := TRUE
CONDUCTIVITY_DIR := stellar

USE_SCREENING := TRUE
SCREEN_METHOD := screen5

USE_SIMPLIFIED_SDC = TRUE

EXTERN_SEARCH += .

Bpack   := ./Make.package
Blocs   := .

include $(MICROPHYSICS_HOME)/unit_test/Make.unit_test
