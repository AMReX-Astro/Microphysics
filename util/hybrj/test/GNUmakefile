PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = FALSE

DIM        = 3

COMP       = gnu

USE_MPI    = FALSE
USE_OMP    = FALSE

USE_REACT = TRUE
USE_SIMPLIFIED_SDC = TRUE
USE_NSE_NET = TRUE
EBASE = main

# define the location of the Microphysics top directory
MICROPHYSICS_HOME  := ../../..

# This sets the EOS directory
EOS_DIR     := helmholtz

# This sets the network directory
NETWORK_DIR := aprox21

SCREEN_METHOD = chabrier1998

CONDUCTIVITY_DIR := stellar

INTEGRATOR_DIR =  VODE

ifeq ($(USE_CUDA), TRUE)
  INTEGRATOR_DIR := VODE
endif

EXTERN_SEARCH += .

Bpack   := ../Make.package ./Make.package
Blocs   := ../ .

include $(MICROPHYSICS_HOME)/unit_test/Make.unit_test
