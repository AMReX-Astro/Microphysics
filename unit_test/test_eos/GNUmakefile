PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = FALSE

DIM        = 3

COMP	   = gnu

USE_MPI    = FALSE
USE_OMP    = FALSE

USE_REACT = FALSE

EBASE = main

USE_EXTRA_THERMO = TRUE

BL_NO_FORT=TRUE

# define the location of the Microphysics top directory
MICROPHYSICS_HOME  := ../..

# This sets the EOS directory
EOS_DIR     := gamma_law_general

# This sets the network directory
NETWORK_DIR := aprox19

# This isn't actually used but we need VODE to compile with CUDA
INTEGRATOR_DIR := VODE

CONDUCTIVITY_DIR := stellar

EXTERN_SEARCH += .

Bpack   := ./Make.package
Blocs   := .

include $(MICROPHYSICS_HOME)/unit_test/Make.unit_test


