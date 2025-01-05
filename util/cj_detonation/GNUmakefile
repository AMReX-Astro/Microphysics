PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = FALSE

DIM        = 3

COMP	   = gnu

USE_MPI    = FALSE
USE_OMP    = FALSE

USE_REACT = TRUE

EBASE = main

# define the location of the Microphysics top directory
MICROPHYSICS_HOME  := ../..

# This sets the EOS directory
EOS_DIR     := helmholtz

# This sets the network directory
NETWORK_DIR := subch_simple

CONDUCTIVITY_DIR := stellar

EXTERN_SEARCH += .

Bpack   := ./Make.package ../../unit_test/Make.package
Blocs   := .  ../../unit_test

include $(MICROPHYSICS_HOME)/unit_test/Make.unit_test

clean::
	$(RM) *.dat

