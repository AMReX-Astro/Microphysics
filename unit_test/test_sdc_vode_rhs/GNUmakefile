PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = FALSE

DIM        = 3

COMP	   = gnu

USE_MPI    = FALSE
USE_OMP    = FALSE

USE_REACT = TRUE
USE_SIMPLIFIED_SDC = TRUE

EBASE = main

USE_EXTRA_THERMO = TRUE

BL_NO_FORT = TRUE

# define the location of the CASTRO top directory
MICROPHYSICS_HOME  := ../..

# This sets the EOS directory in Castro/EOS -- note: gamma_law will not work,
# you'll need to use gamma_law_general
EOS_DIR     := helmholtz

NETWORK_DIR := aprox19

CONDUCTIVITY_DIR := stellar
INTEGRATOR_DIR := VODE
DEFINES += -DSDC_EVOLVE_ENERGY

EXTERN_SEARCH += . ..

Bpack   := ./Make.package
Blocs   := .

include $(MICROPHYSICS_HOME)/unit_test/Make.unit_test

