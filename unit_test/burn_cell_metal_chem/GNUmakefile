PRECISION  = DOUBLE
PROFILE    = FALSE

# Set DEBUG to TRUE if debugging
DEBUG      = TRUE

DIM        = 1

COMP       = gnu

USE_MPI    = FALSE
USE_OMP    = FALSE
#set USE_CUDA to TRUE to compile and run on GPUs
USE_CUDA   = FALSE
USE_REACT = TRUE

# Set USE_MICROPHYSICS_DEBUG to TRUE if debugging
USE_MICROPHYSICS_DEBUG = TRUE

EBASE = main

# define the location of the Microphysics top directory
MICROPHYSICS_HOME  := ../..

# This sets the EOS directory
EOS_DIR     := metal_chem

# This sets the network directory
NETWORK_DIR := metal_chem

CONDUCTIVITY_DIR := stellar

INTEGRATOR_DIR =  VODE

EXTERN_SEARCH += .

Bpack   := ./Make.package
Blocs   := .

include $(MICROPHYSICS_HOME)/unit_test/Make.unit_test
