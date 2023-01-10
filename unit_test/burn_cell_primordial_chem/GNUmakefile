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

USE_CXX_REACTIONS ?= TRUE

ifeq ($(USE_CXX_REACTIONS),TRUE)
  DEFINES += -DCXX_REACTIONS
endif

# define the location of the CASTRO top directory
MICROPHYSICS_HOME  := ../..

# This sets the EOS directory in Castro/EOS -- note: gamma_law will not work,
# you'll need to use gamma_law_general
EOS_DIR     := primordial_chem

# This sets the network directory in Castro/Networks
NETWORK_DIR := primordial_chem

CONDUCTIVITY_DIR := stellar

INTEGRATOR_DIR =  VODE

ifeq ($(USE_CUDA), TRUE)
  INTEGRATOR_DIR := VODE
endif

EXTERN_SEARCH += .

Bpack   := ./Make.package
Blocs   := .

include $(MICROPHYSICS_HOME)/unit_test/Make.unit_test
