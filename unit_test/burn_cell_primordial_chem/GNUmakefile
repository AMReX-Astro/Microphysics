PRECISION  = DOUBLE
PROFILE    = FALSE

# Set DEBUG to TRUE if debugging
DEBUG      = FALSE

DIM        = 1

COMP       = gnu
EXTRACXXFLAGS += -Wno-psabi

USE_MPI    = FALSE
USE_OMP    = FALSE
#set USE_CUDA to TRUE to compile and run on GPUs
USE_CUDA   = FALSE
USE_REACT = TRUE

# Set USE_MICROPHYSICS_DEBUG to TRUE if debugging
USE_MICROPHYSICS_DEBUG = FALSE

EBASE = main

# define the location of the Microphysics top directory
MICROPHYSICS_HOME  := ../..

# This sets the EOS directory
EOS_DIR     := primordial_chem

# This sets the network directory
NETWORK_DIR := primordial_chem

CONDUCTIVITY_DIR := stellar

INTEGRATOR_DIR =  BackwardEuler

EXTERN_SEARCH += .

Bpack   := ./Make.package
Blocs   := .

include $(MICROPHYSICS_HOME)/unit_test/Make.unit_test

ifeq ($(shell uname),Darwin)
AMREX_LINKER := ./post_link_install_name.sh $(AMREX_LINKER)
endif
