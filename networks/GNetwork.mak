# this is the information to be included in a BoxLib F90 makefile
# (e.g. GMaestro.mak) to integrate the network stuff into the build
# This will add to a make variable named MICROPHYS_CORE

include $(NETWORK_TOP_DIR)/$(strip $(NETWORK_DIR))/NETWORK_REQUIRES

# network
NET_DIRS := $(NETWORK_TOP_DIR)
NET_DIRS += $(NETWORK_TOP_DIR)/$(NETWORK_DIR)

# The default integrator is specified by INTEGRATOR_DIR.
# We compile all integrators though.

INTEGRATOR_DIR ?= VODE

INT_DIRS :=

INT_DIRS += $(MICROPHYSICS_HOME)/integration/VODE
#INT_DIRS += $(MICROPHYSICS_HOME)/integration/VODE90
INT_DIRS += $(MICROPHYSICS_HOME)/integration/BS
INT_DIRS += $(MICROPHYSICS_HOME)/integration/VBDF
INT_DIRS += $(MICROPHYSICS_HOME)/integration

INT_DIRS += $(MICROPHYSICS_HOME)/integration/VODE/vode_source
#INT_DIRS += $(MICROPHYSICS_HOME)/integration/VODE90/vode90_source

# we'll assume that all integrators need the linear algebra packages
INT_DIRS += $(MICROPHYSICS_HOME)/util/

ifdef SYSTEM_BLAS
  libraries += -lblas
else
  INT_DIRS += $(MICROPHYSICS_HOME)/util/BLAS
endif

INT_DIRS += $(MICROPHYSICS_HOME)/util/LINPACK


ifeq ($(USE_RATES), TRUE)
  NET_DIRS += $(MICROPHYSICS_HOME)/rates
endif

ifeq ($(USE_SCREENING), TRUE)
  NET_DIRS += $(MICROPHYSICS_HOME)/screening
endif

ifeq ($(USE_NEUTRINOS), TRUE)
  NET_DIRS += $(MICROPHYSICS_HOME)/neutrinos
endif


MICROPHYS_CORE += $(NET_DIRS) $(INT_DIRS)
