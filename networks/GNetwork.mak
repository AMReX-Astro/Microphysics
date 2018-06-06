# this is the information to be included in a BoxLib F90 makefile
# (e.g. GMaestro.mak) to integrate the network stuff into the build
# This will add to a make variable named MICROPHYS_CORE

include $(NETWORK_TOP_DIR)/$(strip $(NETWORK_DIR))/NETWORK_REQUIRES

# network
NET_DIRS := $(NETWORK_TOP_DIR)
NET_DIRS += $(NETWORK_TOP_DIR)/$(NETWORK_DIR)

# the integrator is specified by INTEGRATOR_DIR.  We set the default to VODE
# here
INTEGRATOR_DIR ?= VODE

# Arbitrarily number the integrators so we can preprocessor tests in the code.

INTEGRATOR_NUM := -1

ifeq ($(INTEGRATOR_DIR),VODE)
  INTEGRATOR_NUM := 0
else ifeq ($(INTEGRATOR_DIR),BS)
  INTEGRATOR_NUM := 1
else ifeq ($(INTEGRATOR_DIR),VBDF)
  INTEGRATOR_NUM := 2
else ifeq ($(INTEGRATOR_DIR),VODE90)
  INTEGRATOR_NUM := 3
else ifeq ($(INTEGRATOR_DIR),XNET)
  INTEGRATOR_NUM := 4
  FPP_DEFINES += -DXNET
endif

FPP_DEFINES += -DINTEGRATOR=$(INTEGRATOR_NUM)

INT_DIRS :=

ifeq ($(INTEGRATOR_DIR),VODE)

  # Include VODE and BS

  INT_DIRS += $(MICROPHYSICS_HOME)/integration/VODE
  INT_DIRS += $(MICROPHYSICS_HOME)/integration/VODE/vode_source

  INT_DIRS += $(MICROPHYSICS_HOME)/integration/BS

else ifeq ($(INTEGRATOR_DIR),BS)

  # Include BS and VODE

  INT_DIRS += $(MICROPHYSICS_HOME)/integration/VODE
  INT_DIRS += $(MICROPHYSICS_HOME)/integration/VODE/vode_source

  INT_DIRS += $(MICROPHYSICS_HOME)/integration/BS

else

  INT_DIRS += $(MICROPHYSICS_HOME)/integration/$(INTEGRATOR_DIR)

endif

INT_DIRS += $(MICROPHYSICS_HOME)/integration

# we'll assume that all integrators need the linear algebra packages
INT_DIRS += $(MICROPHYSICS_HOME)/util/

ifndef SYSTEM_BLAS
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
