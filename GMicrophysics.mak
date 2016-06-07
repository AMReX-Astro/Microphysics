# A set of useful macros for putting together one of the initial model
# generator routines

# include the main Makefile stuff
include $(BOXLIB_HOME)/Tools/F_mk/GMakedefs.mak


#-----------------------------------------------------------------------------
# core BoxLib directories
BOXLIB_CORE := Src/F_BaseLib


#-----------------------------------------------------------------------------
# EOS
EOS_TOP_DIR := $(MICROPHYSICS_DIR)/eos

# the helmeos has a table
ifeq ($(findstring helmholtz, $(EOS_DIR)), helmholtz)
  EOS_PATH := $(EOS_TOP_DIR)/helmholtz
  ALL: table
endif

table:
	@if [ ! -f helm_table.dat ]; then echo ${bold}Linking helm_table.dat${normal}; ln -s $(EOS_PATH)/helm_table.dat .;  fi

EOS_DIRS := $(EOS_TOP_DIR)
EOS_DIRS += $(EOS_TOP_DIR)/$(EOS_DIR)


#-----------------------------------------------------------------------------
# network
NETWORK_TOP_DIR := $(MICROPHYSICS_DIR)/networks
NET_DIRS := $(NETWORK_TOP_DIR)
NET_DIRS += $(NETWORK_TOP_DIR)/$(NETWORK_DIR)

ifneq ($(USE_SCREENING), FALSE)
  NET_DIRS += $(MICROPHYSICS_DIR)/screening
endif

ifneq ($(USE_RATES), FALSE)
  NET_DIRS += $(MICROPHYSICS_DIR)/rates
endif

ifneq ($(USE_NEUTRINOS), FALSE)
  NET_DIRS += $(MICROPHYSICS_DIR)/neutrinos
endif


#-----------------------------------------------------------------------------
# integrator

# the integrator is specified by INTEGRATOR_DIR.  We set the default to VODE
# here
INTEGRATOR_DIR ?= VODE
INT_DIRS := $(MICROPHYSICS_DIR)/integration
INT_DIRS += $(MICROPHYSICS_DIR)/integration/$(INTEGRATOR_DIR)

# we'll assume that all integrators need the linear algebra packages
INT_DIRS += $(MICROPHYSICS_DIR)/util/
INT_DIRS += $(MICROPHYSICS_DIR)/util/BLAS
INT_DIRS += $(MICROPHYSICS_DIR)/util/LINPACK


# add in the network, EOS, and conductivity
MICROPHYS_CORE += $(EOS_DIRS) $(NET_DIRS) $(INT_DIRS)


#-----------------------------------------------------------------------------
# unit testing directories
UNIT_DIR := $(MICROPHYSICS_DIR)/unit_test
UNIT_DIR += $(MICROPHYSICS_DIR)/interfaces
UNIT_DIR += $(TEST_DIR)     # set by the test itself


#-----------------------------------------------------------------------------
# core BoxLib directories
Fmpack := $(foreach dir, $(BOXLIB_CORE), $(BOXLIB_HOME)/$(dir)/GPackage.mak)
Fmlocs := $(foreach dir, $(BOXLIB_CORE), $(BOXLIB_HOME)/$(dir))
Fmincs :=

# auxillary directories
Fmpack += $(foreach dir, $(MICROPHYS_CORE), $(dir)/GPackage.mak)
Fmpack += $(foreach dir, $(UNIT_DIR), $(dir)/GPackage.mak)

Fmlocs += $(foreach dir, $(MICROPHYS_CORE), $(dir))
Fmlocs += $(foreach dir, $(UNIT_DIR), $(dir))


# include the necessary GPackage.mak files that define this setup
include $(Fmpack)


# we need a probin.f90, since the various microphysics routines can
# have runtime parameters
f90sources += probin.f90

PROBIN_TEMPLATE := $(MICROPHYSICS_DIR)/unit_test/dummy.probin.template
PROBIN_PARAMETER_DIRS = $(MICROPHYSICS_DIR)/unit_test/
EXTERN_PARAMETER_DIRS += $(MICROPHYS_CORE)


PROBIN_PARAMETERS := $(shell $(BOXLIB_HOME)/Tools/F_scripts/findparams.py $(PROBIN_PARAMETER_DIRS))
EXTERN_PARAMETERS := $(shell $(BOXLIB_HOME)/Tools/F_scripts/findparams.py $(EXTERN_PARAMETER_DIRS))

probin.f90: $(PROBIN_PARAMETERS) $(EXTERN_PARAMETERS) $(PROBIN_TEMPLATE)
	@echo " "
	@echo "${bold}WRITING probin.f90${normal}"
	$(BOXLIB_HOME)/Tools/F_scripts/write_probin.py \
           -t $(PROBIN_TEMPLATE) -o probin.f90 -n probin \
           --pa "$(PROBIN_PARAMETERS)" --pb "$(EXTERN_PARAMETERS)"
	@echo " "



# vpath defines the directories to search for the source files

#  VPATH_LOCATIONS to first search in the problem directory
#  Note: GMakerules.mak will include '.' at the start of the
VPATH_LOCATIONS += $(Fmlocs)


# list of directories to put in the Fortran include path
FINCLUDE_LOCATIONS += $(Fmincs)


init_1d.$(suf).exe: $(objects)
	$(LINK.f90) -o init_1d.$(suf).exe $(objects) $(libraries)
	@echo SUCCESS


# include the fParallel Makefile rules
include $(BOXLIB_HOME)/Tools/F_mk/GMakerules.mak


#-----------------------------------------------------------------------------
# for debugging.  To see the value of a Makefile variable,
# e.g. Fmlocs, simply do "make print-Fmlocs".  This will
# print out the value.
print-%: ; @echo $* is $($*)


#-----------------------------------------------------------------------------
# cleaning.  Add more actions to 'clean' and 'realclean' to remove
# probin.f90 and build_info.f90 -- this is where the '::' in make comes
# in handy
clean::
	$(RM) probin.f90
	$(RM) build_info.f90
