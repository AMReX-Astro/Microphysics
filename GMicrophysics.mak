# A set of useful macros for putting together one of the initial model
# generator routines in the FBoxLib framework

# include the main Makefile stuff
include $(FBOXLIB_HOME)/Tools/F_mk/GMakedefs.mak


#-----------------------------------------------------------------------------
# core FBoxLib directories
FBOXLIB_CORE := Src/BaseLib


#-----------------------------------------------------------------------------
ifdef ACC
  ifeq ($(ACC), t)
    FPP_DEFINES += -DACC -DUSE_ACC
  endif
endif

ifdef SDC
  ifeq ($(SDC), t)
    FPP_DEFINES += -DSDC
  endif
endif

FPP_DEFINES += -DAMREX_DEVICE=""

#-----------------------------------------------------------------------------
# EOS
EOS_TOP_DIR := $(MICROPHYSICS_HOME)/EOS

ifdef EXTRA_THERMO
  ifeq ($(EXTRA_THERMO), t)
    FPP_DEFINES += -DEXTRA_THERMO
  endif
endif

# the helmeos has a table
ifeq ($(findstring helmholtz, $(EOS_DIR)), helmholtz)
  EOS_PATH := $(EOS_TOP_DIR)/helmholtz
  ALL: table
endif

table:
	@if [ ! -f helm_table.dat ]; then echo ${bold}Linking helm_table.dat${normal}; ln -s $(EOS_PATH)/helm_table.dat .;  fi

EOS_DIRS := $(EOS_TOP_DIR)
EOS_DIRS += $(EOS_TOP_DIR)/$(EOS_DIR)

MICROPHYS_CORE += $(EOS_DIRS) 


#-----------------------------------------------------------------------------
# network stuff -- specify your particlar network via NETWORK_DIR
# this will increment MICROPHYS_CORE
NETWORK_TOP_DIR := $(MICROPHYSICS_HOME)/networks
include $(NETWORK_TOP_DIR)/GNetwork.mak

# URCA network has tables
ifeq ($(findstring URCA-simple, $(NETWORK_DIR)), URCA-simple)
  ALL: urcatables
endif

urcatables:
	@if [ ! -f 23Ne-23Na_betadecay.dat ]; then echo ${bold}Linking 23Ne-23Na_betadecay.dat${normal}; ln -s $(NETWORK_TOP_DIR)/$(NETWORK_DIR)/23Ne-23Na_betadecay.dat .;  fi
	@if [ ! -f 23Na-23Ne_electroncapture.dat ]; then echo ${bold}Linking 23Na-23Ne_electroncapture.dat${normal}; ln -s $(NETWORK_TOP_DIR)/$(NETWORK_DIR)/23Na-23Ne_electroncapture.dat .;  fi




#-----------------------------------------------------------------------------
# unit testing directories
UNIT_DIR := $(MICROPHYSICS_HOME)/unit_test
UNIT_DIR += $(MICROPHYSICS_HOME)/interfaces
UNIT_DIR += $(TEST_DIR)     # set by the test itself


#-----------------------------------------------------------------------------
# core FBoxLib directories
Fmpack := $(foreach dir, $(FBOXLIB_CORE), $(FBOXLIB_HOME)/$(dir)/GPackage.mak)
Fmlocs := $(foreach dir, $(FBOXLIB_CORE), $(FBOXLIB_HOME)/$(dir))
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

PROBIN_TEMPLATE := $(MICROPHYSICS_HOME)/unit_test/dummy.probin.template
PROBIN_PARAMETER_DIRS += $(MICROPHYSICS_HOME)/unit_test/  
EXTERN_PARAMETER_DIRS += $(MICROPHYS_CORE)


PROBIN_PARAMETERS := $(shell $(FBOXLIB_HOME)/Tools/F_scripts/findparams.py $(PROBIN_PARAMETER_DIRS))
EXTERN_PARAMETERS := $(shell $(FBOXLIB_HOME)/Tools/F_scripts/findparams.py $(EXTERN_PARAMETER_DIRS))

probin.f90: $(PROBIN_PARAMETERS) $(EXTERN_PARAMETERS) $(PROBIN_TEMPLATE)
	@echo " "
	@echo "${bold}WRITING probin.f90${normal}"
	$(FBOXLIB_HOME)/Tools/F_scripts/write_probin.py \
           -t $(PROBIN_TEMPLATE) -o probin.f90 -n probin \
           --pa "$(PROBIN_PARAMETERS)" --pb "$(EXTERN_PARAMETERS)"
	@echo " "



# vpath defines the directories to search for the source files

#  VPATH_LOCATIONS to first search in the problem directory
#  Note: GMakerules.mak will include '.' at the start of the
VPATH_LOCATIONS += $(Fmlocs)


# list of directories to put in the Fortran include path
FINCLUDE_LOCATIONS += $(Fmincs)

#-----------------------------------------------------------------------------
# build_info stuff
deppairs: build_info.f90

build_info.f90: 
	@echo " "
	@echo "${bold}WRITING build_info.f90${normal}"
	$(FBOXLIB_HOME)/Tools/F_scripts/makebuildinfo.py \
           --modules "$(Fmdirs) $(MICROPHYS_CORE) $(UNIT_DIR)" \
           --FCOMP "$(COMP)" \
           --FCOMP_version "$(FCOMP_VERSION)" \
           --f90_compile_line "$(COMPILE.f90)" \
           --f_compile_line "$(COMPILE.f)" \
           --C_compile_line "$(COMPILE.c)" \
           --link_line "$(LINK.f90)" \
           --fboxlib_home "$(FBOXLIB_HOME)" \
           --source_home "$(MICROPHYSICS_HOME)" \
           --network "$(NETWORK_DIR)" \
           --integrator "$(INTEGRATOR_DIR)" \
           --eos "$(EOS_DIR)"
	@echo " "

$(odir)/build_info.o: build_info.f90
	$(COMPILE.f90) $(OUTPUT_OPTION) build_info.f90
	rm -f build_info.f90



#-----------------------------------------------------------------------------
# for debugging.  To see the value of a Makefile variable,
# e.g. Fmlocs, simply do "make print-Fmlocs".  This will
# print out the value.
print-%: ; @echo $* is $($*)


#-----------------------------------------------------------------------------
# cleaning.  Add more actions to 'clean' and 'realclean' to remove
# probin.f90 and build_info.f90 -- this is where the '::' in make comes
# in handy
clean ::
	$(RM) probin.f90
	$(RM) build_info.f90


realclean ::
	$(RM) helm_table.dat
