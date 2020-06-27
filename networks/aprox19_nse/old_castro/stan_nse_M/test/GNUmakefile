NDEBUG := t
NDEBUG :=  
MPI    :=
OMP    :=

COMP := gfortran
#COMP := Intel

MKVERBOSE := t

programs += test.f90

FPARALLEL := ../../../../..

include $(BOXLIB_HOME)/Tools/F_mk/GMakedefs.mak

EOS_DIR := extern/EOS/helmeos
NETWORK_DIR := $(UCSC_ASTRO_DIR)/networks/stan_nse_M

BOXLIB_CORE := Src/F_BaseLib
# EXTERN_CORE := extern/BLAS \
# 	       extern/VODE \
# 	       $(EOS_DIR) \
# 	       $(NETWORK_DIR) 
EXTERN_CORE := $(NETWORK_DIR)

#Fmincludes := $(EOS_DIR)

Fmpack := $(foreach dir, $(EXTERN_CORE), $(dir)/Make.package)
Fmlocs := $(foreach dir, $(EXTERN_CORE), $(dir))
Fmincs := $(foreach dir, $(Fmincludes), $(FPARALLEL)/$(dir))

Fmpack += $(foreach dir, $(BOXLIB_CORE), $(BOXLIB_HOME)/$(dir)/GPackage.mak)
Fmlocs += $(foreach dir, $(BOXLIB_CORE), $(BOXLIB_HOME)/$(dir))

include $(Fmpack)
VPATH_LOCATIONS += $(Fmlocs) $(NETWORK_DIR)
FINCLUDE_LOCATIONS += $(Fmincs) $(NETWORK_DIR)


all: $(pnames)

include $(BOXLIB_HOME)/Tools/F_mk/GMakerules.mak

%.$(suf).exe:%.f90 $(objects)
ifdef MKVERBOSE
	$(LINK.f90) -o $@ $< $(objects) $(libraries)
else
	@echo "Linking $@ ... "
	@$(LINK.f90) -o $@ $< $(objects) $(libraries)
endif



