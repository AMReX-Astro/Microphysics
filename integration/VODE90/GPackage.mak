ifdef SDC
  f90sources += vode90_integrator_sdc.f90
  f90sources += vode90_rhs_sdc.f90
  f90sources += vode90_type_sdc.f90
else
  f90sources += vode90_integrator.f90
  f90sources += vode90_rhs.f90
  f90sources += vode90_type.f90
endif

VODE90_SOURCE_DIR = $(MICROPHYSICS_HOME)/integration/VODE90/vode90_source/
include $(VODE90_SOURCE_DIR)/GPackage.mak

INCLUDE_LOCATIONS += $(VODE90_SOURCE_DIR)
VPATH_LOCATIONS   += $(VODE90_SOURCE_DIR)
EXTERN_CORE       += $(VODE90_SOURCE_DIR)
