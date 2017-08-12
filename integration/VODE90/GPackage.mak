ifdef SDC
  f90sources += actual_integrator_sdc.f90
  f90sources += vode_rhs_sdc.f90
  f90sources += vode_type_sdc.f90
else
  f90sources += actual_integrator.f90
  f90sources += vode_rhs.f90
  f90sources += vode_type.f90
endif

VODE_SOURCE_DIR = $(MICROPHYSICS_HOME)/integration/VODE90/vode_source/
include $(VODE_SOURCE_DIR)/GPackage.mak

INCLUDE_LOCATIONS += $(VODE_SOURCE_DIR)
VPATH_LOCATIONS   += $(VODE_SOURCE_DIR)
EXTERN_CORE       += $(VODE_SOURCE_DIR)
