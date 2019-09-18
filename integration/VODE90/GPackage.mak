ifdef SDC
  f90sources += actual_integrator_sdc.f90
  f90sources += vode_rhs_sdc.f90
  f90sources += vode_type_sdc.f90
else
  F90sources += actual_integrator.F90
  F90sources += vode_rhs.F90
  F90sources += vode_type.F90
  F90sources += cuvode_parameters.F90
endif

F90sources += vode_rpar.F90

VODE_SOURCE_DIR = $(MICROPHYSICS_HOME)/integration/VODE90/cuVODE/source/
include $(VODE_SOURCE_DIR)/GPackage.mak

INCLUDE_LOCATIONS += $(VODE_SOURCE_DIR)
VPATH_LOCATIONS   += $(VODE_SOURCE_DIR)
EXTERN_CORE       += $(VODE_SOURCE_DIR)
