f90sources += actual_integrator.f90
f90sources += vode_convert.f90
f90sources += vode_rhs.f90

VODE_SOURCE_DIR = $(MICROPHYSICS_DIR)/integration/VODE/vode_source

FINCLUDE_LOCATIONS    += $(VODE_SOURCE_DIR)
VPATH_LOCATIONS       += $(VODE_SOURCE_DIR)
EXTERN_PARAMETER_DIRS += $(VODE_SOURCE_DIR)

include $(VODE_SOURCE_DIR)/GPackage.mak
