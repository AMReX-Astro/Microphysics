
INTEGRATOR_DIR ?= VODE

ifeq ($(INTEGRATOR_DIR), VODE)
  INTEGRATOR_PATH := $(MICROPHYSICS_DIR)/integration/VODE
endif

FINCLUDE_LOCATIONS    += $(INTEGRATOR_PATH)
VPATH_LOCATIONS       += $(INTEGRATOR_PATH)
EXTERN_PARAMETER_DIRS += $(INTEGRATOR_PATH)

include $(INTEGRATION_PATH)/GPackage.mak

f90sources += integrator.f90
f90sources += integration_data.f90
f90sources += temperature_integration.f90
f90sources += numerical_jacobian.f90
