ifeq ($(USE_SDC), TRUE)
  F90sources += actual_integrator_sdc.F90
else
  F90sources += actual_integrator.F90
endif
