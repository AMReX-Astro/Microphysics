ifdef SDC
  F90sources += integrator_sdc.F90
else
  F90sources += integrator.F90
endif
f90sources += integration_data.f90
F90sources += integrator_scaling.F90
