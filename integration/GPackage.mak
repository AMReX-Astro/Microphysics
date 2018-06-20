ifdef SDC
  F90sources += integrator_sdc.F90
  f90sources += numerical_jacobian_sdc.f90
else
  F90sources += integrator.F90
  F90sources += numerical_jacobian.F90
endif
f90sources += integration_data.f90
F90sources += temperature_integration.F90
F90sources += rpar.F90
