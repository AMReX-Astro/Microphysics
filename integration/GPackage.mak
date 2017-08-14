ifdef SDC
  F90sources += integrator_sdc.F90
  f90sources += numerical_jacobian_sdc.f90
else
  F90sources += integrator.F90
  f90sources += numerical_jacobian.f90
endif
f90sources += integration_data.f90
f90sources += temperature_integration.f90
F90sources += rpar.F90
