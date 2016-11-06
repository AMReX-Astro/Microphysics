ifdef SDC
  f90sources += integrator_sdc.f90
  f90sources += numerical_jacobian_sdc.f90
else
  f90sources += integrator.f90
  f90sources += numerical_jacobian.f90
endif
f90sources += integration_data.f90
f90sources += temperature_integration.f90
F90sources += rpar.F90
