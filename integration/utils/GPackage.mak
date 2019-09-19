ifdef SDC
  f90sources += numerical_jacobian_sdc.f90
else
  F90sources += numerical_jacobian.F90
endif
F90sources += network_rhs.F90
F90sources += temperature_integration.F90
F90sources += jacobian_sparsity.F90
f90sources += nonaka_plot.f90
