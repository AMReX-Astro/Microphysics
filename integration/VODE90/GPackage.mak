ifdef SDC
  f90sources += actual_integrator_sdc.f90
  f90sources += vode_rhs_sdc.f90
  f90sources += vode_type_sdc.f90
else
  F90sources += actual_integrator.F90
  F90sources += vode_rhs.F90
  F90sources += vode_type.F90
endif
