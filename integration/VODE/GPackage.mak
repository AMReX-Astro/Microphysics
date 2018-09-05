ifdef SDC
  F90sources += vode_integrator_sdc.F90
  f90sources += vode_rhs_sdc.f90
  f90sources += vode_type_sdc.f90
else
  F90sources += vode_integrator.F90
  f90sources += vode_rhs.f90
  f90sources += vode_type.f90
endif
