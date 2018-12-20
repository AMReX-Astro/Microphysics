ifdef SDC
  f90sources += vode_integrator_sdc.f90
  f90sources += vode_rhs_sdc.f90
  f90sources += vode_type_sdc.f90
else
  f90sources += vode_integrator.f90
  f90sources += vode_rhs.f90
  f90sources += vode_type.f90
endif
