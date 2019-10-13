ifdef SIMPLIFIED_SDC
  f90sources += vode_integrator_sdc.f90
  f90sources += vode_rhs_sdc.f90
  f90sources += vode_type_sdc.f90
else
  f90sources += vode_integrator.f90
  F90sources += vode_rhs.F90
  F90sources += vode_type.F90
endif

F90sources += vode_rpar.F90
