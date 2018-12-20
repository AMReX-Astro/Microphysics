ifdef SDC
  F90sources += bs_integrator_sdc.F90
  F90sources += bs_type_sdc.F90
  F90sources += bs_rhs_sdc.F90
  F90sources += bs_jac_sdc.F90
else
  F90sources += bs_integrator.F90
  F90sources += bs_type.F90
  F90sources += bs_rhs.F90
  F90sources += bs_jac.F90
endif
F90sources += stiff_ode.F90
