f90sources += actual_network.f90
f90sourcse += actual_network_data.f90

ifneq ($(USE_REACT), FALSE)
  f90sources += actual_burner.f90
  f90sources += actual_burner_data.f90
  f90sources += actual_rhs.f90
  f90sources += rates_module.f90
endif
