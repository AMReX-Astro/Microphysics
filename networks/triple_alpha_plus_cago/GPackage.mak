f90sources += actual_network.f90

ifneq ($(USE_REACT), FALSE)
  f90sources += actual_burner.f90
  f90sources += actual_rhs.f90
  f90sources += dydt.f90
  f90sources += screen_module.f90
  f90sources += rates_module.f90
endif
