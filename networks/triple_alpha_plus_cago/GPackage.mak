F90sources += actual_network.F90

ifneq ($(USE_REACT), FALSE)
  F90sources += actual_burner.F90
  F90sources += actual_rhs.F90
  F90sources += dydt.F90
  F90sources += screen_module.F90
  F90sources += rates_module.F90
endif
