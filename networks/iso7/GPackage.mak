f90sources += actual_network.f90
f90sources += actual_network_data.f90

ifneq ($(USE_REACT), FALSE)
  f90sources += actual_burner.f90
  f90sources += actual_rhs.f90

  USE_RATES       = TRUE
  USE_SCREENING   = TRUE
  USE_NEUTRINOS   = TRUE
endif
