f90sources += actual_network.f90

ifneq ($(USE_REACT), FALSE)
  F90sources += actual_burner.F90
  f90sources += actual_rhs.f90
endif
