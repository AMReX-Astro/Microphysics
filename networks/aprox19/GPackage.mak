f90sources += actual_network.f90

ifneq ($(USE_REACT), FALSE)
  f90sources += actual_burner.f90
  f90sources += actual_rhs.f90
endif
