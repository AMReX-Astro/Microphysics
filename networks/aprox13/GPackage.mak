f90sources += actual_network.f90

ifneq ($(USE_REACT), FALSE)
ifndef SDC
  f90sources += actual_burner.f90
endif
  f90sources += actual_rhs.f90
endif
