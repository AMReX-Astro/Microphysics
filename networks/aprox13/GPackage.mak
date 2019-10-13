F90sources += actual_network.F90

ifneq ($(USE_REACT), FALSE)
ifndef SIMPLIFIED_SDC
  F90sources += actual_burner.F90
endif
  F90sources += actual_rhs.F90
endif
