f90sources += xrb_simple_network.f90
f90sources += network_indices.f90

ifneq ($(USE_REACT), FALSE)
  f90sources += xrb_simple_burner.f90
  f90sources += burner_aux.f90
  f90sources += f_rhs.f90
endif
