f90sources += actual_network.f90

ifneq ($(USE_REACT), FALSE)
  f90sources += actual_burner.f90
  f90sources += f_rhs.f90
  f90sources += f_rhs_instantaneous.f90
  f90sources += burner_aux.f90

  USE_SCREENING = TRUE
endif
