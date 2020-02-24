subroutine do_f90_parameters() bind(C, name="do_f90_parameters")

  use extern_probin_module

  implicit none

  !$gpu

  print *, "in Fortran"
  print *, "eos_input_is_constant = ", eos_input_is_constant

end subroutine do_f90_parameters
