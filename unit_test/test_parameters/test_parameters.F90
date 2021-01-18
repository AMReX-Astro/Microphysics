subroutine do_f90_parameters() bind(C, name="do_f90_parameters")

  use extern_probin_module

  implicit none

  !$gpu

  print *, "in Fortran"
  print *, "  eos_input_is_constant = ", eos_input_is_constant
  print *, "  test_string = ", test_string
  print *, "  dens_min = ", dens_min
  print *, "  nonaka_file = ", nonaka_file
  print *, " "

end subroutine do_f90_parameters
