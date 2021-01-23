subroutine init_unit_test_fortran(name, namlen) bind(C, name="init_unit_test_fortran")

  use amrex_fort_module, only: rt => amrex_real
  use extern_probin_module
  use microphysics_module

  implicit none

  integer, intent(in) :: namlen
  integer, intent(in) :: name(namlen)

  call runtime_init(name, namlen)

end subroutine init_unit_test_fortran

subroutine init_microphysics_fortran() bind(C, name="init_microphysics_fortran")

  use amrex_fort_module, only: rt => amrex_real
  use extern_probin_module
  use microphysics_module

  implicit none

  call microphysics_init(small_temp, small_dens)

end subroutine init_microphysics_fortran
