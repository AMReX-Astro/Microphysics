subroutine init_runtime_parameters(name, namlen) bind(C, name="init_runtime_parameters")

  use amrex_fort_module, only: rt => amrex_real
  use extern_probin_module

  implicit none

  integer, intent(in) :: namlen
  integer, intent(in) :: name(namlen)

  call runtime_init(name, namlen)

end subroutine init_runtime_parameters

subroutine init_fortran_microphysics() bind(C, name="init_fortran_microphysics")

  use amrex_fort_module, only: rt => amrex_real
  use microphysics_module
  use extern_probin_module, only : small_temp, small_dens

  implicit none

  call microphysics_init(small_temp, small_dens)

end subroutine init_fortran_microphysics
