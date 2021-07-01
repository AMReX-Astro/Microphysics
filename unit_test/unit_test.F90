subroutine init_microphysics_fortran() bind(C, name="init_microphysics_fortran")

  use amrex_fort_module, only: rt => amrex_real
  use extern_probin_module
  use microphysics_module

  implicit none

  call microphysics_init(small_temp, small_dens)

end subroutine init_microphysics_fortran
