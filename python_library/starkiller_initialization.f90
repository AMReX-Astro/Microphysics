module starkiller_initialization_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  logical, save :: initialized = .false.

contains

  subroutine starkiller_initialize(probin_file)

    use amrex_fort_module, only: rt => amrex_real
    use extern_probin_module, only: small_temp, small_dens
    use microphysics_module

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    character(len=256) :: probin_file

    if (.not. initialized) then

       call runtime_init(probin_file)

       call microphysics_init(small_temp, small_dens)

       initialized = .true.

    endif

  end subroutine starkiller_initialize

end module starkiller_initialization_module
