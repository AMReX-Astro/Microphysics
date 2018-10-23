module starkiller_initialization_module

  implicit none

contains

  subroutine starkiller_initialize(probin_file)

    use amrex_fort_module, only: rt => amrex_real
    use extern_probin_module, only: small_temp, small_dens
    use microphysics_module

    implicit none

    character(len=256) :: probin_file

    call runtime_init(probin_file)

    call microphysics_init(small_temp, small_dens)

  end subroutine starkiller_initialize

end module starkiller_initialization_module
