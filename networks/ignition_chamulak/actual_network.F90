module actual_network

  use network_properties
  use amrex_fort_module, only : rt => amrex_real
  implicit none

  character (len=32), parameter :: network_name = "ignition_chamulak"

  ! M12_chamulak is the effective number of C12 nuclei destroyed per
  ! reaction
  real(rt)        , parameter :: M12_chamulak = 2.93e0_rt

  integer, parameter :: ic12  = 1
  integer, parameter :: io16  = 2
  integer, parameter :: iash  = 3

  real(rt)        , allocatable :: bion(:), mion(:), wion(:)

#ifdef AMREX_USE_CUDA
  attributes(managed) :: bion, mion, wion
#endif

  integer, parameter :: nrates = 1
  integer, parameter :: num_rate_groups = 4

contains

  subroutine actual_network_init

    implicit none

    call network_properties_init()
    allocate(bion(nspec))
    allocate(mion(nspec))
    allocate(wion(nspec))

    ! the ash from C12 burning according to Chamulak et al. is a mixture
    ! of C13, O16, Ne20, and Na23.   Ne20 + alpha results 60% of the time,
    ! while Na23 + p is the other 40%.   Fusing 6 C12 will result in
    ! 1.2 Na23, 1.2 O16 (from the alpha), 1.8 Ne20, and 1.8 C13.
    ! The ash state will have an A and Z corresponding to this mixture.

  end subroutine actual_network_init



  subroutine actual_network_finalize

    implicit none

    call network_properties_finalize()

    if (allocated(bion)) then
       deallocate(bion)
    endif
    if (allocated(mion)) then
       deallocate(mion)
    endif
    if (allocated(wion)) then
       deallocate(wion)
    endif

  end subroutine actual_network_finalize

end module actual_network
