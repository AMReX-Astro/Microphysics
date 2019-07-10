module actual_network

  implicit none

  ! M12_chamulak is the effective number of C12 nuclei destroyed per
  ! reaction
  double precision, parameter :: M12_chamulak = 2.93d0

  integer, parameter :: nspec = 3

  ! we are only explicitly evolving C12
  integer, parameter :: nspec_evolve = 1
  integer, parameter :: naux  = 0

  integer, parameter :: ic12  = 1
  integer, parameter :: io16  = 2
  integer, parameter :: iash  = 3

  character (len=16), save :: spec_names(nspec)
  character (len= 5), save :: short_spec_names(nspec)
  character (len= 5), save :: short_aux_names(naux)

  double precision, allocatable :: aion(:), zion(:), nion(:)
  double precision, allocatable :: bion(:), mion(:), wion(:)

#ifdef AMREX_USE_CUDA
  attributes(managed) :: aion, zion, nion, bion, mion, wion
#endif

  integer, parameter :: nrates = 1
  integer, parameter :: num_rate_groups = 4

contains

  subroutine actual_network_init

    implicit none

    spec_names(ic12)  = "carbon-12"
    spec_names(io16)  = "oxygen-16"
    spec_names(iash)  = "ash"

    short_spec_names(ic12) = "C12"
    short_spec_names(io16) = "O16"
    short_spec_names(iash) = "ash"

    allocate(aion(nspec))
    allocate(zion(nspec))
    allocate(nion(nspec))
    allocate(bion(nspec))
    allocate(mion(nspec))
    allocate(wion(nspec))

    ! the ash from C12 burning according to Chamulak et al. is a mixture
    ! of C13, O16, Ne20, and Na23.   Ne20 + alpha results 60% of the time,
    ! while Na23 + p is the other 40%.   Fusing 6 C12 will result in
    ! 1.2 Na23, 1.2 O16 (from the alpha), 1.8 Ne20, and 1.8 C13.
    ! The ash state will have an A and Z corresponding to this mixture.
    aion(ic12)  = 12.0d0
    aion(io16)  = 16.0d0
    aion(iash)  = 18.0d0

    zion(ic12)  = 6.0d0
    zion(io16)  = 8.0d0
    zion(iash)  = 8.8d0

  end subroutine actual_network_init



  subroutine actual_network_finalize

    implicit none

    if (allocated(aion)) then
       deallocate(aion)
    endif
    if (allocated(zion)) then
       deallocate(zion)
    endif
    if (allocated(nion)) then
       deallocate(nion)
    endif
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
