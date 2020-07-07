module actual_network

  use network_properties
  use amrex_fort_module, only : rt => amrex_real

  implicit none

  integer, parameter :: ihe4  = 1
  integer, parameter :: ic12  = 2
  integer, parameter :: io16  = 3
  integer, parameter :: ife56 = 4

  real(rt)        , allocatable :: ebin(:)

#ifdef AMREX_USE_CUDA
  attributes(managed) :: ebin
#endif

  character (len=32), parameter :: network_name = "triple_alpha_plus_cago"

  ! Rates data

  integer, parameter :: nrates = 2
  integer, parameter :: num_rate_groups = 2

  integer, parameter :: ir3a   = 1
  integer, parameter :: ircago = 2

  character (len=10), save :: reac_names(nrates)

contains

  subroutine actual_network_init()

    call network_properties_init()
    allocate(ebin(nspec))


    ! our convention is that binding energy is negative.  The following are
    ! the binding energies per unit mass (erg / g) obtained by converting
    ! the energies in MeV to erg then multiplying by (N_A / aion) where
    ! N_A = 6.0221415e23 is Avogadro's number
    ebin(ihe4)  = -6.8253797e18_rt    !  28.39603 MeV / nucleon
    ebin(ic12)  = -7.4103097e18_rt    !  92.16294 MeV / nucleon
    ebin(io16)  = -7.6959581e18_rt    ! 127.62093 MeV / nucleon
    ebin(ife56) = -8.4813001e18_rt    ! 492.25389 MeV / nucleon

    ! Reaction rate names
    reac_names(ir3a)   = "3agc"   !     3 He4 --> C12
    reac_names(ircago) = "cago"   ! C12 + He4 --> O16

    !$acc update device(aion, zion, ebin)

  end subroutine actual_network_init



  function network_reaction_index(name)

    character(len=*) :: name
    integer :: network_reaction_index, n

    network_reaction_index = -1

    do n = 1, nrates
       if (name == reac_names(n)) then
          network_reaction_index = n
          exit
       endif
    enddo

    return
  end function network_reaction_index



  subroutine actual_network_finalize

    implicit none

    call network_properties_finalize()
    if (allocated(ebin)) then
       deallocate(ebin)
    endif

  end subroutine actual_network_finalize

end module actual_network
