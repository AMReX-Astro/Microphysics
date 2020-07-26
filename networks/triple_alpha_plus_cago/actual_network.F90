module actual_network

  use network_properties
  use amrex_fort_module, only : rt => amrex_real
  use fundamental_constants_module

  implicit none

  integer, parameter :: ihe4  = 1
  integer, parameter :: ic12  = 2
  integer, parameter :: io16  = 3
  integer, parameter :: ife56 = 4

  real(rt)        , allocatable :: bion(:)

#ifdef AMREX_USE_CUDA
  attributes(managed) :: bion
#endif

  character (len=32), parameter :: network_name = "triple_alpha_plus_cago"

  ! Rates data

  integer, parameter :: nrates = 2
  integer, parameter :: num_rate_groups = 2

  integer, parameter :: ir3a   = 1
  integer, parameter :: ircago = 2

  character (len=10), save :: reac_names(nrates)

  real, parameter :: conv_factor = MeV2eV * ev2erg * n_A

contains

  subroutine actual_network_init()

    call network_properties_init()
    allocate(bion(nspec))


    ! The following are the binding energies of the nuclei in MeV.
    ! The binding energy per unit mass (erg / g) is obtained by converting
    ! the energies in MeV to erg then multiplying by (N_A / aion) where
    ! N_A = 6.0221415e23 is Avogadro's number
    bion(ihe4)  = 28.29603_rt ! MeV / nucleus
    bion(ic12)  = 92.16294_rt ! MeV / nucleus
    bion(io16)  = 127.62093_rt ! MeV / nucleus
    bion(ife56) = 492.25389_rt ! MeV / nucleus

    ! Reaction rate names
    reac_names(ir3a)   = "3agc"   !     3 He4 --> C12
    reac_names(ircago) = "cago"   ! C12 + He4 --> O16

    !$acc update device(bion)

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
    if (allocated(bion)) then
       deallocate(bion)
    endif

  end subroutine actual_network_finalize

end module actual_network
