module actual_network

  use bl_types

  implicit none

  integer, parameter :: nspec = 4
  integer, parameter :: nspec_evolve = 3
  integer, parameter :: naux  = 0

  integer, parameter :: ihe4  = 1
  integer, parameter :: ic12  = 2
  integer, parameter :: io16  = 3
  integer, parameter :: ife56 = 4

  character (len=16), save :: spec_names(nspec)
  character (len= 5), save :: short_spec_names(nspec)
  character (len= 5), save :: short_aux_names(naux)

  double precision, save :: aion(nspec), zion(nspec), ebin(nspec)

  !$acc declare create(aion, zion, ebin)
  
  character (len=22), parameter :: network_name = "triple_alpha_plus_cago"

  ! Rates data

  integer, parameter :: nrates = 2
  integer, parameter :: num_rate_groups = 2

  integer, parameter :: ir3a   = 1
  integer, parameter :: ircago = 2

  character (len=10), save :: reac_names(nrates)

contains

  subroutine actual_network_init()

    ! set the names
    spec_names(ihe4)  = "helium-4"
    spec_names(ic12)  = "carbon-12"
    spec_names(io16)  = "oxygen-16"
    spec_names(ife56) = "iron-56"

    short_spec_names(ihe4)  = "He4"
    short_spec_names(ic12)  = "C12"
    short_spec_names(io16)  = "O16"
    short_spec_names(ife56) = "Fe56"


    ! set the species properties
    aion(ihe4)  =  4.0_dp_t
    aion(ic12)  = 12.0_dp_t
    aion(io16)  = 16.0_dp_t
    aion(ife56) = 56.0_dp_t

    zion(ihe4)  =  2.0_dp_t
    zion(ic12)  =  6.0_dp_t
    zion(io16)  =  8.0_dp_t
    zion(ife56) = 26.0_dp_t

    ! our convention is that binding energy is negative.  The following are 
    ! the binding energies per unit mass (erg / g) obtained by converting
    ! the energies in MeV to erg then multiplying by (N_A / aion) where 
    ! N_A = 6.0221415e23 is Avogadro's number
    ebin(ihe4)  = -6.8253797e18_dp_t    !  28.39603 MeV / nucleon
    ebin(ic12)  = -7.4103097e18_dp_t    !  92.16294 MeV / nucleon
    ebin(io16)  = -7.6959581e18_dp_t    ! 127.62093 MeV / nucleon
    ebin(ife56) = -8.4813001e18_dp_t    ! 492.25389 MeV / nucleon

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

    ! Nothing to do here.

  end subroutine actual_network_finalize

end module actual_network
