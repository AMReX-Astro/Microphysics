! the network module provides the information about the species we are 
! advecting:
!
! nspec      -- the number of species in the network
! nrat       -- the number of reactions in this network
!
! nevolve    -- the number of species to actually evolve
!
! aion       -- atomic number
! zion       -- proton number
! ebin       -- nuclear binding energy
!
! spec_names -- the name of the isotope
! short_spec_names -- abbreviated names
!
! reac_names -- the name of the reaction
!
!
! This module contains three routines:
!
!  network_init()         -- initialize the isotope properties
!
!  network_species_index  -- return the index of the species given its name
!
!  network_reaction_index -- return the index of the reaction given its name
!

module actual_network

  use bl_types

  implicit none

  character (len=*), parameter :: network_name = "triple_alpha_plus_cago"

  integer, parameter :: nspec = 4, nrat = 2, nevolve = 3
  integer, parameter :: naux  = 0

  character (len=16), save :: spec_names(nspec)
  character (len= 5), save :: short_spec_names(nspec)
  character (len= 5), save :: short_aux_names(naux)

  character (len=10), save :: reac_names(nrat)

  real(kind=dp_t), save :: aion(nspec), zion(nspec), ebin(nspec)

  logical, save :: network_initialized = .false.

contains

  subroutine actual_network_init()

    use network_indices
    use rpar_indices

    ! set the names
    spec_names(ihe4_)  = "helium-4"
    spec_names(ic12_)  = "carbon-12"
    spec_names(io16_)  = "oxygen-16"
    spec_names(ife56_) = "iron-56"

    short_spec_names(ihe4_)  = "He4"
    short_spec_names(ic12_)  = "C12"
    short_spec_names(io16_)  = "O16"
    short_spec_names(ife56_) = "Fe56"

    reac_names(ir3a_)   = "3agc"   !     3 He4 --> C12
    reac_names(ircago_) = "cago"   ! C12 + He4 --> O16

    ! set the species properties
    aion(ihe4_)  =  4.0_dp_t
    aion(ic12_)  = 12.0_dp_t
    aion(io16_)  = 16.0_dp_t
    aion(ife56_) = 56.0_dp_t

    zion(ihe4_)  =  2.0_dp_t
    zion(ic12_)  =  6.0_dp_t
    zion(io16_)  =  8.0_dp_t
    zion(ife56_) = 26.0_dp_t

    ! our convention is that binding energy is negative.  The following are 
    ! the binding energies per unit mass (erg / g) obtained by converting
    ! the energies in MeV to erg then multiplying by (N_A / aion) where 
    ! N_A = 6.0221415e23 is Avogadro's number
    ebin(ihe4_)  = -6.8253797e18_dp_t    !  28.39603 MeV / nucleon
    ebin(ic12_)  = -7.4103097e18_dp_t    !  92.16294 MeV / nucleon
    ebin(io16_)  = -7.6959581e18_dp_t    ! 127.62093 MeV / nucleon
    ebin(ife56_) = -8.4813001e18_dp_t    ! 492.25389 MeV / nucleon

    ! rpar is VODE's way of passing information into the RHS and
    ! jacobian routines.  Here we initialize some indices to make
    ! sense of what is stored in the rpar() array.
    call init_rpar_indices(nrat, nspec)


    ! done initializing
    network_initialized = .true.

  end subroutine actual_network_init

  
  function network_reaction_index(name)
    
    character(len=*) :: name
    integer :: network_reaction_index, n

    network_reaction_index = -1

    do n = 1, nrat
       if (name == reac_names(n)) then
          network_reaction_index = n
          exit
       endif
    enddo

    return
  end function network_reaction_index


  subroutine network_finalize()
    ! stub for MAESTRO
  end subroutine network_finalize

end module actual_network
