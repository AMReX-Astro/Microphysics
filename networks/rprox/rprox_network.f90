module actual_network

  use bl_types

  implicit none

  character (len=32), parameter :: network_name = "rprox"

  integer, parameter :: nspec = 10
  integer, parameter :: naux  = 0
  integer, parameter :: nrat  = 18

  character (len=16), save :: spec_names(nspec)
  character (len= 5), save :: short_spec_names(nspec)
  character (len= 5), save :: short_aux_names(naux)

  character (len=10), save :: reac_names(nrat)

  double precision, save :: aion(nspec), zion(nspec), ebin(nspec)

  integer :: ih1, ihe4, ic12, io14, io15, io16, if17, img22, is30, ini56
  integer :: irlambCNO, irag15o, irr1, irag16o, irpg16o, irpg17f, &
             irgp17f, irpg22g, irlambda2, irap14o, irs1, irlambda1, ir3a, irpg12c, &
             irwk14o, irwk17f, irwk15o, irLweak, irla2

  double precision, parameter :: MeV2erg = 1.60217646e-6
  double precision, parameter :: N_A = 6.0221415e23
  
contains

  subroutine actual_network_init

    use bl_constants_module
    use rpar_indices

    ! set the indices; ordering base on rprox.f
    ic12 = 1
    io14 = 2
    io15 = 3
    io16 = 4
    if17 = 5
    img22 = 6
    is30 = 7
    ini56 = 8
    ihe4 = 9
    ih1 = 10

    irlambCNO = 1
    irag15o = 2
    irr1 = 3
    irag16o = 4
    irpg16o = 5
    irpg17f = 6
    irgp17f = 7
    irlambda2 = 8
    irap14o = 9
    irs1 = 10
    irlambda1 = 11
    ir3a = 12
    irpg12c = 13
    irwk14o = 14
    irwk17f = 15
    irwk15o = 16
    irLweak = 17
    irla2 = 18

    ! set the names
    spec_names(ic12) = "carbon-12"
    spec_names(io14) = "oxygen-14"
    spec_names(io15) = "oxygen-15"
    spec_names(io16) = "oxygen-16"
    spec_names(if17) = "flourine-17"
    spec_names(img22) = "magnesium-22"
    spec_names(is30) = "sulfur-30"
    spec_names(ini56) = "nickel-56"
    spec_names(ihe4) = "helium-4"
    spec_names(ih1) = "hydrogen-1"

    short_spec_names(ic12) = "C12"
    short_spec_names(io14) = "O14"
    short_spec_names(io15) = "O15"
    short_spec_names(io16) = "O16"
    short_spec_names(if17) = "F17"
    short_spec_names(img22) = "Mg22"
    short_spec_names(is30) = "S30"
    short_spec_names(ini56) = "Ni56"
    short_spec_names(ihe4) = "He4"
    short_spec_names(ih1) = "H1"    

    reac_names(irlambCNO) = "rlambdaCNO"
    reac_names(irag15o) = "rag15o"
    reac_names(irr1) = "rr1"
    reac_names(irag16o) = "rag16o"
    reac_names(irpg16o) = "rpg16o"
    reac_names(irpg17f) = "rpg17f"
    reac_names(irgp17f) = "rgp17f"
    reac_names(irlambda2) = "rlambda2"
    reac_names(irap14o) = "rap14o"
    reac_names(irs1) = "rs1"
    reac_names(irlambda1) = "rlambda1"
    reac_names(ir3a) = "r3a"
    reac_names(irpg12c) = "rpg12c"
    reac_names(irwk14o) = "wk14o"
    reac_names(irwk17f) = "wk17f"
    reac_names(irwk15o) = "wk15o"
    reac_names(irLweak) = "Lweak"
    reac_names(irla2) = "la2"

    ! set the species properties
    aion(ic12) = TWELVE
    aion(io14) = 14.0_dp_t
    aion(io15) = 15.0_dp_t
    aion(io16) = 16.0_dp_t
    aion(if17) = 17.0_dp_t
    aion(img22) = 22.0_dp_t
    aion(is30) = 30.0_dp_t
    aion(ini56) = 56.0_dp_t
    aion(ihe4) = FOUR
    aion(ih1) = ONE

    zion(ic12) = SIX
    zion(io14) = EIGHT
    zion(io15) = EIGHT
    zion(io16) = EIGHT
    zion(if17) = NINE
    zion(img22) = TWELVE
    zion(is30) = 16.0_dp_t
    zion(ini56) = 28.0_dp_t
    zion(ihe4) = TWO
    zion(ih1) = ONE

    ! our convention is that binding energy is negative.  The
    ! following are the binding energy (MeV).
    ebin(ic12) =  92.16279_dp_t
    ebin(io14) =  98.7325_dp_t
    ebin(io15) = 111.9569_dp_t
    ebin(io16) = 127.6207_dp_t
    ebin(if17) = 128.2211_dp_t
    ebin(img22) = 168.5768_dp_t
    ebin(is30) = 243.6866_dp_t
    ebin(ini56) = 483.995_dp_t
    ebin(ihe4) = 28.29599_dp_t
    ebin(ih1) = ZERO

    ! convert to erg / g by multiplying by N_A / aion and converting to erg
    ebin = -ebin * N_A * MeV2erg / aion

    ! rpar is VODE's way of passing information into the RHS and                
    ! jacobian routines.  Here we initialize some indices to make               
    ! sense of what is stored in the rpar() array.                              
    call init_rpar_indices(nrat, nspec)

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

  end function network_reaction_index

end module actual_network
