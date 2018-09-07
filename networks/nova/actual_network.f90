module actual_network
  use physical_constants, only: ERG_PER_MeV
  use bl_types
  
  implicit none

  public num_rate_groups

  double precision, parameter :: avo = 6.0221417930d23
  double precision, parameter :: c_light = 2.99792458d10
  double precision, parameter :: enuc_conv2 = -avo*c_light*c_light

  double precision, parameter :: ev2erg  = 1.60217648740d-12
  double precision, parameter :: mev2erg = ev2erg*1.0d6
  double precision, parameter :: mev2gr  = mev2erg/c_light**2

  double precision, parameter :: mass_neutron  = 1.67492721184d-24
  double precision, parameter :: mass_proton   = 1.67262163783d-24
  double precision, parameter :: mass_electron = 9.10938215450d-28

  integer, parameter :: nrates = 0
  integer, parameter :: num_rate_groups = 4

  ! Evolution and auxiliary
  integer, parameter :: nspec_evolve = 0
  integer, parameter :: naux  = 0

  ! Number of nuclear species in the network
  integer, parameter :: nspec = 0

  ! Number of reaclib rates
  integer, parameter :: nrat_reaclib = 0
  
  ! Number of tabular rates
  integer, parameter :: nrat_tabular = 0

  ! Binding Energies Per Nucleon (MeV)
  double precision :: ebind_per_nucleon(nspec)

  ! aion: Nucleon mass number A
  ! zion: Nucleon atomic number Z
  ! nion: Nucleon neutron number N
  ! bion: Binding Energies (ergs)

  ! Nuclides

  ! Reactions

  ! reactvec indices
  integer, parameter :: i_rate        = 1
  integer, parameter :: i_drate_dt    = 2
  integer, parameter :: i_scor        = 3
  integer, parameter :: i_dscor_dt    = 4
  integer, parameter :: i_dqweak      = 5
  integer, parameter :: i_epart       = 6

  character (len=16), save :: spec_names(nspec) 
  character (len= 5), save :: short_spec_names(nspec)
  character (len= 5), save :: short_aux_names(naux)

  double precision :: aion(nspec), zion(nspec), bion(nspec)
  double precision :: nion(nspec), mion(nspec), wion(nspec)

  !$acc declare create(aion, zion, bion, nion, mion, wion)

contains

  subroutine actual_network_init()
    
    implicit none
    
    integer :: i







    do i = 1, nspec
       bion(i) = ebind_per_nucleon(i) * aion(i) * ERG_PER_MeV
    end do

    ! Set the mass
    mion(:) = nion(:) * mass_neutron + zion(:) * (mass_proton + mass_electron) &
         - bion(:)/(c_light**2)

    ! Molar mass
    wion(:) = avo * mion(:)

    ! Common approximation
    !wion(:) = aion(:)

    !$acc update device(aion, zion, bion, nion, mion, wion)
  end subroutine actual_network_init

  subroutine actual_network_finalize()
    ! STUB FOR MAESTRO
  end subroutine actual_network_finalize
  
  subroutine ener_gener_rate(dydt, enuc)
    ! Computes the instantaneous energy generation rate
    !$acc routine seq
  
    implicit none

    double precision :: dydt(nspec), enuc

    ! This is basically e = m c**2

    enuc = sum(dydt(:) * mion(:)) * enuc_conv2

  end subroutine ener_gener_rate

end module actual_network
