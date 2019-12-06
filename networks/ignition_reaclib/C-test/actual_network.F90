module actual_network

  use physical_constants, only: ERG_PER_MeV
  use microphysics_type_module
  
  implicit none

  public

  real(rt), parameter :: avo = 6.0221417930e23_rt
  real(rt), parameter :: c_light = 2.99792458e10_rt
  real(rt), parameter :: enuc_conv2 = -avo*c_light*c_light

  real(rt), parameter :: ev2erg  = 1.60217648740e-12_rt
  real(rt), parameter :: mev2erg = ev2erg*1.0e6_rt
  real(rt), parameter :: mev2gr  = mev2erg/c_light**2

  real(rt), parameter :: mass_neutron  = 1.67492721184e-24_rt
  real(rt), parameter :: mass_proton   = 1.67262163783e-24_rt
  real(rt), parameter :: mass_electron = 9.10938215450e-28_rt

  integer, parameter :: nrates = 1
  integer, parameter :: num_rate_groups = 4

  ! Evolution and auxiliary
  integer, parameter :: nspec_evolve = 3
  integer, parameter :: naux  = 0

  ! Number of nuclear species in the network
  integer, parameter :: nspec = 3

  ! Number of reaclib rates
  integer, parameter :: nrat_reaclib = 1
  
  ! Number of tabular rates
  integer, parameter :: nrat_tabular = 0

  ! Binding Energies Per Nucleon (MeV)
  real(rt) :: ebind_per_nucleon(nspec)

  ! aion: Nucleon mass number A
  ! zion: Nucleon atomic number Z
  ! nion: Nucleon neutron number N
  ! bion: Binding Energies (ergs)

  ! Nuclides
  integer, parameter :: jhe4   = 1
  integer, parameter :: jc12   = 2
  integer, parameter :: jo16   = 3

  ! Reactions
  integer, parameter :: k_he4_c12__o16   = 1

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

  real(rt), allocatable, save :: aion(:), zion(:), bion(:)
  real(rt), allocatable, save :: nion(:), mion(:), wion(:)

#ifdef AMREX_USE_CUDA
  attributes(managed) :: aion, zion, bion, nion, mion, wion
#endif

  !$acc declare create(aion, zion, bion, nion, mion, wion)

contains

  subroutine actual_network_init()
    
    implicit none
    
    integer :: i

    ! Allocate ion info arrays
    allocate(aion(nspec))
    allocate(zion(nspec))
    allocate(bion(nspec))
    allocate(nion(nspec))
    allocate(mion(nspec))
    allocate(wion(nspec))

    spec_names(jhe4)   = "helium-4"
    spec_names(jc12)   = "carbon-12"
    spec_names(jo16)   = "oxygen-16"

    short_spec_names(jhe4)   = "he4"
    short_spec_names(jc12)   = "c12"
    short_spec_names(jo16)   = "o16"

    ebind_per_nucleon(jhe4)   = 7.07391500000000e+00_rt
    ebind_per_nucleon(jc12)   = 7.68014400000000e+00_rt
    ebind_per_nucleon(jo16)   = 7.97620600000000e+00_rt

    aion(jhe4)   = 4.00000000000000e+00_rt
    aion(jc12)   = 1.20000000000000e+01_rt
    aion(jo16)   = 1.60000000000000e+01_rt

    zion(jhe4)   = 2.00000000000000e+00_rt
    zion(jc12)   = 6.00000000000000e+00_rt
    zion(jo16)   = 8.00000000000000e+00_rt

    nion(jhe4)   = 2.00000000000000e+00_rt
    nion(jc12)   = 6.00000000000000e+00_rt
    nion(jo16)   = 8.00000000000000e+00_rt

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
    ! Deallocate storage arrays
    deallocate(aion)
    deallocate(zion)
    deallocate(bion)
    deallocate(nion)
    deallocate(mion)
    deallocate(wion)
  end subroutine actual_network_finalize


  subroutine ener_gener_rate(dydt, enuc)
    ! Computes the instantaneous energy generation rate

    !$acc routine seq
  
    implicit none

    real(rt) :: dydt(nspec), enuc

    !$gpu    

    ! This is basically e = m c**2

    enuc = sum(dydt(:) * mion(:)) * enuc_conv2

  end subroutine ener_gener_rate

end module actual_network
