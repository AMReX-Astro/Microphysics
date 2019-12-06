module actual_network

  use amrex_fort_module, only : rt => amrex_real
  use physical_constants, only: ERG_PER_MeV
  use microphysics_type_module
  
  implicit none

  public num_rate_groups

  real(rt), parameter :: avo = 6.0221417930e23_rt
  real(rt), parameter :: c_light = 2.99792458e10_rt
  real(rt), parameter :: enuc_conv2 = -avo*c_light*c_light

  real(rt), parameter :: ev2erg  = 1.60217648740e-12_rt
  real(rt), parameter :: mev2erg = ev2erg*1.0e6_rt
  real(rt), parameter :: mev2gr  = mev2erg/c_light**2

  real(rt), parameter :: mass_neutron  = 1.67492721184e-24_rt
  real(rt), parameter :: mass_proton   = 1.67262163783e-24_rt
  real(rt), parameter :: mass_electron = 9.10938215450e-28_rt

  integer, parameter :: nrates = 8
  integer, parameter :: num_rate_groups = 4

  ! Evolution and auxiliary
  integer, parameter :: nspec_evolve = 11
  integer, parameter :: naux  = 0

  ! Number of nuclear species in the network
  integer, parameter :: nspec = 11

  ! Number of reaclib rates
  integer, parameter :: nrat_reaclib = 8
  
  ! Number of tabular rates
  integer, parameter :: nrat_tabular = 0

  ! Binding Energies Per Nucleon (MeV)
  real(rt) :: ebind_per_nucleon(nspec)

  ! aion: Nucleon mass number A
  ! zion: Nucleon atomic number Z
  ! nion: Nucleon neutron number N
  ! bion: Binding Energies (ergs)

  ! Nuclides
  integer, parameter :: jp   = 1
  integer, parameter :: jhe4   = 2
  integer, parameter :: jc12   = 3
  integer, parameter :: jc14   = 4
  integer, parameter :: jn13   = 5
  integer, parameter :: jn14   = 6
  integer, parameter :: jo16   = 7
  integer, parameter :: jo18   = 8
  integer, parameter :: jf18   = 9
  integer, parameter :: jne20   = 10
  integer, parameter :: jne21   = 11

  ! Reactions
  integer, parameter :: k_he4_he4_he4__c12   = 1
  integer, parameter :: k_he4_c12__o16   = 2
  integer, parameter :: k_he4_n14__f18   = 3
  integer, parameter :: k_he4_f18__p_ne21   = 4
  integer, parameter :: k_p_c12__n13   = 5
  integer, parameter :: k_he4_n13__p_o16   = 6
  integer, parameter :: k_he4_o16__ne20   = 7
  integer, parameter :: k_he4_c14__o18   = 8

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

  real(rt), allocatable :: aion(:), zion(:), bion(:)
  real(rt), allocatable :: nion(:), mion(:), wion(:)

#ifdef AMREX_USE_CUDA
  attributes(managed) :: aion, zion, bion, nion, mion, wion
#endif

  !$acc declare create(aion, zion, bion, nion, mion, wion)

contains

  subroutine actual_network_init()
    
    implicit none
    
    integer :: i

    spec_names(jp)   = "hydrogen-1"
    spec_names(jhe4)   = "helium-4"
    spec_names(jc12)   = "carbon-12"
    spec_names(jc14)   = "carbon-14"
    spec_names(jn13)   = "nitrogen-13"
    spec_names(jn14)   = "nitrogen-14"
    spec_names(jo16)   = "oxygen-16"
    spec_names(jo18)   = "oxygen-18"
    spec_names(jf18)   = "fluorine-18"
    spec_names(jne20)   = "neon-20"
    spec_names(jne21)   = "neon-21"

    short_spec_names(jp)   = "h1"
    short_spec_names(jhe4)   = "he4"
    short_spec_names(jc12)   = "c12"
    short_spec_names(jc14)   = "c14"
    short_spec_names(jn13)   = "n13"
    short_spec_names(jn14)   = "n14"
    short_spec_names(jo16)   = "o16"
    short_spec_names(jo18)   = "o18"
    short_spec_names(jf18)   = "f18"
    short_spec_names(jne20)   = "ne20"
    short_spec_names(jne21)   = "ne21"

    allocate(aion(nspec))
    allocate(zion(nspec))
    allocate(nion(nspec))
    allocate(bion(nspec))
    allocate(mion(nspec))
    allocate(wion(nspec))
    
    ebind_per_nucleon(jp)   = 0.00000000000000e+00_rt
    ebind_per_nucleon(jhe4)   = 7.07391500000000e+00_rt
    ebind_per_nucleon(jc12)   = 7.68014400000000e+00_rt
    ebind_per_nucleon(jc14)   = 7.52031900000000e+00_rt
    ebind_per_nucleon(jn13)   = 7.23886300000000e+00_rt
    ebind_per_nucleon(jn14)   = 7.47561400000000e+00_rt
    ebind_per_nucleon(jo16)   = 7.97620600000000e+00_rt
    ebind_per_nucleon(jo18)   = 7.76709700000000e+00_rt
    ebind_per_nucleon(jf18)   = 7.63163800000000e+00_rt
    ebind_per_nucleon(jne20)   = 8.03224000000000e+00_rt
    ebind_per_nucleon(jne21)   = 7.97171300000000e+00_rt

    aion(jp)   = 1.00000000000000e+00_rt
    aion(jhe4)   = 4.00000000000000e+00_rt
    aion(jc12)   = 1.20000000000000e+01_rt
    aion(jc14)   = 1.40000000000000e+01_rt
    aion(jn13)   = 1.30000000000000e+01_rt
    aion(jn14)   = 1.40000000000000e+01_rt
    aion(jo16)   = 1.60000000000000e+01_rt
    aion(jo18)   = 1.80000000000000e+01_rt
    aion(jf18)   = 1.80000000000000e+01_rt
    aion(jne20)   = 2.00000000000000e+01_rt
    aion(jne21)   = 2.10000000000000e+01_rt

    zion(jp)   = 1.00000000000000e+00_rt
    zion(jhe4)   = 2.00000000000000e+00_rt
    zion(jc12)   = 6.00000000000000e+00_rt
    zion(jc14)   = 6.00000000000000e+00_rt
    zion(jn13)   = 7.00000000000000e+00_rt
    zion(jn14)   = 7.00000000000000e+00_rt
    zion(jo16)   = 8.00000000000000e+00_rt
    zion(jo18)   = 8.00000000000000e+00_rt
    zion(jf18)   = 9.00000000000000e+00_rt
    zion(jne20)   = 1.00000000000000e+01_rt
    zion(jne21)   = 1.00000000000000e+01_rt

    nion(jp)   = 0.00000000000000e+00_rt
    nion(jhe4)   = 2.00000000000000e+00_rt
    nion(jc12)   = 6.00000000000000e+00_rt
    nion(jc14)   = 8.00000000000000e+00_rt
    nion(jn13)   = 6.00000000000000e+00_rt
    nion(jn14)   = 7.00000000000000e+00_rt
    nion(jo16)   = 8.00000000000000e+00_rt
    nion(jo18)   = 1.00000000000000e+01_rt
    nion(jf18)   = 9.00000000000000e+00_rt
    nion(jne20)   = 1.00000000000000e+01_rt
    nion(jne21)   = 1.10000000000000e+01_rt

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
