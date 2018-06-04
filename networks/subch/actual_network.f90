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
  double precision :: ebind_per_nucleon(nspec)

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

  double precision :: aion(nspec), zion(nspec), bion(nspec)
  double precision :: nion(nspec), mion(nspec), wion(nspec)

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

    ebind_per_nucleon(jp)   = 0.00000000000000d+00
    ebind_per_nucleon(jhe4)   = 7.07391500000000d+00
    ebind_per_nucleon(jc12)   = 7.68014400000000d+00
    ebind_per_nucleon(jc14)   = 7.52031900000000d+00
    ebind_per_nucleon(jn13)   = 7.23886300000000d+00
    ebind_per_nucleon(jn14)   = 7.47561400000000d+00
    ebind_per_nucleon(jo16)   = 7.97620600000000d+00
    ebind_per_nucleon(jo18)   = 7.76709700000000d+00
    ebind_per_nucleon(jf18)   = 7.63163800000000d+00
    ebind_per_nucleon(jne20)   = 8.03224000000000d+00
    ebind_per_nucleon(jne21)   = 7.97171300000000d+00

    aion(jp)   = 1.00000000000000d+00
    aion(jhe4)   = 4.00000000000000d+00
    aion(jc12)   = 1.20000000000000d+01
    aion(jc14)   = 1.40000000000000d+01
    aion(jn13)   = 1.30000000000000d+01
    aion(jn14)   = 1.40000000000000d+01
    aion(jo16)   = 1.60000000000000d+01
    aion(jo18)   = 1.80000000000000d+01
    aion(jf18)   = 1.80000000000000d+01
    aion(jne20)   = 2.00000000000000d+01
    aion(jne21)   = 2.10000000000000d+01

    zion(jp)   = 1.00000000000000d+00
    zion(jhe4)   = 2.00000000000000d+00
    zion(jc12)   = 6.00000000000000d+00
    zion(jc14)   = 6.00000000000000d+00
    zion(jn13)   = 7.00000000000000d+00
    zion(jn14)   = 7.00000000000000d+00
    zion(jo16)   = 8.00000000000000d+00
    zion(jo18)   = 8.00000000000000d+00
    zion(jf18)   = 9.00000000000000d+00
    zion(jne20)   = 1.00000000000000d+01
    zion(jne21)   = 1.00000000000000d+01

    nion(jp)   = 0.00000000000000d+00
    nion(jhe4)   = 2.00000000000000d+00
    nion(jc12)   = 6.00000000000000d+00
    nion(jc14)   = 8.00000000000000d+00
    nion(jn13)   = 6.00000000000000d+00
    nion(jn14)   = 7.00000000000000d+00
    nion(jo16)   = 8.00000000000000d+00
    nion(jo18)   = 1.00000000000000d+01
    nion(jf18)   = 9.00000000000000d+00
    nion(jne20)   = 1.00000000000000d+01
    nion(jne21)   = 1.10000000000000d+01

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

end module actual_network
