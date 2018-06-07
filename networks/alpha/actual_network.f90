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

  integer, parameter :: nrates = 34
  integer, parameter :: num_rate_groups = 4

  ! Evolution and auxiliary
  integer, parameter :: nspec_evolve = 14
  integer, parameter :: naux  = 0

  ! Number of nuclear species in the network
  integer, parameter :: nspec = 14

  ! Number of reaclib rates
  integer, parameter :: nrat_reaclib = 34
  
  ! Number of tabular rates
  integer, parameter :: nrat_tabular = 0

  ! Binding Energies Per Nucleon (MeV)
  double precision :: ebind_per_nucleon(nspec)

  ! aion: Nucleon mass number A
  ! zion: Nucleon atomic number Z
  ! nion: Nucleon neutron number N
  ! bion: Binding Energies (ergs)

  ! Nuclides
  integer, parameter :: jhe4   = 1
  integer, parameter :: jc12   = 2
  integer, parameter :: jo16   = 3
  integer, parameter :: jne20   = 4
  integer, parameter :: jmg24   = 5
  integer, parameter :: jsi28   = 6
  integer, parameter :: js32   = 7
  integer, parameter :: jar36   = 8
  integer, parameter :: jca40   = 9
  integer, parameter :: jti44   = 10
  integer, parameter :: jcr48   = 11
  integer, parameter :: jfe52   = 12
  integer, parameter :: jni56   = 13
  integer, parameter :: jzn60   = 14

  ! Reactions
  integer, parameter :: k_o16__he4_c12   = 1
  integer, parameter :: k_ne20__he4_o16   = 2
  integer, parameter :: k_mg24__he4_ne20   = 3
  integer, parameter :: k_si28__he4_mg24   = 4
  integer, parameter :: k_s32__he4_si28   = 5
  integer, parameter :: k_ar36__he4_s32   = 6
  integer, parameter :: k_ca40__he4_ar36   = 7
  integer, parameter :: k_ti44__he4_ca40   = 8
  integer, parameter :: k_cr48__he4_ti44   = 9
  integer, parameter :: k_fe52__he4_cr48   = 10
  integer, parameter :: k_ni56__he4_fe52   = 11
  integer, parameter :: k_zn60__he4_ni56   = 12
  integer, parameter :: k_c12__he4_he4_he4   = 13
  integer, parameter :: k_he4_c12__o16   = 14
  integer, parameter :: k_he4_o16__ne20   = 15
  integer, parameter :: k_he4_ne20__mg24   = 16
  integer, parameter :: k_he4_mg24__si28   = 17
  integer, parameter :: k_he4_si28__s32   = 18
  integer, parameter :: k_he4_s32__ar36   = 19
  integer, parameter :: k_he4_ar36__ca40   = 20
  integer, parameter :: k_he4_ca40__ti44   = 21
  integer, parameter :: k_he4_ti44__cr48   = 22
  integer, parameter :: k_he4_cr48__fe52   = 23
  integer, parameter :: k_he4_fe52__ni56   = 24
  integer, parameter :: k_he4_ni56__zn60   = 25
  integer, parameter :: k_c12_c12__he4_ne20   = 26
  integer, parameter :: k_c12_o16__he4_mg24   = 27
  integer, parameter :: k_o16_o16__he4_si28   = 28
  integer, parameter :: k_he4_ne20__c12_c12   = 29
  integer, parameter :: k_c12_ne20__he4_si28   = 30
  integer, parameter :: k_he4_mg24__c12_o16   = 31
  integer, parameter :: k_he4_si28__c12_ne20   = 32
  integer, parameter :: k_he4_si28__o16_o16   = 33
  integer, parameter :: k_he4_he4_he4__c12   = 34

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

    spec_names(jhe4)   = "helium-4"
    spec_names(jc12)   = "carbon-12"
    spec_names(jo16)   = "oxygen-16"
    spec_names(jne20)   = "neon-20"
    spec_names(jmg24)   = "magnesium-24"
    spec_names(jsi28)   = "silicon-28"
    spec_names(js32)   = "sulfur-32"
    spec_names(jar36)   = "argon-36"
    spec_names(jca40)   = "calcium-40"
    spec_names(jti44)   = "titanium-44"
    spec_names(jcr48)   = "chromium-48"
    spec_names(jfe52)   = "iron-52"
    spec_names(jni56)   = "nickel-56"
    spec_names(jzn60)   = "zinc-60"

    short_spec_names(jhe4)   = "he4"
    short_spec_names(jc12)   = "c12"
    short_spec_names(jo16)   = "o16"
    short_spec_names(jne20)   = "ne20"
    short_spec_names(jmg24)   = "mg24"
    short_spec_names(jsi28)   = "si28"
    short_spec_names(js32)   = "s32"
    short_spec_names(jar36)   = "ar36"
    short_spec_names(jca40)   = "ca40"
    short_spec_names(jti44)   = "ti44"
    short_spec_names(jcr48)   = "cr48"
    short_spec_names(jfe52)   = "fe52"
    short_spec_names(jni56)   = "ni56"
    short_spec_names(jzn60)   = "zn60"

    ebind_per_nucleon(jhe4)   = 7.07391500000000d+00
    ebind_per_nucleon(jc12)   = 7.68014400000000d+00
    ebind_per_nucleon(jo16)   = 7.97620600000000d+00
    ebind_per_nucleon(jne20)   = 8.03224000000000d+00
    ebind_per_nucleon(jmg24)   = 8.26070900000000d+00
    ebind_per_nucleon(jsi28)   = 8.44774400000000d+00
    ebind_per_nucleon(js32)   = 8.49312900000000d+00
    ebind_per_nucleon(jar36)   = 8.51990900000000d+00
    ebind_per_nucleon(jca40)   = 8.55130300000000d+00
    ebind_per_nucleon(jti44)   = 8.53352000000000d+00
    ebind_per_nucleon(jcr48)   = 8.57226900000000d+00
    ebind_per_nucleon(jfe52)   = 8.60957400000000d+00
    ebind_per_nucleon(jni56)   = 8.64277900000000d+00
    ebind_per_nucleon(jzn60)   = 8.58305000000000d+00

    aion(jhe4)   = 4.00000000000000d+00
    aion(jc12)   = 1.20000000000000d+01
    aion(jo16)   = 1.60000000000000d+01
    aion(jne20)   = 2.00000000000000d+01
    aion(jmg24)   = 2.40000000000000d+01
    aion(jsi28)   = 2.80000000000000d+01
    aion(js32)   = 3.20000000000000d+01
    aion(jar36)   = 3.60000000000000d+01
    aion(jca40)   = 4.00000000000000d+01
    aion(jti44)   = 4.40000000000000d+01
    aion(jcr48)   = 4.80000000000000d+01
    aion(jfe52)   = 5.20000000000000d+01
    aion(jni56)   = 5.60000000000000d+01
    aion(jzn60)   = 6.00000000000000d+01

    zion(jhe4)   = 2.00000000000000d+00
    zion(jc12)   = 6.00000000000000d+00
    zion(jo16)   = 8.00000000000000d+00
    zion(jne20)   = 1.00000000000000d+01
    zion(jmg24)   = 1.20000000000000d+01
    zion(jsi28)   = 1.40000000000000d+01
    zion(js32)   = 1.60000000000000d+01
    zion(jar36)   = 1.80000000000000d+01
    zion(jca40)   = 2.00000000000000d+01
    zion(jti44)   = 2.20000000000000d+01
    zion(jcr48)   = 2.40000000000000d+01
    zion(jfe52)   = 2.60000000000000d+01
    zion(jni56)   = 2.80000000000000d+01
    zion(jzn60)   = 3.00000000000000d+01

    nion(jhe4)   = 2.00000000000000d+00
    nion(jc12)   = 6.00000000000000d+00
    nion(jo16)   = 8.00000000000000d+00
    nion(jne20)   = 1.00000000000000d+01
    nion(jmg24)   = 1.20000000000000d+01
    nion(jsi28)   = 1.40000000000000d+01
    nion(js32)   = 1.60000000000000d+01
    nion(jar36)   = 1.80000000000000d+01
    nion(jca40)   = 2.00000000000000d+01
    nion(jti44)   = 2.20000000000000d+01
    nion(jcr48)   = 2.40000000000000d+01
    nion(jfe52)   = 2.60000000000000d+01
    nion(jni56)   = 2.80000000000000d+01
    nion(jzn60)   = 3.00000000000000d+01

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
