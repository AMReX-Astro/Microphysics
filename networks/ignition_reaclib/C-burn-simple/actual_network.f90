module actual_network
  use physical_constants, only: ERG_PER_MeV
  use bl_types
  use actual_network_data
  
  implicit none

  double precision, parameter, private :: avo = 6.0221417930d23
  double precision, parameter, private :: c_light = 2.99792458d10

  double precision, parameter, private :: ev2erg  = 1.60217648740d-12
  double precision, parameter, private :: mev2erg = ev2erg*1.0d6
  double precision, parameter, private :: mev2gr  = mev2erg/c_light**2

  double precision, parameter, private :: mn = 1.67492721184d-24
  double precision, parameter, private :: mp = 1.67262163783d-24
  double precision, parameter, private :: me = 9.1093821545d-28

contains

  subroutine actual_network_init()

    implicit none
    
    integer :: i

    spec_names(jn)   = "neutron"
    spec_names(jp)   = "hydrogen"
    spec_names(jhe4)   = "helium-4"
    spec_names(jc12)   = "carbon-12"
    spec_names(jne20)   = "neon-20"
    spec_names(jna23)   = "sodium-23"
    spec_names(jmg23)   = "magnesium-23"

    short_spec_names(jn)   = "n"
    short_spec_names(jp)   = "p"
    short_spec_names(jhe4)   = "he4"
    short_spec_names(jc12)   = "c12"
    short_spec_names(jne20)   = "ne20"
    short_spec_names(jna23)   = "na23"
    short_spec_names(jmg23)   = "mg23"

    ebind_per_nucleon(jn)   = 0.00000000000000d+00
    ebind_per_nucleon(jp)   = 0.00000000000000d+00
    ebind_per_nucleon(jhe4)   = 7.07391500000000d+00
    ebind_per_nucleon(jc12)   = 7.68014400000000d+00
    ebind_per_nucleon(jne20)   = 8.03224000000000d+00
    ebind_per_nucleon(jna23)   = 8.11149300000000d+00
    ebind_per_nucleon(jmg23)   = 7.90110400000000d+00

    aion(jn)   = 1.00000000000000d+00
    aion(jp)   = 1.00000000000000d+00
    aion(jhe4)   = 4.00000000000000d+00
    aion(jc12)   = 1.20000000000000d+01
    aion(jne20)   = 2.00000000000000d+01
    aion(jna23)   = 2.30000000000000d+01
    aion(jmg23)   = 2.30000000000000d+01

    zion(jn)   = 0.00000000000000d+00
    zion(jp)   = 1.00000000000000d+00
    zion(jhe4)   = 2.00000000000000d+00
    zion(jc12)   = 6.00000000000000d+00
    zion(jne20)   = 1.00000000000000d+01
    zion(jna23)   = 1.10000000000000d+01
    zion(jmg23)   = 1.20000000000000d+01

    nion(jn)   = 1.00000000000000d+00
    nion(jp)   = 0.00000000000000d+00
    nion(jhe4)   = 2.00000000000000d+00
    nion(jc12)   = 6.00000000000000d+00
    nion(jne20)   = 1.00000000000000d+01
    nion(jna23)   = 1.20000000000000d+01
    nion(jmg23)   = 1.10000000000000d+01

    do i = 1, nspec
       bion(i) = ebind_per_nucleon(i) * aion(i) * ERG_PER_MeV
    end do

    ! Set the mass
    mion(:) = nion(:) * mn + zion(:) * (mp + me) &
         - bion(:)/(c_light**2)

    ! Molar mass
    wion(:) = avo * mion(:)

    ! Common approximation
    !wion(:) = aion(:)

    !$acc update device(aion, zion, bion, nion, mion, wion)
  end subroutine actual_network_init

  subroutine network_finalize()
    ! stub for Maestro
  end subroutine network_finalize
  
end module actual_network
