module actual_network

  use amrex_fort_module, only : rt => amrex_real

  implicit none

  integer, parameter :: nspec = 7
  integer, parameter :: nspec_evolve = 7
  integer, parameter :: naux  = 0
  
  integer, parameter :: ihe4  = 1
  integer, parameter :: ic12  = 2
  integer, parameter :: io16  = 3
  integer, parameter :: ine20 = 4
  integer, parameter :: img24 = 5
  integer, parameter :: isi28 = 6
  integer, parameter :: ini56 = 7

  double precision, save :: aion(nspec), zion(nspec), nion(nspec)
  double precision, save :: bion(nspec), mion(nspec), wion(nspec)

  character (len=16), save :: spec_names(nspec)
  character (len= 5), save :: short_spec_names(nspec)
  character (len= 5), save :: short_aux_names(naux)

  character (len=32), parameter :: network_name = "iso7"

  ! Some fundamental physical constants

  double precision, parameter :: avo = 6.0221417930d23
  double precision, parameter :: c_light = 2.99792458d10

  double precision, parameter :: ev2erg  = 1.60217648740d-12
  double precision, parameter :: mev2erg = ev2erg*1.0d6
  double precision, parameter :: mev2gr  = mev2erg/c_light**2

  double precision, parameter :: mn = 1.67492721184d-24
  double precision, parameter :: mp = 1.67262163783d-24
  double precision, parameter :: me = 9.1093821545d-28

  ! Conversion factor for the nuclear energy generation rate.

  double precision, parameter :: enuc_conv2 = -avo*c_light*c_light

  ! Rates data

  integer, parameter :: nrates  = 17
  integer, parameter :: num_rate_groups = 4

  integer, parameter :: ircag   = 1
  integer, parameter :: iroga   = 2
  integer, parameter :: ir3a    = 3
  integer, parameter :: irg3a   = 4
  integer, parameter :: ir1212  = 5
  integer, parameter :: ir1216  = 6
  integer, parameter :: ir1616  = 7
  integer, parameter :: iroag   = 8
  integer, parameter :: irnega  = 9
  integer, parameter :: irneag  = 10
  integer, parameter :: irmgga  = 11
  integer, parameter :: irmgag  = 12
  integer, parameter :: irsiga  = 13
  integer, parameter :: ircaag  = 14
  integer, parameter :: irtiga  = 15
  integer, parameter :: irsi2ni = 16
  integer, parameter :: irni2si = 17

  character (len=20), save :: ratenames(nrates)

contains

  subroutine actual_network_init

    implicit none

    short_spec_names(ihe4)  = 'he4'
    short_spec_names(ic12)  = 'c12'
    short_spec_names(io16)  = 'o16'
    short_spec_names(ine20) = 'ne20'
    short_spec_names(img24) = 'mg24'
    short_spec_names(isi28) = 'si28'
    short_spec_names(ini56) = 'ni56'

    spec_names(ihe4)  = "helium-4"
    spec_names(ic12)  = "carbon-12"
    spec_names(io16)  = "oxygen-16"
    spec_names(ine20) = "neon-20"
    spec_names(img24) = "magnesium-24"
    spec_names(isi28) = "silicon-28"
    spec_names(ini56) = "nickel-56"


    ! Set the number of nucleons in the element
    aion(ihe4)  = 4.0d0
    aion(ic12)  = 12.0d0
    aion(io16)  = 16.0d0
    aion(ine20) = 20.0d0
    aion(img24) = 24.0d0
    aion(isi28) = 28.0d0
    aion(ini56) = 56.0d0

    ! Set the number of protons in the element
    zion(ihe4)  = 2.0d0
    zion(ic12)  = 6.0d0
    zion(io16)  = 8.0d0
    zion(ine20) = 10.0d0
    zion(img24) = 12.0d0
    zion(isi28) = 14.0d0
    zion(ini56) = 28.0d0

    ! Set the binding energy of the element
    bion(ihe4)  = 28.29603d0
    bion(ic12)  = 92.16294d0
    bion(io16)  = 127.62093d0
    bion(ine20) = 160.64788d0
    bion(img24) = 198.25790d0
    bion(isi28) = 236.53790d0
    bion(ini56) = 484.00300d0

    ! Set the number of neutrons
    nion(:) = aion(:) - zion(:)

    ! Set the mass
    mion(:) = nion(:) * mn + zion(:) * (mp + me) - bion(:) * mev2gr

    ! Molar mass
    wion(:) = avo * mion(:)

    ! Common approximation
    wion(:) = aion(:)

    ! set the names of the reaction rates
    ratenames(ircag)   = 'rcag '
    ratenames(iroga)   = 'roga '
    ratenames(ir3a)    = 'r3a  '
    ratenames(irg3a)   = 'rg3a '    ! inverse rate
    ratenames(ir1212)  = 'r1212'
    ratenames(ir1216)  = 'r1216'
    ratenames(ir1616)  = 'r1616'
    ratenames(iroag)   = 'roag '
    ratenames(irnega)  = 'rnega'
    ratenames(irneag)  = 'rneag'
    ratenames(irmgga)  = 'rmgga'
    ratenames(irmgag)  = 'rmgag'
    ratenames(irsiga)  = 'rsiga'
    ratenames(ircaag)  = 'rcaag'
    ratenames(irtiga)  = 'rtiga'

    ratenames(irsi2ni) = 'rsi2ni'
    ratenames(irni2si) = 'rni2si'

  end subroutine actual_network_init



  subroutine actual_network_finalize

    implicit none

    ! Nothing to do here.

  end subroutine actual_network_finalize

end module actual_network
