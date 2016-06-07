module actual_network

  use bl_types

  implicit none

  double precision, parameter, private :: clight = 2.99792458d10
  double precision, parameter, private :: ev2erg  = 1.60217648740d-12
  double precision, parameter, private :: mev2erg = ev2erg*1.0d6
  double precision, parameter, private :: mev2gr  = mev2erg/clight**2

  double precision, parameter, private :: mn = 1.67492721184d-24
  double precision, parameter, private :: mp = 1.67262163783d-24
  double precision, parameter, private :: me = 9.1093821545d-28

  integer, parameter :: nspec = 3
  integer, parameter :: nspec_evolve = 1
  integer, parameter :: naux  = 0

  integer, parameter :: ic12  = 1
  integer, parameter :: io16  = 2
  integer, parameter :: img24 = 3

  character (len=16), save :: spec_names(nspec) 
  character (len= 5), save :: short_spec_names(nspec)
  character (len= 5), save :: short_aux_names(naux)

  double precision :: aion(nspec), zion(nspec), bion(nspec)
  double precision :: nion(nspec), mion(nspec), wion(nspec)

  !$acc declare create(aion, zion, bion, nion, mion, wion)

  integer, parameter :: nrates = 1

  ! Conversion factor for the nuclear energy generation rate.

  double precision, parameter :: avo = 6.0221417930d23
  double precision, parameter :: enuc_conv2 = -avo*clight*clight

contains

  subroutine actual_network_init

    implicit none

    spec_names(ic12)  = "carbon-12"
    spec_names(io16)  = "oxygen-16"
    spec_names(img24) = "magnesium-24"

    short_spec_names(ic12)  = "C12"
    short_spec_names(io16)  = "O16"
    short_spec_names(img24) = "Mg24"

    aion(ic12)  = 12.0d0
    aion(io16)  = 16.0d0
    aion(img24) = 24.0d0

    zion(ic12)  = 6.0d0
    zion(io16)  = 8.0d0
    zion(img24) = 12.0d0

    ! Binding energies per nucleon in MeV
    bion(ic12)  = 92.16294d0
    bion(io16)  = 127.62093d0
    bion(img24) = 198.2579d0

    ! Set the number of neutrons
    nion(:) = aion(:) - zion(:)

    ! Set the mass
    mion(:) = nion(:) * mn + zion(:) * (mp + me) - bion(:) * mev2gr

    ! Molar mass
    wion(:) = avo * mion(:)

    ! Common approximation
    wion(:) = aion(:)

    !$acc update device(aion, zion, bion, nion, mion, wion)

  end subroutine actual_network_init



  subroutine actual_network_finalize

    implicit none

    ! Nothing to do here.

  end subroutine actual_network_finalize

end module actual_network
