module actual_network

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

  subroutine actual_network_init

    use rpar_indices

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

  end subroutine actual_network_init

  subroutine network_finalize()
    ! stub for Maestro
  end subroutine network_finalize

end module actual_network
