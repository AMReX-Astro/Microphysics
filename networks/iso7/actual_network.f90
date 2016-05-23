module actual_network

  use bl_types
  use actual_network_data

  implicit none

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

  end subroutine actual_network_init

end module actual_network
