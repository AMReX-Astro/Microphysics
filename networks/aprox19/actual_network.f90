module actual_network

  use bl_types
  use actual_network_data

  implicit none

contains
  
  subroutine actual_network_init

    implicit none
    
    ! The following comes directly from init_aprox19

    short_spec_names(ih1)   = 'h1'
    short_spec_names(ihe3)  = 'he3'  
    short_spec_names(ihe4)  = 'he4'
    short_spec_names(ic12)  = 'c12'
    short_spec_names(in14)  = 'n14'
    short_spec_names(io16)  = 'o16'
    short_spec_names(ine20) = 'ne20'
    short_spec_names(img24) = 'mg24'
    short_spec_names(isi28) = 'si28'
    short_spec_names(is32)  = 's32'
    short_spec_names(iar36) = 'ar36'
    short_spec_names(ica40) = 'ca40'
    short_spec_names(iti44) = 'ti44'
    short_spec_names(icr48) = 'cr48'
    short_spec_names(ife52) = 'fe52'
    short_spec_names(ife54) = 'fe54'
    short_spec_names(ini56) = 'ni56'
    short_spec_names(ineut) = 'neut'
    short_spec_names(iprot) = 'prot'

    spec_names(ih1)   = "hydrogen-1"
    spec_names(ihe3)  = "helium-3"
    spec_names(ihe4)  = "helium-4"
    spec_names(ic12)  = "carbon-12"
    spec_names(in14)  = "nitrogen-14"
    spec_names(io16)  = "oxygen-16"
    spec_names(ine20) = "neon-20"
    spec_names(img24) = "magnesium-24"
    spec_names(isi28) = "silicon-28"
    spec_names(is32)  = "sulfur-32"
    spec_names(iar36) = "argon-36"
    spec_names(ica40) = "calcium-40"
    spec_names(iti44) = "titanium-44"
    spec_names(icr48) = "chromium-48"
    spec_names(ife52) = "iron-52"
    spec_names(ife54) = "iron-54"
    spec_names(ini56) = "nickel-56"
    spec_names(ineut) = "neutron"
    spec_names(iprot) = "proton"
    

    ! Set the number of nucleons in the element
    aion(ih1)   = 1.0d0
    aion(ihe3)  = 3.0d0
    aion(ihe4)  = 4.0d0
    aion(ic12)  = 12.0d0
    aion(in14)  = 14.0d0
    aion(io16)  = 16.0d0
    aion(ine20) = 20.0d0
    aion(img24) = 24.0d0
    aion(isi28) = 28.0d0
    aion(is32)  = 32.0d0
    aion(iar36) = 36.0d0
    aion(ica40) = 40.0d0
    aion(iti44) = 44.0d0
    aion(icr48) = 48.0d0
    aion(ife52) = 52.0d0
    aion(ife54) = 54.0d0
    aion(ini56) = 56.0d0
    aion(ineut) = 1.0d0
    aion(iprot) = 1.0d0

    ! Set the number of protons in the element
    zion(ih1)   = 1.0d0
    zion(ihe3)  = 2.0d0
    zion(ihe4)  = 2.0d0
    zion(ic12)  = 6.0d0
    zion(in14)  = 7.0d0
    zion(io16)  = 8.0d0
    zion(ine20) = 10.0d0
    zion(img24) = 12.0d0
    zion(isi28) = 14.0d0
    zion(is32)  = 16.0d0
    zion(iar36) = 18.0d0
    zion(ica40) = 20.0d0
    zion(iti44) = 22.0d0
    zion(icr48) = 24.0d0
    zion(ife52) = 26.0d0
    zion(ife54) = 26.0d0
    zion(ini56) = 28.0d0
    zion(ineut) = 0.0d0
    zion(iprot) = 1.0d0    

    ! Set the binding energy of the element
    bion(ih1)   = 0.0d0
    bion(ihe3)  = 7.71819d0
    bion(ihe4)  = 28.29603d0
    bion(ic12)  = 92.16294d0
    bion(in14)  = 104.65998d0
    bion(io16)  = 127.62093d0
    bion(ine20) = 160.64788d0
    bion(img24) = 198.25790d0
    bion(isi28) = 236.53790d0
    bion(is32)  = 271.78250d0
    bion(iar36) = 306.72020d0
    bion(ica40) = 342.05680d0
    bion(iti44) = 375.47720d0
    bion(icr48) = 411.46900d0
    bion(ife52) = 447.70800d0
    bion(ife54) = 471.7696d0
    bion(ini56) = 484.00300d0
    bion(ineut) = 0.0d0
    bion(iprot) = 0.0d0

    ! Set the number of neutrons
    nion(:) = aion(:) - zion(:)

    ! Set the mass
    mion(:) = nion(:) * mn + zion(:) * (mp + me) - bion(:) * mev2gr

    ! Molar mass
    wion(:) = avo * mion(:)

    ! Common approximation
    wion(:) = aion(:)

  end subroutine actual_network_init



  subroutine actual_network_finalize

    implicit none

    ! Nothing to do here.

  end subroutine actual_network_finalize

end module actual_network
