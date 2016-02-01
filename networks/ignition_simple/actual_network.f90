module actual_network

  use bl_types
  use actual_network_data

  implicit none

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

    aion(ic12)  = 12.0_dp_t
    aion(io16)  = 16.0_dp_t
    aion(img24) = 24.0_dp_t
    
    zion(ic12)  = 6.0_dp_t
    zion(io16)  = 8.0_dp_t
    zion(img24) = 12.0_dp_t

    ! our convention is that the binding energies are negative.  We convert
    ! from the MeV values that are traditionally written in astrophysics 
    ! papers by multiplying by 1.e6 eV/MeV * 1.60217646e-12 erg/eV.  The
    ! MeV values are per nucleus, so we divide by aion to make it per
    ! nucleon and we multiple by Avogardo's # (6.0221415e23) to get the 
    ! value in erg/g
    ebin(ic12)  = -7.4103097e18_dp_t     !  92.16294 MeV
    ebin(io16)  = -7.6959672e18_dp_t     ! 127.62093 MeV
    ebin(img24) = -7.9704080e18_dp_t     ! 198.2579  MeV

  end subroutine actual_network_init

  subroutine network_finalize()
    ! stub for Maestro
  end subroutine network_finalize

end module actual_network
