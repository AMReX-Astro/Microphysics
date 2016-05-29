module actual_network

  use bl_types
  use actual_network_data

  implicit none

  double precision, parameter :: MeV2erg = 1.60217646e-6
  double precision, parameter :: N_A = 6.0221415e23
  
contains

  subroutine actual_network_init

    use bl_constants_module

    ! set the names
    spec_names(ic12)  = "carbon-12"
    spec_names(io14)  = "oxygen-14"
    spec_names(io15)  = "oxygen-15"
    spec_names(io16)  = "oxygen-16"
    spec_names(if17)  = "flourine-17"
    spec_names(img22) = "magnesium-22"
    spec_names(is30)  = "sulfur-30"
    spec_names(ini56) = "nickel-56"
    spec_names(ihe4)  = "helium-4"
    spec_names(ih1)   = "hydrogen-1"

    short_spec_names(ic12)  = "C12"
    short_spec_names(io14)  = "O14"
    short_spec_names(io15)  = "O15"
    short_spec_names(io16)  = "O16"
    short_spec_names(if17)  = "F17"
    short_spec_names(img22) = "Mg22"
    short_spec_names(is30)  = "S30"
    short_spec_names(ini56) = "Ni56"
    short_spec_names(ihe4)  = "He4"
    short_spec_names(ih1)   = "H1"    

    ! set the species properties
    aion(ic12)  = TWELVE
    aion(io14)  = 14.0_dp_t
    aion(io15)  = 15.0_dp_t
    aion(io16)  = 16.0_dp_t
    aion(if17)  = 17.0_dp_t
    aion(img22) = 22.0_dp_t
    aion(is30)  = 30.0_dp_t
    aion(ini56) = 56.0_dp_t
    aion(ihe4)  = FOUR
    aion(ih1)   = ONE

    zion(ic12)  = SIX
    zion(io14)  = EIGHT
    zion(io15)  = EIGHT
    zion(io16)  = EIGHT
    zion(if17)  = NINE
    zion(img22) = TWELVE
    zion(is30)  = 16.0_dp_t
    zion(ini56) = 28.0_dp_t
    zion(ihe4)  = TWO
    zion(ih1)   = ONE

    ! Our convention is that binding energy is negative.  The
    ! following are the binding energies in MeV.
    ebin(ic12)  = 92.16279_dp_t
    ebin(io14)  = 98.7325_dp_t
    ebin(io15)  = 111.9569_dp_t
    ebin(io16)  = 127.6207_dp_t
    ebin(if17)  = 128.2211_dp_t
    ebin(img22) = 168.5768_dp_t
    ebin(is30)  = 243.6866_dp_t
    ebin(ini56) = 483.995_dp_t
    ebin(ihe4)  = 28.29599_dp_t
    ebin(ih1)   = ZERO

    ! convert to erg / g by multiplying by N_A / aion and converting to erg
    ebin = -ebin * N_A * MeV2erg / aion

  end subroutine actual_network_init



  subroutine actual_network_finalize

    implicit none

    ! Nothing to do here.

  end subroutine actual_network_finalize

end module actual_network
