module actual_network

  use microphysics_type_module, only: rt, ZERO, ONE, TWO, FOUR, EIGHT, TEN

  integer, parameter :: nspec = 7
  integer, parameter :: nspec_evolve = 7
  integer, parameter :: naux = 0
  integer, parameter :: nrates = 6  ! including constant weak rates
  integer, parameter :: num_rate_groups = 1

  real(rt), save :: aion(nspec), zion(nspec), bion(nspec)
  real(rt) :: rates(nrates)

  character (len=16), save :: spec_names(nspec)
  character (len=5), save :: short_spec_names(nspec)
  character (len=5), save :: short_aux_names(naux)

  character (len=32), parameter :: network_name = "xrb_simple"

  integer, parameter :: ih1 = 1
  integer, parameter :: ihe4 = 2
  integer, parameter :: io14 = 3
  integer, parameter :: io15 = 4
  integer, parameter :: ine18 = 5
  integer, parameter :: isi25 = 6
  integer, parameter :: ife56 = 7

  integer, parameter :: ir3a    = 1
  integer, parameter :: irag15  = 2
  integer, parameter :: irap14  = 3
  integer, parameter :: irap18  = 4
  integer, parameter :: irwk14o = 5
  integer, parameter :: irwk15o = 6

  ! hard-coded weak rates from Stan; should double check with newer estimates
  real(rt) :: wk14o = 9.832e-3_rt
  real(rt) :: wk15o = 5.682e-3_rt

contains

  subroutine actual_network_init()
    
    real(rt), parameter :: MeV2erg = 1.60217646e-6_rt, N_A = 6.0221415e23_rt

    short_spec_names(ih1) = "h1"
    short_spec_names(ihe4) = "he4"
    short_spec_names(io14) = "o14"
    short_spec_names(io15) = "o15"
    short_spec_names(ine18) = "ne18"
    short_spec_names(isi25) = "si25"
    short_spec_names(ife56) = "fe56"

    spec_names(ih1) = "hydrogen-1"
    spec_names(ihe4) = "helium-4"
    spec_names(io14) = "oxygen-14"
    spec_names(io15) = "oxygen-15"
    spec_names(ine18) = "neon-18"
    spec_names(isi25) = "silicon-25"
    spec_names(ife56) = "iron-56"

    aion(ih1) = ONE
    aion(ihe4) = FOUR
    aion(io14) = 14.0_rt
    aion(io15) = 15.0_rt
    aion(ine18) = 18.0_rt
    aion(isi25) = 25.0_rt
    aion(ife56) = 56.0_rt

    zion(ih1) = ONE
    zion(ihe4) = TWO
    zion(io14) = EIGHT
    zion(io15) = EIGHT
    zion(ine18) = TEN
    zion(isi25) = 14.0_rt
    zion(ife56) = 26.0_rt

    ! Binding Energy in MeV
    ! from Stan:
    ! be from http://t2.lanl.gov/nis/data/astro/molnix96/quemd.php
    bion(ih1) = ZERO
    bion(ihe4) = 28.296006_rt
    bion(io14) = 98.733369_rt
    bion(io15) = 111.956652_rt
    bion(ine18) = 132.144124_rt
    bion(isi25) = 187.005541_rt
    bion(ife56) = 492.2450_rt   ! older value, but not important -- this is inert

    ! convert to erg/g by multiplying by N_A / aion and converting to erg
    bion = -bion * N_A * MeV2erg / aion

  end subroutine actual_network_init



  subroutine actual_network_finalize()

    implicit none

  end subroutine actual_network_finalize
  
end module actual_network
