module actual_network

  use network_properties
  use amrex_fort_module, only : rt => amrex_real

  implicit none

  integer, parameter :: nrates = 6  ! including constant weak rates
  integer, parameter :: num_rate_groups = 1

  real(rt)        , save :: bion(nspec)
  real(rt)         :: rates(nrates)

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
  real(rt)        , parameter :: wk14o = 9.832e-3_rt
  real(rt)        , parameter :: wk15o = 5.682e-3_rt

contains

  subroutine actual_network_init()
    
    use amrex_constants_module, only: ZERO, ONE, TWO, FOUR, EIGHT, TEN

    real(rt), parameter :: MeV2erg = 1.60217646e-6_rt, N_A = 6.0221415e23_rt

    call network_properties_init()

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

    call network_properties_finalize()

  end subroutine actual_network_finalize
  
end module actual_network
