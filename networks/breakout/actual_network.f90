! network module
!
! nspec -> number of species
! naux  -> number of auxiliary variables also used in the EOS
!
! aion -> atomic number
! zion -> proton number
!

module actual_network

  use amrex_fort_module, only : rt => amrex_real

  implicit none

  integer, parameter :: nspec = 1
  integer, parameter :: naux  = 2
  integer, parameter :: nspec_evolve = 1
  integer, parameter :: nrates = 0
  integer, parameter :: num_rate_groups = 0

  character (len=16), save :: spec_names(nspec)
  character (len= 5), save :: short_spec_names(nspec)

  character (len=16), save ::  aux_names(naux)
  character (len= 5), save :: short_aux_names(naux)

  real(rt), save :: aion(nspec), zion(nspec), ebin(nspec)

contains

  subroutine actual_network_init()

    implicit none

    spec_names(1) = "X"
    short_spec_names(1) = "X"

    aux_names(1) = "Ye"
    short_aux_names(1) = "Ye"

    aux_names(2) = "invmu"   ! 1/mu in P = rho*R*T/mu
    short_aux_names(2) = "invmu"

    aion(:) = 1.e0_rt
    zion(:) = 1.e0_rt
    ebin(:) = 0.e0_rt

  end subroutine actual_network_init



  subroutine actual_network_finalize()

    implicit none

  end subroutine actual_network_finalize

end module actual_network
