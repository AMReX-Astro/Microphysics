! network module
!
! nspec -> number of species
! naux  -> number of auxiliary variables also used in the EOS
!
! aion -> atomic number
! zion -> proton number
!

module actual_network

  use network_properties
  use amrex_fort_module, only : rt => amrex_real

  implicit none

  integer, parameter :: nspec_evolve = 1
  integer, parameter :: nrates = 0
  integer, parameter :: num_rate_groups = 0

  real(rt), save :: ebin(nspec)

contains

  subroutine actual_network_init()

    implicit none

    call network_properties_init()

    ebin(:) = 0.e0_rt

  end subroutine actual_network_init



  subroutine actual_network_finalize()

    implicit none

    call network_properties_finalize()

  end subroutine actual_network_finalize

end module actual_network
