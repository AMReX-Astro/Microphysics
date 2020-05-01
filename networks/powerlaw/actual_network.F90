module actual_network

  use network_properties
  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer, parameter :: ifuel_  = 1
  integer, parameter :: iash_   = 2
  integer, parameter :: iinert_ = 3


  real(rt), allocatable, save :: ebin(:)

  integer, parameter :: nrates = 0
  integer, parameter :: num_rate_groups = 0

  character (len=32), parameter :: network_name = "powerlaw"

#ifdef AMREX_USE_CUDA
  attributes(managed) :: ebin
#endif

contains
  
  subroutine actual_network_init

    use extern_probin_module, only: specific_q_burn
    use amrex_constants_module, only: ZERO
    use fundamental_constants_module, only: N_A

    call network_properties_init()

    ! Binding energies in erg / g

    allocate(ebin(3))

    ebin(ifuel_)  = ZERO
    ebin(iash_)   = specific_q_burn
    ebin(iinert_) = ZERO

  end subroutine actual_network_init



  subroutine actual_network_finalize

    implicit none

    call network_properties_finalize()

    deallocate(ebin)

  end subroutine actual_network_finalize

end module actual_network
