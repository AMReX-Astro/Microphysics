! Common variables and routines for burners
! that use VODE for their integration.

module vode_integrator_module

  use network
  use burn_type_module
  use amrex_fort_module, only: rt => amrex_real
  use amrex_error_module
  use integration_data

  implicit none
  
contains

  subroutine vode_integrator_init()

    implicit none

  end subroutine vode_integrator_init


  ! Main interface
  subroutine vode_integrator(state_in, state_out, dt, time, status)

    implicit none

    ! Input arguments

    type (burn_t), intent(in   ) :: state_in
    type (burn_t), intent(inout) :: state_out
    real(rt),    intent(in   ) :: dt, time
    type (integration_status_t), intent(inout) :: status

    call amrex_error("Error: VODE is no longer supported for Fortran")

  end subroutine vode_integrator

end module vode_integrator_module
