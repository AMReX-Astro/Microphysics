! Common variables and routines for burners
! that use VODE for their integration.

module actual_integrator_module

  use burn_type_module
  use microphysics_type_module

  implicit none
  
contains

  subroutine actual_integrator_init()

    implicit none

  end subroutine actual_integrator_init


  ! Main interface
  subroutine actual_integrator(state_in, state_out, dt, time)

    !$acc routine seq

    implicit none

    ! Input arguments

    type (burn_t), intent(in   ) :: state_in
    type (burn_t), intent(inout) :: state_out
    real(rt),    intent(in   ) :: dt, time

    !$gpu
    
  end subroutine actual_integrator

end module actual_integrator_module
