! The true SDC integration scheme does not use this integrator
! framework

module actual_integrator_module

  use burn_type_module

  implicit none

contains

  subroutine actual_integrator_init()

    implicit none

  end subroutine actual_integrator_init


  subroutine actual_integrator(state_in, state_out, dt, time)

    implicit none

    ! Input arguments

    type (burn_t), intent(in   ) :: state_in
    type (burn_t), intent(inout) :: state_out
    real(rt),    intent(in   ) :: dt, time

  end subroutine actual_integrator

end module actual_integrator_module
