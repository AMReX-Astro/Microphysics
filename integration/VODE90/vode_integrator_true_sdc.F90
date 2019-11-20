! The true SDC integration scheme does not use this integrator
! framework

module vode_integrator_module

  use burn_type_module

  implicit none

contains

  subroutine vode_integrator_init()

    implicit none

  end subroutine vode_integrator_init


  subroutine vode_integrator(state_in, state_out, dt, time, status)

    use integration_data

    implicit none

    ! Input arguments

    type (burn_t), intent(in   ) :: state_in
    type (burn_t), intent(inout) :: state_out
    real(rt),    intent(in   ) :: dt, time
    type (integration_status_t), intent(inout) :: status

  end subroutine vode_integrator

end module vode_integrator_module
