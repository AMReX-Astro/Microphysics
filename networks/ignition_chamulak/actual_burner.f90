module actual_burner_module

  use burn_type_module
  use actual_network

contains

  subroutine actual_burner_init()

    use integrator_module, only: integrator_init

    implicit none

    call integrator_init()

  end subroutine actual_burner_init



  subroutine actual_burner(state_in, state_out, dt, time)

    use integrator_module, only: integrator

    implicit none

    type (burn_t),    intent(in   ) :: state_in
    type (burn_t),    intent(inout) :: state_out
    double precision, intent(in   ) :: dt, time

    call integrator(state_in, state_out, dt, time)

  end subroutine actual_burner

end module actual_burner_module
