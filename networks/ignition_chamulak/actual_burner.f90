module actual_burner_module

  use bl_constants_module
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

    state_out % xn(io16_) = state_in % xn(io16_)
    state_out % xn(iash_) = (ONE - state_out % xn(ic12_) - state_out % xn(io16_))

  end subroutine actual_burner

end module actual_burner_module
