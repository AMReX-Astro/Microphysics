module actual_burner_module

  use bl_types
  use bl_constants_module
  use bl_error_module
  use eos_module
  use network
  use burn_type_module

  implicit none

contains

  subroutine actual_burner(state_in, state_out, dt, time)

    use integrator_module, only: integrator

    implicit none

    type (burn_t),       intent(in   ) :: state_in
    type (burn_t),       intent(inout) :: state_out
    double precision,    intent(in   ) :: dt, time

    call integrator(state_in, state_out, dt, time)

  end subroutine actual_burner



  subroutine actual_burner_init()

    use integrator_module, only: integrator_init
    use integration_data, only: ener_scale

    implicit none

    ener_scale = c_light * c_light

    call integrator_init()

  end subroutine actual_burner_init

end module actual_burner_module
