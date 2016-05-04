module actual_burner_module

  use bl_types
  use bl_constants_module
  use network
  use burn_type_module
  use actual_burner_data

  implicit none

contains

  subroutine actual_burner_init()

    use integrator_module, only: integrator_init
    use integration_data, only: temp_scale, ener_scale

    implicit none

    call integrator_init()

    temp_scale = 1.0d10
    ener_scale = c_light * c_light

    !$acc update device(temp_scale, ener_scale)

  end subroutine actual_burner_init



  subroutine actual_burner(state_in, state_out, dt, time)

    !$acc routine seq

    use integrator_module, only: integrator

    implicit none

    type (burn_t),    intent(in   ) :: state_in
    type (burn_t),    intent(inout) :: state_out
    double precision, intent(in   ) :: dt, time

    call integrator(state_in, state_out, dt, time)

  end subroutine actual_burner

end module actual_burner_module
