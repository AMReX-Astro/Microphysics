module integrator_module

  implicit none

  public

contains

  subroutine integrator_init()

    use actual_integrator_module, only: actual_integrator_init
    use integration_data, only: temp_scale, ener_scale, dens_scale

    implicit none

    call actual_integrator_init()

    !$acc update device(temp_scale, ener_scale, dens_scale)

  end subroutine integrator_init



  subroutine integrator(state_in, state_out, dt, time)

    !$acc routine seq

    use actual_integrator_module, only: actual_integrator
    use bl_error_module, only: bl_error
    use burn_type_module, only: burn_t

    implicit none

    type (burn_t),       intent(in   ) :: state_in
    type (burn_t),       intent(inout) :: state_out
    double precision,    intent(in   ) :: dt, time

    call actual_integrator(state_in, state_out, dt, time)

  end subroutine integrator

end module integrator_module
