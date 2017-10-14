module integrator_module

  implicit none

  public

contains

  subroutine integrator_init()

#if (INTEGRATOR == 0 || INTEGRATOR == 1)
    use vode_integrator_module, only: vode_integrator_init
    use bs_integrator_module, only: bs_integrator_init
#else
    use actual_integrator_module, only: actual_integrator_init
#endif

    implicit none

#if (INTEGRATOR == 0 || INTEGRATOR == 1)
    call vode_integrator_init()
    call bs_integrator_init()
#else
    call actual_integrator_init()
#endif

  end subroutine integrator_init



  subroutine integrator(state_in, state_out, dt, time)

    !$acc routine seq

#if (INTEGRATOR == 0 || INTEGRATOR == 1)
    use vode_integrator_module, only: vode_integrator
    use bs_integrator_module, only: bs_integrator
#else
    use actual_integrator_module, only: actual_integrator_init
#endif

    use bl_error_module, only: bl_error
    use bl_types, only: dp_t
    use integration_data, only: integration_status_t
    use sdc_type_module, only: sdc_t

    implicit none

    type (sdc_t),  intent(in   ) :: state_in
    type (sdc_t),  intent(inout) :: state_out
    real(dp_t),    intent(in   ) :: dt, time

    type (integration_status_t) :: status

#if INTEGRATOR == 0
    call vode_integrator(state_in, state_out, dt, time, status)
#elif INTEGRATOR == 1
    call bs_integrator(state_in, state_out, dt, time, status)
#else
    call bl_error("Unknown integrator.")
#endif

  end subroutine integrator

end module integrator_module
