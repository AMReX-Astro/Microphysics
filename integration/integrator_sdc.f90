module integrator_module

  implicit none

  public

contains

  subroutine integrator_init()

    use actual_integrator_module, only: actual_integrator_init
    use integration_data, only: aionInv
    use actual_network, only: aion
    use bl_constants_module, only: ONE

    implicit none

    call actual_integrator_init()

    aionInv = ONE / aion

    !$acc update device(aionInv)

  end subroutine integrator_init



  subroutine integrator(state_in, state_out, dt, time)

    !$acc routine seq

    use actual_integrator_module, only: actual_integrator
    use bl_error_module, only: bl_error
    use bl_types, only: dp_t
    use sdc_type_module, only: sdc_t

    implicit none

    type (sdc_t),  intent(in   ) :: state_in
    type (sdc_t),  intent(inout) :: state_out
    real(dp_t),    intent(in   ) :: dt, time

    call actual_integrator(state_in, state_out, dt, time)

  end subroutine integrator

end module integrator_module
