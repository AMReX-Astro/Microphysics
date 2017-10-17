module actual_burner_module

  use bl_types
  use bl_constants_module
  use network

#ifndef SDC
  use burn_type_module
#else
  use sdc_type_module
#endif

  implicit none

contains

  subroutine actual_burner_init()

    use integrator_module, only: integrator_init

    implicit none

    call integrator_init()

  end subroutine actual_burner_init



  subroutine actual_burner(state_in, state_out, dt, time)

    !$acc routine seq

    use integrator_module, only: integrator

    implicit none

#ifndef SDC
    type (burn_t),    intent(in   ) :: state_in
    type (burn_t),    intent(inout) :: state_out
#else
    type (sdc_t),    intent(in   )  :: state_in
    type (sdc_t),    intent(inout)  :: state_out
#endif
    double precision, intent(in   ) :: dt, time

    call integrator(state_in, state_out, dt, time)

  end subroutine actual_burner

end module actual_burner_module
