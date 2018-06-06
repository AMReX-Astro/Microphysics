module actual_burner_module

  use bl_types
  use bl_constants_module
  use bl_error_module
  use eos_module
  use network
  use burn_type_module

  implicit none

  interface actual_burner
     module procedure actual_burner1
     module procedure actual_burnerv
  end interface actual_burner

contains

  subroutine actual_burner1(state_in, state_out, dt, time)

    !$acc routine seq

    use integrator_module, only: integrator

    implicit none

    type (burn_t),       intent(in   ) :: state_in
    type (burn_t),       intent(inout) :: state_out
    double precision,    intent(in   ) :: dt, time

    call integrator(state_in, state_out, dt, time)

  end subroutine actual_burner1

  subroutine actual_burnerv(state_in, state_out, dt, time)

    use integrator_module, only: integrator

    implicit none

    type (burn_t),       intent(in   ) :: state_in(:)
    type (burn_t),       intent(inout) :: state_out(:)
    double precision,    intent(in   ) :: dt, time

    call integrator(state_in, state_out, dt, time)

  end subroutine actual_burnerv


  subroutine actual_burner_init()

    use integrator_module, only: integrator_init

    implicit none

    call integrator_init()

  end subroutine actual_burner_init

end module actual_burner_module
