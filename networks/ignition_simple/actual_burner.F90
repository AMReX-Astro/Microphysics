module actual_burner_module

  use amrex_constants_module
  use network
  use burn_type_module
  use microphysics_type_module

  implicit none

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
    real(rt), intent(in   ) :: dt, time

    !$gpu

    call integrator(state_in, state_out, dt, time)

  end subroutine actual_burner

end module actual_burner_module
