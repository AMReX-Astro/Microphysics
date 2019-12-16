module actual_burner_module

  use amrex_constants_module
  use amrex_fort_module, only : rt => amrex_real

  use network
  use burn_type_module
  use extern_probin_module, only: specific_q_burn

  use amrex_fort_module, only : rt => amrex_real
  implicit none

contains

  subroutine actual_burner_init()

    use integrator_module, only: integrator_init

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    call integrator_init()

  end subroutine actual_burner_init

  subroutine actual_burner(state_in, state_out, dt, time)

    use integrator_module, only: integrator

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    type (burn_t),    intent(in   ) :: state_in
    type (burn_t),    intent(inout) :: state_out
    real(rt)        , intent(in   ) :: dt, time

    call integrator(state_in, state_out, dt, time)

  end subroutine actual_burner

end module actual_burner_module
