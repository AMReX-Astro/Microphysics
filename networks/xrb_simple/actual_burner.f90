module actual_burner_module

  ! From Woosley:
  ! approximation network for the early stages of runaway in hydrogen-rich
  ! bursts during their convective stage
  !
  ! 6 network 14o, 15o, 18ne, 25si, alpha, p
  ! - s. woosley 08/26/2015

contains

  subroutine actual_burner_init()

    use integrator_module, only: integrator_init

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    call integrator_init()

  end subroutine actual_burner_init



  subroutine actual_burner(state_in, state_out, dt, time)

    use amrex_fort_module, only : rt => amrex_real
    use integrator_module, only: integrator
    use burn_type_module, only: burn_t

    implicit none

    type (burn_t) :: state_in, state_out
    real(rt)      :: dt, time

    call integrator(state_in, state_out, dt, time)

  end subroutine actual_burner

end module actual_burner_module
