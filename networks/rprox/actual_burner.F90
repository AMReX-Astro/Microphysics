! This module contains a version of Wallace and Woosley's (ApJS 45,389
! (1981)) rprox reaction network burner.

module actual_burner_module

  use amrex_constants_module
  use amrex_error_module
  use amrex_fort_module, only : rt => amrex_real
  use eos_module
  use eos_type_module
  use network
  use burn_type_module

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
    real(rt)        , intent(in   ) :: dt, time

    real(rt)         :: T9

    !$gpu

    T9 = state_in % T * 1.e-9_rt

    ! Only burn if 0.2 < T9 < 2.5 or X(H1) > 0.05.
    ! The last restriction is a kludge based on the last paragraph of WW81.

    if ((T9 .gt. 0.2e0_rt .and. T9 .lt. 2.5e0_rt) .or. state_in % xn(ih1) > 0.05e0_rt) then

       ! Call the integration routine.

       call integrator(state_in, state_out, dt, time)

    endif

  end subroutine actual_burner

end module actual_burner_module
