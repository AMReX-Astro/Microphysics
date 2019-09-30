module react_zones_module

  use amrex_fort_module, only : rt => amrex_real
  use amrex_constants_module

  implicit none

contains

  subroutine react_test(dt) bind(C, name="react_test")

    use sdc_ode_module, only: ode, sdc_t

    implicit none

    real(rt), intent(in), value :: dt

    type(sdc_t) :: sdc_state
    integer :: ierr

    ! Set the absolute tolerances
    sdc_state % atol(1) = 1.d-8
    sdc_state % atol(2) = 1.d-14
    sdc_state % atol(3) = 1.d-6

    ! Set the relative tolerances
    sdc_state % rtol(:) = 1.d-4

    ! Initialize the integration time and set the final time to dt
    sdc_state % t = ZERO
    sdc_state % tmax = dt

    ! Initialize the initial conditions
    sdc_state % y(:) = 0.0_rt
    sdc_state % y(1) = 1.0_rt

    ! Call the integration routine.
    call ode(sdc_state, ierr)

    ! Check if the integration failed
    if (ierr  < 0) then
       print *, 'ERROR: integration failed', ierr
       stop
    endif

    ! print the final result
    print *, sdc_state % y(:)

  end subroutine react_test

end module react_zones_module
