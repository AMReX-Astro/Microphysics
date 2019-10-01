module react_zones_module

  use amrex_fort_module, only : rt => amrex_real
  use amrex_constants_module

  implicit none

contains

  subroutine react_test() bind(C, name="react_test")

    use sdc_ode_module, only: ode, sdc_t

    implicit none

    type(sdc_t) :: sdc_state
    integer :: ierr
    integer :: i

    ! Set the absolute tolerances
    sdc_state % atol(1) = 1.d-10
    sdc_state % atol(2) = 1.d-14
    sdc_state % atol(3) = 1.d-8

    ! Set the relative tolerances
    sdc_state % rtol(:) = 1.d-7

    ! Initialize the integration time and set the final time to dt
    sdc_state % t = ZERO

    ! Initialize the initial conditions
    sdc_state % y(:) = 0.0_rt
    sdc_state % y(1) = 1.0_rt

    sdc_state % tmax = 0.4_rt

    ! for a timestep estimation
    sdc_state % dt_ini = -1.0_rt

    do i = 1, 12

       ! Call the integration routine.
       call ode(sdc_state, ierr)

       ! Check if the integration failed
       if (ierr  < 0) then
          print *, "ERROR: integration failed", ierr
          print *, sdc_state % t
          stop
       endif

       print *, sdc_state % t, sdc_state % n, sdc_state % y(:)

       sdc_state % tmax = 10.0_rt * sdc_state % tmax

    end do

    print *, "number of steps = ", sdc_state % n
  end subroutine react_test

end module react_zones_module
