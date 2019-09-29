module sdc_ode_module

  use amrex_constants_module
  use amrex_fort_module, only : rt => amrex_real
  use rhs_module
  use sdc_sizes_module
  use sdc_quadrature_module

  implicit none

  real(rt), parameter, private :: dt_min = 1.d-24
  real(rt), parameter, private :: dt_ini = 1.d-16
  real(rt), parameter, private :: SMALL = 1.d-30


  ! error codes
  integer, parameter :: IERR_NONE = 0
  integer, parameter :: IERR_DT_TOO_SMALL = -100
  integer, parameter :: IERR_TOO_MANY_STEPS = -101
  integer, parameter :: IERR_DT_UNDERFLOW = -102
  integer, parameter :: IERR_NO_CONVERGENCE = -103

  integer, parameter :: IERR_LU_DECOMPOSITION_ERROR = -200

  integer, parameter :: MAX_TRY = 50

contains

  ! integrate from t to tmax

  subroutine safety_check(y_old, y, retry)
    !$acc routine seq

    use extern_probin_module, only: safety_factor

    real(rt), intent(in) :: y_old(sdc_neqs), y(sdc_neqs)
    logical, intent(out) :: retry

    real(rt) :: ratio

    retry = .false.

    ratio = abs(y(net_itemp)/y_old(net_itemp))
    if (ratio > safety_factor) then ! .or. ratio < ONE/safety_factor) then
       retry = .true.
    endif

    ratio = abs(y(net_ienuc)/y_old(net_ienuc))
    if (ratio > safety_factor) then ! .or. ratio < ONE/safety_factor) then
       retry = .true.
    endif

    ! not sure what check to do on species

  end subroutine safety_check


  subroutine ode(sdc, t, tmax, eps, ierr)

    ! this is a basic driver for the ODE integration, based on the NR
    ! routine.  This calls an integration method to take a single step
    ! and return an estimate of the net step size needed to achieve
    ! our desired tolerance.

    !$acc routine seq

    use extern_probin_module, only: ode_max_steps, use_timestep_estimator, &
                                    scaling_method, ode_scale_floor, ode_method
#ifndef ACC
    use amrex_error_module, only: amrex_error
#endif

    type (sdc_t), intent(inout) :: sdc

    real(rt), intent(inout) :: t
    real(rt), intent(in) :: tmax
    real(rt), intent(in) :: eps
    integer, intent(out) :: ierr

    logical :: finished

    integer :: n

    ! initialize

    sdc % t = t
    sdc % tmax = tmax

    finished = .false.
    ierr = IERR_NONE

    if (use_timestep_estimator) then
       call f_rhs(sdc)
       call initial_timestep(sdc)
    else
       sdc % dt = dt_ini
    endif

    do n = 1, ode_max_steps

       ! make sure we don't overshoot the ending time
       if (sdc % t + sdc % dt > tmax) sdc % dt = tmax - sdc % t

       ! take a step -- this routine will update the solution array,
       ! advance the time, and also give an estimate of the next step
       ! size
       call single_step_sdc(sdc, eps, ierr)

       if (ierr /= IERR_NONE) then
          exit
       end if

       ! finished?
       if (sdc % t - tmax >= ZERO) then
          finished = .true.
       endif

       sdc % dt = sdc % dt_next

       if (sdc % dt < dt_min) then
          ierr = IERR_DT_TOO_SMALL
          exit
       endif

       if (finished) exit

    enddo

    sdc % n = n

    if (.not. finished .and. ierr == IERR_NONE) then
       ierr = IERR_TOO_MANY_STEPS
    endif

  end subroutine ode


  subroutine initial_timestep(sdc)

    ! this is a version of the timestep estimation algorithm used by
    ! VODE

    !$acc routine seq

    type (sdc_t), intent(inout) :: sdc

    type (sdc_t) :: sdc_temp
    real(rt) :: h, h_old, hL, hU, ddydtt(sdc_neqs), eps, ewt(sdc_neqs), yddnorm
    integer :: n

    sdc_temp = sdc

    ! Initial lower and upper bounds on the timestep

    hL = 100.0d0 * epsilon(ONE) * max(abs(sdc % t), abs(sdc % tmax))
    hU = 0.1d0 * abs(sdc % tmax - sdc % t)

    ! Initial guess for the iteration

    h = sqrt(hL * hU)
    h_old = 10.0 * h

    ! Iterate on ddydtt = (RHS(t + h, y + h * dydt) - dydt) / h

    do n = 1, 4

       h_old = h

       ! Get the error weighting -- this is similar to VODE's dewset
       ! routine

       ewt = sdc % rtol(:) * abs(sdc % y) + sdc % atol(:)

       ! Construct the trial point.

       sdc_temp % t = sdc % t + h
#ifdef SDC
       sdc_temp % y = sdc % y + h * sdc % ydot
#else
       sdc_temp % y = sdc % y + h * sdc % burn_s % ydot
#endif

       ! Call the RHS, then estimate the finite difference.
       call f_rhs(sdc_temp)
#ifdef SDC
       ddydtt = (sdc_temp % ydot - sdc % ydot) / h
#else
       ddydtt = (sdc_temp % burn_s % ydot - sdc % burn_s % ydot) / h
#endif

       yddnorm = sqrt( sum( (ddydtt*ewt)**2 ) / sdc_neqs )

       if (yddnorm*hU*hU > TWO) then
          h = sqrt(TWO / yddnorm)
       else
          h = sqrt(h * hU)
       endif

       if (h_old < TWO * h .and. h_old > HALF * h) exit

    enddo

    ! Save the final timestep, with a bias factor.

    sdc % dt = h / TWO
    sdc % dt = min(max(h, hL), hU)

  end subroutine initial_timestep


  subroutine single_step_sdc(sdc, eps, ierr)

    ! This routine will do the following:
    !
    ! * we come in with an initial state, sdc % y(:)
    !
    ! * we integrate using SDC-4 with a single step of size sdc % dt
    !
    ! * we do the same integration taking two steps of size sdc % dt / 2
    !
    ! * we compute an error based on this
    !
    !   - if the error is small enough, then:
    !
    !     * we are happy and we store the new solution back in sdc % y(:)
    !
    !     * estimate the next timestep from the expected truncation error
    !
    !   - otherwise, the error was too large, so we redo the step

    use amrex_error_module, only: amrex_error

    implicit none

    type (sdc_t) :: sdc
    real(rt), intent(in) :: eps
    integer, intent(out) :: ierr

    real(rt) :: y_save(sdc_neqs), y_s(sdc_neqs), y_d(sdc_neqs)

    ! storage for the solution on each time node
    real(rt) :: y_node(0:SDC_NODES-1, sdc_neqs)

    ! storage for the previous iteration's RHS
    real(rt) :: f_old(0:SDC_NODES-1, sdc_neqs)

    real(rt) :: f_source(sdc_neqs), C(sdc_neqs)

    real(rt) :: dt, dt_m
    real(rt) :: t_start
    real(rt) :: err_max

    type (sdc_t) :: sdc_temp

    logical :: converged

    integer :: i, k, m, n

    ! this is the number of time step attempts we try
    integer, parameter :: max_timestep_attempts = 10

    dt = sdc % dt
    y_save(:) = sdc % y(:)

    ! get the jacobian
    call jac(sdc)

    converged = .false.

    do n = 1, max_timestep_attempts

       ! each iteration is a new attempt at taking a step, so reset
       ! errors at the start of the attempt
       ierr = IERR_NONE

       t_start = sdc % t

       sdc % t_new = sdc % t + dt

       ! we start integrating from y_save(:)
       sdc_temp = sdc
       sdc_temp % y(:) = y_save(:)

       do m = 0, SDC_NODES-1
          y_node(m, :) = y_save(:)
       end do

       ! compute an estimate for the RHS at the "old" iteration
       call f_rhs(sdc_temp)

       do m = 0, SDC_NODES-1
          f_old(m, :) = sdc_temp % burn_s % ydot(:)
       end do

       do k = 1, SDC_MAX_ITERATIONS


          ! loop over time nodes
          do m = 0, SDC_NODES-2

             ! load the initial data for node m
             sdc_temp % y(:) = y_node(m, :)

             dt_m = dt * (dt_sdc(m+1) - dt_sdc(m))

             ! compute the integral term
             call sdc_C_source(m, dt, dt_m, f_old, C)

             ! compute the explicit source
             f_source(:) = y_node(m, :) + C(:)

             ! do the nonlinear solve to find the solution at time node m+1
             call sdc_newton_solve(dt_m, sdc_temp, y_node(m+1, :), f_source, k, ierr)

             ! did the solve converge?
             
          end do

          ! recompute f on all time nodes and store
          do m = 0, SDC_NODES-1
             sdc_temp % y(:) = y_node(m, :)
             call f_rhs(sdc_temp)

             f_old(m, :) = sdc_temp % burn_s % ydot(:)
          end do

       end do

       ! did we converge? if so, break out

       ! if we didn't converge, reduce the timestep and try again


    end do


    if (.not. converged .and. ierr == IERR_NONE) then
       ierr = IERR_NO_CONVERGENCE
       print *, "Integration failed due to non-convergence in single_step_sdc"
       call dump_sdc_state(sdc)
       return
    endif

    sdc % t = sdc % t_new
    sdc % dt_did = dt

    ! compute the next timestep
    !sdc % dt_next = 

  end subroutine single_step_sdc

end module sdc_ode_module
