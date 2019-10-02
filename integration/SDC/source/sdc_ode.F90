module sdc_ode_module

  use amrex_constants_module
  use amrex_fort_module, only : rt => amrex_real
  use rhs_module
  use sdc_sizes_module
  use sdc_quadrature_module
  use sdc_parameters_module
  use sdc_rpar_indices, only : N_RPAR_COMPS

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
  integer, parameter :: IERR_TOO_MANY_DT_ATTEMPTS = -104

  integer, parameter :: IERR_LU_DECOMPOSITION_ERROR = -200

  ! this is the number of time step attempts we try
  integer, parameter :: MAX_TIMESTEP_ATTEMPTS = 50

  ! timestep control
  real, parameter :: S1 = 0.9_rt
  real, parameter :: S2 = 10.0_rt

  type :: sdc_t
     real(rt) :: rpar(N_RPAR_COMPS)
     real(rt) :: rtol(SDC_NEQS)
     real(rt) :: atol(SDC_NEQS)
     real(rt) :: y(SDC_NEQS)

     real(rt) :: t
     real(rt) :: tmax
     real(rt) :: dt

     ! if dt_ini < 0 then we estimate the initial timestep
     real(rt) :: dt_ini

     integer :: n

  end type sdc_t

contains

  ! integrate from t to tmax


  subroutine ode(sdc, idiag, ierr)

    ! this is a basic driver for the ODE integration, based on the NR
    ! routine.  This calls an integration method to take a single step
    ! and return an estimate of the net step size needed to achieve
    ! our desired tolerance.

    !$acc routine seq

    use extern_probin_module, only: ode_max_steps
    use amrex_error_module, only: amrex_error

    type (sdc_t), intent(inout) :: sdc
    integer, intent(out) :: ierr
    type (sdc_diag_t), intent(out) :: idiag

    real(rt) :: ydot(SDC_NEQS)
    logical :: finished
    real(rt) :: t_start, tmax
    integer :: n

    ! initialize the RHS call diags
    idiag % count = 0
    idiag % retries = 0

    ! store local copies
    t_start = sdc % t
    tmax = sdc % tmax

    finished = .false.
    ierr = IERR_NONE

    ! estimate the initial timestep
    if (sdc % dt_ini < 0.0_rt) then
       call f_rhs(sdc % t, sdc % y, ydot, sdc % rpar)
       idiag % count = idiag % count + 1
       call initial_timestep(sdc, idiag, ydot)
    else
       sdc % dt = sdc % dt_ini
    end if

    print *, "initial timestep = ", sdc % dt

    do n = 1, ode_max_steps

       ! make sure we don't overshoot the ending time
       if (sdc % t + sdc % dt > tmax) sdc % dt = tmax - sdc % t

       ! take a step -- this routine will update the solution array,
       ! advance the time, and also give an estimate of the next step
       ! size
       call adaptive_step_driver(sdc, idiag, ierr)

       if (ierr /= IERR_NONE) then
          exit
       end if

       ! finished?
       if (sdc % t - tmax >= ZERO) then
          finished = .true.
       endif

       !sdc % dt = sdc % dt_next

       if (sdc % dt < dt_min) then
          ierr = IERR_DT_TOO_SMALL
          exit
       endif

       ! each "step" is actually 2 half steps
       sdc % n = 2*n

       if (finished) exit

    enddo

    ! store an estimate for next time
    sdc % dt_ini = sdc % dt

    if (.not. finished .and. ierr == IERR_NONE) then
       ierr = IERR_TOO_MANY_STEPS
    endif

    if (ierr /= IERR_NONE) then
       print *, "errors encountered"
    end if

  end subroutine ode


  subroutine initial_timestep(sdc, idiag, ydot)

    ! this is a version of the timestep estimation algorithm used by
    ! VODE

    !$acc routine seq

    type (sdc_t), intent(inout) :: sdc
    real(rt), intent(in) :: ydot(SDC_NEQS)
    type (sdc_diag_t), intent(inout) :: idiag

    real(rt) :: h, h_old, hL, hU, ddydtt(SDC_NEQS), eps, ewt(SDC_NEQS), yddnorm
    real(rt) :: t_temp
    real(rt) :: y_temp(SDC_NEQS)
    real(rt) :: ydot_temp(SDC_NEQS)
    integer :: n

    y_temp = sdc % y

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

       t_temp = sdc % t + h
       y_temp = sdc % y + h * ydot

       ! Call the RHS, then estimate the finite difference.
       call f_rhs(t_temp, y_temp, ydot_temp, sdc % rpar)
       idiag % count = idiag % count + 1

       ddydtt = (ydot_temp - ydot) / h


       yddnorm = sqrt( sum( (ddydtt*ewt)**2 ) / SDC_NEQS )

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


  subroutine adaptive_step_driver(sdc, idiag, ierr)

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

    implicit none

    type (sdc_t), intent(inout) :: sdc
    integer, intent(out) :: ierr
    type (sdc_diag_t), intent(inout) :: idiag

    type (sdc_t) :: sdc_s, sdc_d

    real(rt) :: y_save(SDC_NEQS), y_s(SDC_NEQS), y_d(SDC_NEQS)
    real(rt) :: t_save, dt_save
    real(rt) :: eps
    logical :: converged

    real(rt) :: weights(SDC_NEQS)
    real(rt) :: dt, dt_est

    integer :: i

    ! store the initial state
    t_save = sdc % t
    dt_save = sdc % dt
    y_save(:) = sdc % y(:)

    ! main attempt loop
    dt = dt_save
    converged = .false.

    weights(:) = 1.0_rt / (sdc % rtol(:) * abs(y_save(:)) + sdc % atol(:))

    do i = 1, MAX_TIMESTEP_ATTEMPTS

       ! take a single step of size dt
       sdc_s = sdc
       sdc_s % t = t_save
       sdc_s % dt = dt
       sdc_s % y(:) = y_save(:)

       call single_step_sdc(sdc_s, idiag, ierr)

       ! if this was successful, the in all likelihood, the two
       ! integrations with dt/2 will be too.  Otherwise, bail here and
       ! cut the timestep
       if (ierr /= IERR_NONE) then
          ! maybe we can do something more smartly here later
          dt = 0.5d0 * dt
          cycle
       end if

       ! now take two steps of since dt/2
       sdc_d = sdc
       sdc_d % t = t_save
       sdc_d % dt = 0.5d0 * dt
       sdc_d % y(:) = y_save(:)

       call single_step_sdc(sdc_d, idiag, ierr)

       if (ierr /= IERR_NONE) then
          dt = 0.5d0 * dt
          cycle
       end if

       call single_step_sdc(sdc_d, idiag, ierr)

       if (ierr /= IERR_NONE) then
          dt = 0.5d0 * dt
          cycle
       end if

       ! now the single dt and the two dt/2 attempts are at the same
       ! time.  check convergence as |y_s - y_d| < rtol * |y_old| + atol
       eps = sqrt(sum(weights(:) * abs(sdc_s % y(:) - sdc_d % y(:))) / SDC_NEQS)

       ! predict the new timestep -- if we converged, we'll pass this
       ! out.  Otherwise, we will use this for the next attempt.
       dt_est = dt * eps**(-0.2_rt)
       dt = min(max(S1 * dt_est, dt/S2), S2 * dt)

       if (eps < 1.0_rt) then
          ! we converged!!!
          ierr = IERR_NONE
          converged = .true.
          exit
       end if

    end do

    if (.not. converged) then
       ierr = IERR_TOO_MANY_DT_ATTEMPTS
    end if

    ! Store the solution
    sdc % y(:) = sdc_d % y(:)
    sdc % t = sdc_d % t
    sdc % dt = dt

  end subroutine adaptive_step_driver

  subroutine single_step_sdc(sdc, idiag, ierr)

    ! for a given state and timestep (encoded in sdc_t), we do the
    ! update to the new time

    use amrex_error_module, only: amrex_error

    implicit none

    type (sdc_t), intent(inout) :: sdc
    integer, intent(out) :: ierr
    type (sdc_diag_t), intent(inout) :: idiag

    ! storage for the solution on each time node
    real(rt) :: y_node(0:SDC_NODES-1, SDC_NEQS)

    ! storage for the previous iteration's RHS
    real(rt) :: f_old(0:SDC_NODES-1, SDC_NEQS)
    real(rt) :: ydot(SDC_NEQS)
    real(rt) :: Z_source(SDC_NEQS), C(SDC_NEQS)
    real(rt) :: weights(SDC_NEQS)

    real(rt) :: dt, dt_m
    real(rt) :: t_start
    real(rt) :: err_max

    integer :: i, k, m, n
    integer :: inewton_err

    dt = sdc % dt

    t_start = sdc % t

    ! set the error to "all clear"
    ierr = IERR_NONE

    do m = 0, SDC_NODES-1
       y_node(m, :) = sdc % y(:)
    end do

    ! compute an estimate for the RHS at the "old" iteration
    call f_rhs(sdc % t, sdc % y, ydot, sdc % rpar)
    idiag % count = idiag % count + 1

    do m = 0, SDC_NODES-1
       f_old(m, :) = ydot(:)
    end do

    do k = 1, SDC_MAX_ITERATIONS

       ! loop over time nodes
       do m = 0, SDC_NODES-2

          ! our starting point is y_node(m, :)

          dt_m = dt * (dt_sdc(m+1) - dt_sdc(m))

          ! compute the integral term
          call sdc_C_source(m, dt, dt_m, f_old, C)

          ! compute the explicit source
          Z_source(:) = y_node(m, :) + dt_m * C(:)

          ! initial guess
          y_node(m+1, :) = y_node(m, :) + dt_m * f_old(m+1, :)

          ! Mike's initial guess
          y_node(m+1, :) = Z_source(:) + dt_m * f_old(m, :)

          ! do the nonlinear solve to find the solution at time node m+1
          weights(:) = 1.0_rt/(sdc % rtol(:) * abs(y_node(m, :)) + sdc % atol(:))

          call sdc_newton_solve(t_start + dt_sdc(m+1), dt_m, &
                                y_node(m+1, :), Z_source, k, sdc % rpar, weights, &
                                idiag, inewton_err)

          ! did the solve converge?
          if (inewton_err /= NEWTON_SUCCESS) then
             ierr = IERR_NO_CONVERGENCE
             exit
          end if

       end do

       if (ierr == IERR_NONE) then
          ! recompute f on all time nodes and store
          do m = 0, SDC_NODES-1
             call f_rhs(sdc % t + dt*dt_sdc(m), y_node(m, :), ydot, sdc % rpar)
             idiag % count = idiag % count + 1
             f_old(m, :) = ydot(:)
          end do

       else
          ! no reason to keep going -- bail out to the caller and have
          ! it handle things
          exit
       end if

    end do

    if (ierr == IERR_NONE) then
       sdc % t = t_start + dt
       sdc % y = y_node(SDC_NODES-1, :)
    end if

  end subroutine single_step_sdc

end module sdc_ode_module
