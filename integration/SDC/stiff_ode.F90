module stiff_ode

  use amrex_constants_module
  use amrex_fort_module, only : rt => amrex_real
  use burn_type_module
  use sdc_type_module
#ifdef SDC
  use sdc_rhs_module
  use sdc_jac_module
#else
  use rhs_module
  use jac_module
#endif

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

    sdc % eps_old = ZERO

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

       ewt = sdc % rtol(:) * abs(sdc % y) + sdc & atol(:)

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


  subroutine single_step_sdc(sdc, eps, yscal, ierr)

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
    real(rt), intent(in) :: yscal(sdc_neqs)
    integer, intent(out) :: ierr

    real(rt) :: y_save(sdc_neqs), yerr(sdc_neqs), yseq(sdc_neqs)
    real(rt) :: err(KMAXX)

    real(rt) :: dt, fac, scale, red, eps1, work, work_min, xest
    real(rt) :: err_max

    logical :: converged, reduce, skip_loop, retry

    integer :: i, k, n, kk, km, kstop, ierr_temp
    integer, parameter :: max_iters = 10 ! Should not need more than this


    dt = bs % dt
    y_save(:) = bs % y(:)

    ! get the jacobian
    call jac(sdc)

    reduce = .false.

    converged = .false.

    do n = 1, max_iters

       ! setting skip_loop = .true. in the next loop is a GPU-safe-way to
       ! indicate that we are discarding this timestep attempt and will
       ! instead try again in the next `n` iteration
       skip_loop = .false.
       retry = .false.

       ! each iteration is a new attempt at taking a step, so reset
       ! errors at the start of the attempt
       ierr = IERR_NONE

       do k = 1, bs % kmax

          if (.not. skip_loop) then

             bs % t_new = bs % t + dt

             call semi_implicit_extrap(bs, y_save, dt, nseq(k), yseq, ierr_temp)
             ierr = ierr_temp
             if (ierr == IERR_LU_DECOMPOSITION_ERROR) then
                skip_loop = .true.
                red = ONE/nseq(k)
             endif

             call safety_check(y_save, yseq, retry)
             if (retry) then
                skip_loop = .true.
                red = RED_BIG_FACTOR
             endif

             if (ierr == IERR_NONE .and. .not. retry) then
                xest = (dt/nseq(k))**2
                call poly_extrap(k, xest, yseq, bs % y, yerr, t_extrap, qcol)

                if (k /= 1) then
                   err_max = max(SMALL, maxval(abs(yerr(:)/yscal(:))))
                   err_max = err_max / eps
                   km = k - 1
                   err(km) = (err_max/S1)**(1.0/(2*km+1))
                endif

                if (k /= 1 .and. (k >=  bs % kopt-1 .or. bs % first)) then

                   if (err_max < 1) then
                      converged = .true.
                      kstop = k
                      exit
                   else

                      ! reduce stepsize if necessary
                      if (k == bs % kmax .or. k == bs % kopt+1) then
                         red = S2/err(km)
                         reduce = .true.
                         skip_loop = .true.
                      else if (k == bs % kopt) then
                         if (bs % alpha(bs % kopt-1, bs % kopt) < err(km)) then
                            red = ONE/err(km)
                            reduce = .true.
                            skip_loop = .true.
                         endif
                      else if (bs % kopt == bs % kmax) then
                         if (bs % alpha(km, bs % kmax-1) < err(km)) then
                            red = bs % alpha(km, bs % kmax-1)*S2/err(km)
                            reduce = .true.
                            skip_loop = .true.
                         endif
                      else if (bs % alpha(km, bs % kopt) < err(km)) then
                         red = bs % alpha(km, bs % kopt - 1)/err(km)
                         reduce = .true.
                         skip_loop = .true.
                      endif
                   endif

                endif

                kstop = k
             endif
          endif

       enddo

       if (.not. converged) then
          ! note, even if ierr /= IERR_NONE, we still try again, since
          ! we may eliminate LU decomposition errors (singular matrix)
          ! with a smaller timestep
          red = max(min(red, RED_BIG_FACTOR), RED_SMALL_FACTOR)
          dt = dt*red
       else
          exit
       endif

    enddo   ! iteration loop (n) varying dt

    if (.not. converged .and. ierr == IERR_NONE) then
       ierr = IERR_NO_CONVERGENCE
#ifndef ACC
       print *, "Integration failed due to non-convergence in single_step_bs"
       call dump_bs_state(bs)
       return
#endif
    endif

#ifndef ACC
    ! km and kstop should have been set during the main loop.
    ! If they never got updated from the original nonsense values,
    ! that means something went really wrong and we should abort.

    if (km < 0) then
       call amrex_error("Error: km < 0 in subroutine single_step_bs, something has gone wrong.")
    endif

    if (kstop < 0) then
       call amrex_error("Error: kstop < 0 in subroutine single_step_bs, something has gone wrong.")
    endif
#endif

    bs % t = bs % t_new
    bs % dt_did = dt
    bs % first = .false.

    ! optimal convergence properties
    work_min = 1.e35
    do kk = 1, km
       fac = max(err(kk), SCALMX)
       work = fac*bs % a(kk+1)

       if (work < work_min) then
          scale = fac
          work_min = work
          bs % kopt = kk+1
       endif
    enddo

    ! increase in order
    bs % dt_next = dt / scale

    if (bs % kopt >= kstop .and. bs % kopt /= bs % kmax .and. .not. reduce) then
       fac = max(scale/bs % alpha(bs % kopt-1, bs % kopt), SCALMX)
       if (bs % a(bs % kopt+1)*fac <= work_min) then
          bs % dt_next = dt/fac
          bs % kopt = bs % kopt + 1
       endif
    endif

  end subroutine single_step_bs

end module stiff_ode
