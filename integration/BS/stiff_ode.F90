! A variable-sized stepsize driver for the BS integrator,
! based on the NR routine

module stiff_ode

  use bl_constants_module
  use bl_types
  use burn_type_module
  use bs_type_module
  use rhs_module

  implicit none

  real(kind=dp_t), parameter, private :: dt_min = 1.d-24
  real(kind=dp_t), parameter, private :: dt_ini = 1.d-16
  real(kind=dp_t), parameter, private :: SMALL = 1.d-30

contains

  ! integrate from t to tmax

  subroutine ode(bs, t, tmax, eps, ierr)

    !$acc routine seq
    !$acc routine(f_rhs) seq

    use extern_probin_module, only: scaling_method, ode_max_steps, ode_scale_floor, use_timestep_estimator
#ifndef ACC
    use bl_error_module, only: bl_error
#endif

    type (bs_t), intent(inout) :: bs

    real(kind=dp_t), intent(inout) :: t
    real(kind=dp_t), intent(in) :: tmax
    real(kind=dp_t), intent(in) :: eps
    integer, intent(out) :: ierr

    real(kind=dp_t) :: yscal(neqs)
    logical :: finished

    integer :: n

    ! initialize

    bs % t = t
    bs % tmax = tmax

    finished = .false.
    ierr = IERR_NONE

    bs % eps_old = ZERO
    bs % n_rhs = 0
    bs % n_jac = 0

    if (use_timestep_estimator) then
       call f_rhs(bs)
       call initial_timestep(bs)
    else
       bs % dt = dt_ini
    endif

    do n = 1, ode_max_steps

       ! Get the scaling.
       call f_rhs(bs)

       if (scaling_method == 1) then
          yscal(:) = abs(bs % y(:)) + abs(bs % dt * bs % burn_s % ydot(:)) + SMALL

       else if (scaling_method == 2) then
          yscal = max(abs(bs % y(:)), ode_scale_floor)

#ifndef ACC
       else
          call bl_error("Unknown scaling method in subroutine ode.")
#endif
       endif

       ! make sure we don't overshoot the ending time
       if (bs % t + bs % dt > tmax) bs % dt = tmax - bs % t

       ! take a step -- this routine will update the solution vector,
       ! advance the time, and also give an estimate of the next step
       ! size
       call single_step(bs, eps, yscal, ierr)

       ! finished?
       if (bs % t - tmax >= ZERO) then
          finished = .true.
       endif

       bs % dt = bs % dt_next

       if (bs % dt < dt_min) then
          ierr = IERR_DT_TOO_SMALL
       endif

       if (finished) exit

    enddo

    bs % n = n

    if (.not. finished .and. ierr == IERR_NONE) then
       ierr = IERR_TOO_MANY_STEPS
    endif

  end subroutine ode



  subroutine initial_timestep(bs)

    !$acc routine seq
    !$acc routine(f_rhs) seq

    type (bs_t), intent(inout) :: bs

    type (bs_t) :: bs_temp
    real(kind=dp_t) :: h, h_old, hL, hU, ddydtt(neqs), eps, ewt(neqs), yddnorm
    integer :: n

    bs_temp = bs

    eps = maxval(bs % rtol)

    ! Initial lower and upper bounds on the timestep

    hL = 100.0d0 * epsilon(ONE) * max(abs(bs % t), abs(bs % tmax))
    hU = 0.1d0 * abs(bs % tmax - bs % t)

    ! Initial guess for the iteration

    h = sqrt(hL * hU)
    h_old = 10.0 * h

    ! Iterate on ddydtt = (RHS(t + h, y + h * dydt) - dydt) / h

    do n = 1, 4

       h_old = h

       ! Get the error weighting -- this is similar to VODE's dewset
       ! routine

       ewt = eps * abs(bs % y) + SMALL

       ! Construct the trial point.

       bs_temp % t = bs % t + h
       bs_temp % y = bs % y + h * bs % burn_s % ydot

       ! Call the RHS, then estimate the finite difference.
       ! FIXME -- don't we need to disable have_rates here?
       call f_rhs(bs_temp)
       ddydtt = (bs_temp % burn_s % ydot - bs % burn_s % ydot) / h

       yddnorm = sqrt( sum( (ddydtt*ewt)**2 ) / neqs )

       if (yddnorm*hU*hU > TWO) then
          h = sqrt(TWO / yddnorm)
       else
          h = sqrt(h * hU)
       endif

       if (h_old < TWO * h .and. h_old > HALF * h) exit

    enddo

    ! Save the final timestep, with a bias factor.

    bs % dt = h / TWO
    bs % dt = min(max(h, hL), hU)

  end subroutine initial_timestep



  subroutine semi_implicit_extrap(bs, y, dt_tot, N_sub, y_out, ierr)

    !$acc routine seq
    !$acc routine(f_rhs) seq
    !$acc routine(dgesl) seq
    !$acc routine(dgefa) seq

    type (bs_t), intent(inout) :: bs
    real(kind=dp_t), intent(in) :: y(neqs)
    real(kind=dp_t), intent(in) :: dt_tot
    integer, intent(in) :: N_sub
    real(kind=dp_t), intent(out) :: y_out(neqs)
    integer, intent(out) :: ierr

    real(kind=dp_t) :: A(neqs,neqs)
    real(kind=dp_t) :: del(neqs)
    real(kind=dp_t) :: h

    integer :: n

    integer :: ipiv(neqs), ierr_linpack

    type (bs_t) :: bs_temp

    real(kind=dp_t) :: t

    ierr = IERR_NONE

    ! substep size
    h = dt_tot/N_sub

    ! I - h J
    A(:,:) = -h * bs % burn_s % jac(:,:)
    do n = 1, neqs
       A(n,n) = ONE + A(n,n)
    enddo

    ! get the LU decomposition from LINPACK
    call dgefa(A, neqs, neqs, ipiv, ierr_linpack)
    if (ierr_linpack /= 0) then
       ierr = IERR_LU_DECOMPOSITION_ERROR
    endif

    bs_temp = bs
    bs_temp % n_rhs = 0

    ! do an Euler step to get the RHS for the first substep
    t = bs % t
    y_out(:) = h * bs % burn_s % ydot(:)

    ! solve the first step using the LU solver
    call dgesl(A, neqs, neqs, ipiv, y_out, 0)

    del(:) = y_out(:)
    bs_temp % y(:) = y(:) + del(:)

    t = t + h
    bs_temp % t = t
    call f_rhs(bs_temp)

    do n = 2, N_sub
       y_out(:) = h * bs_temp % burn_s % ydot(:) - del(:)

       ! LU solve
       call dgesl(A, neqs, neqs, ipiv, y_out, 0)

       del(:) = del(:) + TWO * y_out(:)
       bs_temp % y = bs_temp % y + del(:)

       t = t + h
       bs_temp % t = t
       call f_rhs(bs_temp)
    enddo

    y_out(:) = h * bs_temp % burn_s % ydot(:) - del(:)

    ! last LU solve
    call dgesl(A, neqs, neqs, ipiv, y_out, 0)

    ! last step
    y_out(:) = bs_temp % y(:) + y_out(:)

    ! Store the number of function evaluations.

    bs % n_rhs = bs % n_rhs + bs_temp % n_rhs

  end subroutine semi_implicit_extrap



  subroutine single_step(bs, eps, yscal, ierr)

    !$acc routine seq
    !$acc routine(jac) seq

#ifndef ACC
    use bl_error_module, only: bl_error
#endif

    implicit none

    type (bs_t) :: bs
    real(kind=dp_t), intent(in) :: eps
    real(kind=dp_t), intent(in) :: yscal(neqs)
    integer, intent(out) :: ierr

    real(kind=dp_t) :: y_save(neqs), yerr(neqs), yseq(neqs)
    real(kind=dp_t) :: err(KMAXX)

    real(kind=dp_t) :: dt, fac, scale, red, eps1, work, work_min, xest
    real(kind=dp_t) :: err_max

    logical :: converged, reduce, loop_flag

    integer :: i, k, n, kk, km, kstop, ierr_temp
    integer, parameter :: max_iters = 10 ! Should not need more than this

    ! for internal storage of the polynomial extrapolation
    real(kind=dp_t) :: t_extrap(KMAXX+1), qcol(neqs, KMAXX+1)

    ! reinitialize
    if (eps /= bs % eps_old) then
       bs % dt_next = -1.d29
       bs % t_new = -1.d29
       eps1 = S1*eps

       bs % a(1) = nseq(1)+1
       do k = 1, KMAXX
          bs % a(k+1) = bs % a(k) + nseq(k+1)
       enddo

       ! compute alpha coefficients (NR 16.4.10)
       do i = 2, KMAXX
          do k = 1, i-1
             bs % alpha(k,i) = &
                  eps1**((bs % a(k+1) - bs % a(i+1)) / &
                         ((bs % a(i+1) - bs % a(1) + ONE)*(2*k+1)))
          enddo
       enddo

       bs % eps_old = eps

       bs % a(1) = neqs + bs % a(1)
       do k = 1, KMAXX
          bs % a(k+1) = bs % a(k) + nseq(k+1)
       enddo

       ! optimal row number
       do k = 2, KMAXX-1
          if (bs % a(k+1) > bs % a(k)* bs % alpha(k-1,k)) exit
       enddo

       bs % kopt = k
       bs % kmax = k

    endif

    dt = bs % dt
    y_save(:) = bs % y(:)

    ! get the jacobian
    call jac(bs)

    if (dt /= bs % dt_next .or. bs % t /= bs % t_new) then
       bs % first = .true.
       bs % kopt = bs % kmax
    endif

    reduce = .false.

    ierr = IERR_NONE

    converged = .false.

    km = -1
    kstop = -1

    do n = 1, max_iters

       loop_flag = .false.

       do k = 1, bs % kmax

          if (.not. loop_flag) then

             bs % t_new = bs % t + dt

             call semi_implicit_extrap(bs, y_save, dt, nseq(k), yseq, ierr_temp)
             ierr = ierr_temp

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
                      loop_flag = .true.
                   else if (k == bs % kopt) then
                      if (bs % alpha(bs % kopt-1, bs % kopt) < err(km)) then
                         red = ONE/err(km)
                         reduce = .true.
                         loop_flag = .true.
                      endif
                   else if (bs % kopt == bs % kmax) then
                      if (bs % alpha(km, bs % kmax-1) < err(km)) then
                         red = bs % alpha(km, bs % kmax-1)*S2/err(km)
                         reduce = .true.
                         loop_flag = .true.
                      endif
                   else if (bs % alpha(km, bs % kopt) < err(km)) then
                      red = bs % alpha(km, bs % kopt - 1)/err(km)
                      reduce = .true.
                      loop_flag = .true.
                   endif
                endif

             endif

             kstop = k

          endif

       enddo

       if (.not. converged .and. IERR == IERR_NONE) then
          red = max(min(red, RED_BIG_FACTOR), RED_SMALL_FACTOR)
          dt = dt*red
       else
          exit
       endif

    enddo   ! while loop

    if (.not. converged) then
       IERR = IERR_NO_CONVERGENCE
    endif

#ifndef ACC
    ! km and kstop should have been set during the main loop.
    ! If they never got updated from the original nonsense values,
    ! that means something went really wrong and we should abort.

    if (km < 0) then
       call bl_error("Error: km < 0 in subroutine single_step, something has gone wrong.")
    endif

    if (kstop < 0) then
       call bl_error("Error: kstop < 0 in subroutine single_step, something has gone wrong.")
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

  end subroutine single_step


  subroutine poly_extrap(iest, test, yest, yz, dy, t, qcol)

    !$acc routine seq

    ! this does polynomial extrapolation according to the Neville
    ! algorithm.  Given test and yest (t and a y-vector), this gives
    ! the value yz and the error in the extrapolation, dy by
    ! building a polynomial through the points, where the order
    ! is iest

    integer, intent(in) :: iest
    real(kind=dp_t), intent(in) :: test, yest(neqs)
    real(kind=dp_t), intent(inout) :: yz(neqs), dy(neqs)

    ! these are for internal storage to save the state between calls
    real(kind=dp_t), intent(inout) :: t(KMAXX+1), qcol(neqs, KMAXX+1)

    integer :: j, k
    real(kind=dp_t) :: delta, f1, f2, q, d(neqs)

    t(iest) = test

    dy(:) = yest(:)
    yz(:) = yest(:)

    if (iest == 1) then
       ! nothing to do -- this is just a constant
       qcol(:,1) = yest(:)
    else
       ! we have more than 1 point, so build higher order
       ! polynomials
       d(:) = yest(:)

       do k = 1, iest-1
          delta = ONE/(t(iest-k)-test)
          f1 = test*delta
          f2 = t(iest-k)*delta

          do j = 1, neqs
             q = qcol(j,k)
             qcol(j,k) = dy(j)
             delta = d(j) - q
             dy(j) = f1*delta
             d(j) = f2*delta
             yz(j) = yz(j) + dy(j)

          enddo

       enddo
       qcol(:,iest) = dy(:)
    endif

  end subroutine poly_extrap

end module stiff_ode

