module stiff_ode

  use bl_constants_module
  use bl_types
  use burn_type_module
  use bs_type_module
  use rhs_module
  use jac_module

  implicit none

  real(kind=dp_t), parameter, private :: dt_min = 1.d-24
  real(kind=dp_t), parameter, private :: dt_ini = 1.d-16
  real(kind=dp_t), parameter, private :: SMALL = 1.d-30


  ! error codes
  integer, parameter :: IERR_NONE = 0
  integer, parameter :: IERR_DT_TOO_SMALL = -100
  integer, parameter :: IERR_TOO_MANY_STEPS = -101
  integer, parameter :: IERR_DT_UNDERFLOW = -102
  integer, parameter :: IERR_NO_CONVERGENCE = -103

  integer, parameter :: IERR_LU_DECOMPOSITION_ERROR = -200


  ! these are parameters for the BS method
  real(kind=dp_t), parameter :: S1 = 0.25_dp_t
  real(kind=dp_t), parameter :: S2 = 0.7_dp_t

  real(kind=dp_t), parameter :: RED_BIG_FACTOR = 0.7_dp_t
  real(kind=dp_t), parameter :: RED_SMALL_FACTOR = 1.e-5_dp_t
  real(kind=dp_t), parameter :: SCALMX = 0.1_dp_t

  ! these are parameters for the Rosenbrock method
  real(kind=dp_t), parameter :: GAMMA = HALF
  real(kind=dp_t), parameter :: A21 = TWO
  real(kind=dp_t), parameter :: A31 = 48.0_dp_t/25.0_dp_t
  real(kind=dp_t), parameter :: A32 = SIX/25.0_dp_t
  ! note: we are using the fact here that for both the original Kaps
  ! and Rentrop params and the parameters from Shampine (that NR
  ! likes) we have A41 = A31, A42 = A32, and A43 = 0, so the last 2
  ! solves use the same intermediate y (see Stoer & Bulirsch, TAM,
  ! p. 492)
  real(kind=dp_t), parameter :: C21 = -EIGHT
  real(kind=dp_t), parameter :: C31 = 372.0_dp_t/25.0_dp_t
  real(kind=dp_t), parameter :: C32 = TWELVE/FIVE
  real(kind=dp_t), parameter :: C41 = -112.0_dp_t/125.0_dp_t
  real(kind=dp_t), parameter :: C42 = -54.0_dp_t/125.0_dp_t
  real(kind=dp_t), parameter :: C43 = -TWO/FIVE
  real(kind=dp_t), parameter :: B1 = 19.0_dp_t/NINE
  real(kind=dp_t), parameter :: B2 = HALF
  real(kind=dp_t), parameter :: B3 = 25.0_dp_t/108.0_dp_t
  real(kind=dp_t), parameter :: B4 = 125.0_dp_t/108.0_dp_t
  real(kind=dp_t), parameter :: E1 = 17.0_dp_t/54.0_dp_t
  real(kind=dp_t), parameter :: E2 = 7.0_dp_t/36.0_dp_t
  real(kind=dp_t), parameter :: E3 = ZERO
  real(kind=dp_t), parameter :: E4 = 125.0_dp_t/108.0_dp_t
  real(kind=dp_t), parameter :: A2X = ONE
  real(kind=dp_t), parameter :: A3X = THREE/FIVE

  real(kind=dp_t), parameter :: SAFETY = 0.9_dp_t
  real(kind=dp_t), parameter :: GROW = 1.5_dp_t
  real(kind=dp_t), parameter :: PGROW = -0.25_dp_t
  real(kind=dp_t), parameter :: SHRINK = HALF
  real(kind=dp_t), parameter :: PSHRINK = -THIRD
  real(kind=dp_t), parameter :: ERRCON = 0.1296_dp_t

  integer, parameter :: MAX_TRY = 50

contains

  ! integrate from t to tmax

  subroutine safety_check(y_old, y, retry)
    !$acc routine seq

    use extern_probin_module, only: safety_factor
    
    real(kind=dp_t), intent(in) :: y_old(bs_neqs), y(bs_neqs)
    logical, intent(out) :: retry

    real(kind=dp_t) :: ratio
    integer :: n
    
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

  
  subroutine ode(bs, t, tmax, eps, ierr)

    ! this is a basic driver for the ODE integration, based on the NR
    ! routine.  This calls an integration method to take a single step
    ! and return an estimate of the net step size needed to achieve
    ! our desired tolerance.

    !$acc routine seq

    use extern_probin_module, only: ode_max_steps, use_timestep_estimator, &
                                    scaling_method, ode_scale_floor, ode_method
#ifndef ACC
    use bl_error_module, only: bl_error
#endif

    type (bs_t), intent(inout) :: bs

    real(kind=dp_t), intent(inout) :: t
    real(kind=dp_t), intent(in) :: tmax
    real(kind=dp_t), intent(in) :: eps
    integer, intent(out) :: ierr

    real(kind=dp_t) :: yscal(bs_neqs)
    logical :: finished

    integer :: n

    ! initialize

    bs % t = t
    bs % tmax = tmax

    finished = .false.
    ierr = IERR_NONE

    bs % eps_old = ZERO

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
#ifdef SDC
          yscal(:) = abs(bs % y(:)) + abs(bs % dt * bs % ydot(:)) + SMALL
#else
          yscal(:) = abs(bs % y(:)) + abs(bs % dt * bs % burn_s % ydot(:)) + SMALL
#endif

       else if (scaling_method == 2) then
          yscal = max(abs(bs % y(:)), ode_scale_floor)

#ifndef ACC
       else
          call bl_error("Unknown scaling_method in ode")
#endif
       endif

       ! make sure we don't overshoot the ending time
       if (bs % t + bs % dt > tmax) bs % dt = tmax - bs % t

       ! take a step -- this routine will update the solution vector,
       ! advance the time, and also give an estimate of the next step
       ! size
       if (ode_method == 1) then
          call single_step_bs(bs, eps, yscal, ierr)
       else if (ode_method == 2) then
          call single_step_rosen(bs, eps, yscal, ierr)
#ifndef ACC
       else
          call bl_error("Unknown ode_method in ode")
#endif
       endif

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

    ! this is a version of the timestep estimation algorithm used by
    ! VODE

    !$acc routine seq

    type (bs_t), intent(inout) :: bs

    type (bs_t) :: bs_temp
    real(kind=dp_t) :: h, h_old, hL, hU, ddydtt(bs_neqs), eps, ewt(bs_neqs), yddnorm
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
#ifdef SDC
       bs_temp % y = bs % y + h * bs % ydot
#else
       bs_temp % y = bs % y + h * bs % burn_s % ydot
#endif

       ! Call the RHS, then estimate the finite difference.
       call f_rhs(bs_temp)
#ifdef SDC
       ddydtt = (bs_temp % ydot - bs % ydot) / h
#else
       ddydtt = (bs_temp % burn_s % ydot - bs % burn_s % ydot) / h
#endif

       yddnorm = sqrt( sum( (ddydtt*ewt)**2 ) / bs_neqs )

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
    !$acc routine(dgesl) seq
    !$acc routine(dgefa) seq

    type (bs_t), intent(inout) :: bs
    real(kind=dp_t), intent(in) :: y(bs_neqs)
    real(kind=dp_t), intent(in) :: dt_tot
    integer, intent(in) :: N_sub
    real(kind=dp_t), intent(out) :: y_out(bs_neqs)
    integer, intent(out) :: ierr

    real(kind=dp_t) :: A(bs_neqs,bs_neqs)
    real(kind=dp_t) :: del(bs_neqs)
    real(kind=dp_t) :: h

    integer :: n

    integer :: ipiv(bs_neqs), ierr_linpack

    type (bs_t) :: bs_temp

    real(kind=dp_t) :: t

    ierr = IERR_NONE

    ! substep size
    h = dt_tot/N_sub

    ! I - h J
#ifdef SDC
    A(:,:) = -h * bs % jac(:,:)
#else
    A(:,:) = -h * bs % burn_s % jac(:,:)
#endif
    do n = 1, bs_neqs
       A(n,n) = ONE + A(n,n)
    enddo

    ! get the LU decomposition from LINPACK
    call dgefa(A, bs_neqs, bs_neqs, ipiv, ierr_linpack)
    if (ierr_linpack /= 0) then
       ierr = IERR_LU_DECOMPOSITION_ERROR
    endif

    if (ierr == IERR_NONE) then
       bs_temp = bs
#ifdef SDC
       bs_temp % n_rhs = 0
#else
       bs_temp % burn_s % n_rhs = 0
#endif

       ! do an Euler step to get the RHS for the first substep
       t = bs % t
#ifdef SDC
       y_out(:) = h * bs % ydot(:)
#else
       y_out(:) = h * bs % burn_s % ydot(:)
#endif

       ! solve the first step using the LU solver
       call dgesl(A, bs_neqs, bs_neqs, ipiv, y_out, 0)

       del(:) = y_out(:)
       bs_temp % y(:) = y(:) + del(:)

       t = t + h
       bs_temp % t = t
       call f_rhs(bs_temp)

       do n = 2, N_sub
#ifdef SDC
          y_out(:) = h * bs_temp % ydot(:) - del(:)
#else
          y_out(:) = h * bs_temp % burn_s % ydot(:) - del(:)
#endif

          ! LU solve
          call dgesl(A, bs_neqs, bs_neqs, ipiv, y_out, 0)

          del(:) = del(:) + TWO * y_out(:)
          bs_temp % y = bs_temp % y + del(:)

          t = t + h
          bs_temp % t = t
          call f_rhs(bs_temp)
       enddo

#ifdef SDC
       y_out(:) = h * bs_temp % ydot(:) - del(:)
#else
       y_out(:) = h * bs_temp % burn_s % ydot(:) - del(:)
#endif

       ! last LU solve
       call dgesl(A, bs_neqs, bs_neqs, ipiv, y_out, 0)

       ! last step
       y_out(:) = bs_temp % y(:) + y_out(:)
    
       ! Store the number of function evaluations.

#ifdef SDC
       bs % n_rhs = bs % n_rhs + bs_temp % n_rhs
#else
       bs % burn_s % n_rhs = bs % burn_s % n_rhs + bs_temp % burn_s % n_rhs
#endif
    else
       y_out(:) = y(:)
    endif
       
  end subroutine semi_implicit_extrap



  subroutine single_step_bs(bs, eps, yscal, ierr)

    !$acc routine seq

#ifndef ACC
    use bl_error_module, only: bl_error
#endif

    implicit none

    type (bs_t) :: bs
    real(kind=dp_t), intent(in) :: eps
    real(kind=dp_t), intent(in) :: yscal(bs_neqs)
    integer, intent(out) :: ierr

    real(kind=dp_t) :: y_save(bs_neqs), yerr(bs_neqs), yseq(bs_neqs)
    real(kind=dp_t) :: err(KMAXX)

    real(kind=dp_t) :: dt, fac, scale, red, eps1, work, work_min, xest
    real(kind=dp_t) :: err_max

    logical :: converged, reduce, skip_loop, retry

    integer :: i, k, n, kk, km, kstop, ierr_temp
    integer, parameter :: max_iters = 10 ! Should not need more than this

    ! for internal storage of the polynomial extrapolation
    real(kind=dp_t) :: t_extrap(KMAXX+1), qcol(bs_neqs, KMAXX+1)

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

       bs % a(1) = bs_neqs + bs % a(1)
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

    converged = .false.

    km = -1
    kstop = -1

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
       call bl_error("Error: non convergence in single_step_bs, something has gone wrong.")
#endif       
    endif

#ifndef ACC
    ! km and kstop should have been set during the main loop.
    ! If they never got updated from the original nonsense values,
    ! that means something went really wrong and we should abort.

    if (km < 0) then
       call bl_error("Error: km < 0 in subroutine single_step_bs, something has gone wrong.")
    endif

    if (kstop < 0) then
       call bl_error("Error: kstop < 0 in subroutine single_step_bs, something has gone wrong.")
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


  subroutine poly_extrap(iest, test, yest, yz, dy, t, qcol)

    !$acc routine seq

    ! this does polynomial extrapolation according to the Neville
    ! algorithm.  Given test and yest (t and a y-vector), this gives
    ! the value yz and the error in the extrapolation, dy by
    ! building a polynomial through the points, where the order
    ! is iest

    integer, intent(in) :: iest
    real(kind=dp_t), intent(in) :: test, yest(bs_neqs)
    real(kind=dp_t), intent(inout) :: yz(bs_neqs), dy(bs_neqs)

    ! these are for internal storage to save the state between calls
    real(kind=dp_t), intent(inout) :: t(KMAXX+1), qcol(bs_neqs, KMAXX+1)

    integer :: j, k
    real(kind=dp_t) :: delta, f1, f2, q, d(bs_neqs)

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

          do j = 1, bs_neqs
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


  subroutine single_step_rosen(bs, eps, yscal, ierr)

    ! this does a single step of the Rosenbrock method.  Note: we
    ! assume here that our RHS is not an explicit function of t, but
    ! only of our integration variable, y

    !$acc routine seq
    !$acc routine(dgesl) seq
    !$acc routine(dgefa) seq

#ifndef ACC
    use bl_error_module, only: bl_error
#endif

    implicit none

    type (bs_t) :: bs
    real(kind=dp_t), intent(in) :: eps
    real(kind=dp_t), intent(in) :: yscal(bs_neqs)
    integer, intent(out) :: ierr
    
    real(kind=dp_t) :: A(bs_neqs,bs_neqs)
    real(kind=dp_t) :: g1(bs_neqs), g2(bs_neqs), g3(bs_neqs), g4(bs_neqs)
    real(kind=dp_t) :: err(bs_neqs)

    real(kind=dp_t) :: h, h_tmp, errmax

    integer :: q, n

    integer :: ipiv(bs_neqs), ierr_linpack

    type (bs_t) :: bs_temp

    logical :: converged

    h = bs % dt

    ! note: we come in already with a RHS evalulation from the driver

    ! get the jacobian
    call jac(bs)

    ierr = IERR_NONE

    converged = .false.

    q = 1
    do while (q <= MAX_TRY .and. .not. converged .and. ierr == IERR_NONE)

       bs_temp = bs

       ! create I/(gamma h) - ydot -- this is the matrix used for all the
       ! linear systems that comprise a single step
#ifdef SDC
       A(:,:) = -bs % jac(:,:)
#else
       A(:,:) = -bs % burn_s % jac(:,:)
#endif
       do n = 1, bs_neqs
          A(n,n) = ONE/(gamma * h) + A(n,n)
       enddo
       
       ! LU decomposition
       call dgefa(A, bs_neqs, bs_neqs, ipiv, ierr_linpack)
       if (ierr_linpack /= 0) then
          ierr = IERR_LU_DECOMPOSITION_ERROR
       endif
       
       ! setup the first RHS and solve the linear system (note: the linear
       ! solve replaces the RHS with the solution in place)
#ifdef SDC
       g1(:) = bs % ydot(:)
#else
       g1(:) = bs % burn_s % ydot(:)
#endif

       call dgesl(A, bs_neqs, bs_neqs, ipiv, g1, 0)

       ! new value of y
       bs_temp % y(:) = bs % y(:) + A21*g1(:)
       bs_temp % t = bs % t + A2X*h
       
       ! get derivatives at this intermediate position and setup the next
       ! RHS
       call f_rhs(bs_temp)

#ifdef SDC
       g2(:) = bs_temp % ydot(:) + C21*g1(:)/h
#else
       g2(:) = bs_temp % burn_s % ydot(:) + C21*g1(:)/h
#endif       
       call dgesl(A, bs_neqs, bs_neqs, ipiv, g2, 0)

       ! new value of y
       bs_temp % y(:) = bs % y(:) + A31*g1(:) + A32*g2(:)
       bs_temp % t = bs % t + A3X*h

       ! get derivatives at this intermediate position and setup the next
       ! RHS
       call f_rhs(bs_temp)

#ifdef SDC
       g3(:) = bs_temp % ydot(:) + (C31*g1(:) + C32*g2(:))/h
#else
       g3(:) = bs_temp % burn_s % ydot(:) + (C31*g1(:) + C32*g2(:))/h
#endif
       
       call dgesl(A, bs_neqs, bs_neqs, ipiv, g3, 0)

       ! our choice of parameters prevents us from needing another RHS 
       ! evaluation here

       ! final intermediate RHS
#ifdef SDC
       g4(:) = bs_temp % ydot(:) + (C41*g1(:) + C42*g2(:) + C43*g3(:))/h
#else
       g4(:) = bs_temp % burn_s % ydot(:) + (C41*g1(:) + C42*g2(:) + C43*g3(:))/h
#endif
       
       call dgesl(A, bs_neqs, bs_neqs, ipiv, g4, 0)

       ! now construct our 4th order estimate of y
       bs_temp % y(:) = bs % y(:) + B1*g1(:) + B2*g2(:) + B3*g3(:) + B4*g4(:)
       bs_temp % t = bs % t + h
       err(:) = E1*g1(:) + E2*g2(:) + E3*g3(:) + E4*g4(:)

       if (bs_temp % t == bs % t) then
          ierr = IERR_DT_UNDERFLOW
       endif

       ! get the error and scale it to the desired tolerance
       errmax = maxval(abs(err(:)/yscal(:)))
       errmax = errmax/eps

       if (errmax <= 1) then
          ! we were successful -- store the solution
          bs % y(:) = bs_temp % y(:)
          bs % t = bs_temp % t

          bs % dt_did = h
          if (errmax > ERRCON) then
             bs % dt_next = SAFETY*h*errmax**PGROW
          else
             bs % dt_next = GROW*h
          endif
          
          converged = .true.

       else if (ierr == IERR_NONE) then
          ! integration did not meet error criteria.  Return h and
          ! try again

          ! this is essentially the step control from Stoer &
          ! Bulircsh (TAM) Eq. 7.2.5.17, as shown on p. 493
          h_tmp = SAFETY*h*errmax**PSHRINK

          h = sign(max(abs(h_tmp), SHRINK*abs(h)), h)
       endif
          
       q = q + 1
       
    enddo
    
    if (.not. converged .and. ierr == IERR_NONE) then
       ierr = IERR_NO_CONVERGENCE
    endif

  end subroutine single_step_rosen

end module stiff_ode

