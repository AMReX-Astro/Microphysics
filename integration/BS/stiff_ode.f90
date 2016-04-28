! A variable-sized stepsize driver for the BS integrator,
! based on the NR routine


module stiff_ode

  use bl_types

  implicit none

  real(kind=dp_t), parameter :: ZERO = 0.0_dp_t
  real(kind=dp_t), parameter :: ONE = 1.0_dp_t
  real(kind=dp_t), parameter :: TWO = 2.0_dp_t


  real(kind=dp_t), parameter :: dt_min = 1.d-20
  real(kind=dp_t), parameter :: dt_ini = 1.d-10
  real(kind=dp_t), parameter :: SMALL = 1.d-30

  integer, parameter :: MAX_STEPS = 10000

  ! BS parameters -- see the discussion in 16.4
  integer, parameter :: KMAXX = 7
  integer, parameter, dimension(KMAXX+1) :: nseq(KMAXX+1) = [2, 6, 10, 14, 22, 34, 50, 70]

  ! error codes
  integer, parameter :: IERR_NONE = 0
  integer, parameter :: IERR_DT_TOO_SMALL = -100
  integer, parameter :: IERR_TOO_MANY_STEPS = -101
  integer, parameter :: IERR_DT_UNDERFLOW = -102

  integer, parameter :: IERR_LU_DECOMPOSITION_ERROR = -200

  real(kind=dp_t), parameter :: S1 = 0.25_dp_t
  real(kind=dp_t), parameter :: S2 = 0.7_dp_t

  real(kind=dp_t), parameter :: REDMIN = 0.7_dp_t
  real(kind=dp_t), parameter :: REDMAX = 1.e-5_dp_t
  real(kind=dp_t), parameter :: SCALMX = 0.1_dp_t

  type integrator_t
     logical :: first
     real(kind=dp_t) :: eps_old
     real(kind=dp_t) :: dt_did
     real(kind=dp_t) :: dt_next
     real(kind=dp_t) :: a(KMAXX+1)
     real(kind=dp_t) :: alpha(KMAXX, KMAXX)
     real(kind=dp_t) :: t_new
     integer :: kmax
     integer :: kopt
  end type integrator_t

contains

  subroutine ode(yinit, neq, t, tmax, eps, f_rhs, ierr)
    ! integrate from t to tmax

    integer, intent(in) :: neq
    real(kind=dp_t), intent(inout) :: yinit(neq)
    real(kind=dp_t), intent(inout) :: t
    real(kind=dp_t), intent(in) :: tmax
    real(kind=dp_t), intent(in) :: eps
    integer, intent(out) :: ierr

    external f_rhs

    real(kind=dp_t) :: y(neq), yscal(neq), dydt(neq)
    real(kind=dp_t) :: dt, dt_next
    logical :: finished

    type (integrator_t) :: int_stat

    integer :: n

    ! initialize
    y(:) = yinit(:)
    dt = dt_ini

    finished = .false.
    ierr = IERR_NONE

    do n = 1, MAX_STEPS

       ! get the scaling
       call f_rhs(neq, t, y, dydt)
       yscal(:) =  abs(y(:)) + abs(dt*dydt(:)) + SMALL

       ! make sure we don't overshoot the ending time
       if (t + dt > tmax) dt = tmax - t

       ! take a step -- this routine will update the solution vector,
       ! advance the time, and also give an estimate of the next step
       ! size
       call single_step(y, dydt, neq, t, dt, eps, yscal, int_stat, f_rhs, ierr)

       if (.not. ierr == IERR_NONE) exit

       ! finished?
       if (t - tmax >= ZERO) then
          yinit(:) = y(:)
          finished = .true.
          exit
       endif

       dt = int_stat % dt_next

       if (dt < dt_min) then
          ierr = IERR_DT_TOO_SMALL
          exit
       endif

    enddo

    if (.not. finished .and. ierr == IERR_NONE) then
       ierr = IERR_TOO_MANY_STEPS
    endif

  end subroutine ode


  subroutine semi_implicit_extrap(y, dydt, J, neq, t0, dt_tot, N_sub, y_out, f_rhs)

    integer, intent(in) :: neq
    real(kind=dp_t), intent(in) :: y(neq), dydt(neq), J(neq, neq)
    real(kind=dp_t), intent(in) :: t0, dt_tot
    integer, intent(in) :: N_sub
    real(kind=dp_t), intent(out) :: y_out(neq)

    external f_rhs

    real(kind=dp_t) :: A(neq,neq)
    real(kind=dp_t) :: del(neq), y_temp(neq), dydt_h(neq)
    real(kind=dp_t) :: h

    integer :: n

    integer :: ipiv(neq), ierr_linpack

    real(kind=dp_t) :: t

    ! substep size
    h = dt_tot/N_sub

    ! I - h J
    A(:,:) = -h*J(:,:)
    do n = 1, neq
       A(n,n) = ONE + A(n,n)
    enddo

    ! get the LU decomposition from LIPACK
    call dgefa(A, neq, neq, ipiv, ierr_linpack)

    ! do an Euler step to get the RHS for the first substep
    t = t0
    y_out(:) = h*dydt(:)

    ! solve the first step using the LU solver
    call dgesl(A, neq, neq, ipiv, y_out, 0)

    del(:) = y_out(:)
    y_temp(:) = y(:) + del(:)

    t = t + h
    call f_rhs(neq, t, y_temp, dydt_h)

    do n = 2, N_sub
       y_out(:) = h*dydt_h(:) - del(:)

       ! LU solve
       call dgesl(A, neq, neq, ipiv, y_out, 0)

       del(:) = del(:) + TWO*y_out(:)
       y_temp(:) = y_temp(:) + del(:)

       t = t + h
       call f_rhs(neq, t, y_temp, dydt_h)
    enddo

    y_out(:) = h*dydt_h(:) - del(:)

    ! last LU solve
    call dgesl(A, neq, neq, ipiv, y_out, 0)

    ! last step
    y_out(:) = y_temp(:) + y_out(:)

  end subroutine semi_implicit_extrap


  subroutine single_step(y, dydt, neq, t, dt_try, eps, yscal, int_stat, f_rhs, ierr)

    integer, intent(in) :: neq
    real(kind=dp_t), intent(inout) :: y(neq)
    real(kind=dp_t), intent(out) :: dydt(neq)
    real(kind=dp_t), intent(inout) :: t
    real(kind=dp_t), intent(in) :: dt_try
    real(kind=dp_t), intent(in) :: eps
    real(kind=dp_t), intent(in) :: yscal(neq)
    type(integrator_t), intent(inout) :: int_stat
    integer, intent(out) :: ierr

    external f_rhs

    real(kind=dp_t) :: y_save(neq), yerr(neq), yseq(neq)
    real(kind=dp_t) :: err(KMAXX)

    real(kind=dp_t) :: dfdy(neq, neq)

    real(kind=dp_t) :: dt, fac, scale, red, eps1, work, work_min, xest
    real(kind=dp_t) :: err_max

    logical :: converged, reduce

    integer :: i, k, kk, km, kopt

    ! for internal storage of the polynomial extrapolation
    real(kind=dp_t) :: t_extrap(KMAXX+1), qcol(neq, KMAXX+1)

    ! reinitialize
    if (eps /= int_stat % eps_old) then
       int_stat % dt_next = -1.d29
       int_stat % t_new = -1.d29
       eps1 = S1*eps

       int_stat % a(1) = nseq(1)+1
       do k = 1, KMAXX
          int_stat % a(k+1) = int_stat % a(k) + nseq(k+1)
       enddo

       ! compute alpha coefficients (NR 16.4.10)
       do i = 2, KMAXX
          do k = 1, i-1
             int_stat % alpha(k,i) = &
                  eps1**((int_stat % a(k+1) - int_stat % a(i+1)) / &
                         ((int_stat % a(i+1) - int_stat % a(1) + ONE)*(2*k+1)))
          enddo
       enddo

       int_stat % eps_old = eps

       int_stat % a(1) = neq + int_stat % a(1)
       do k = 1, KMAXX
          int_stat % a(k+1) = int_stat % a(k) + nseq(k+1)
       enddo


       ! optimal row number
       do kopt = 2, KMAXX-1
          if (int_stat % a(kopt+1) > int_stat % a(kopt)* int_stat % alpha(kopt-1,kopt)) exit
       enddo
       int_stat % kopt = kopt
       int_stat % kmax = int_stat % kopt


    endif

    dt = dt_try
    y_save(:) = y(:)

    ! get the jacobian
    call jac(neq, t, y, dfdy)

    if (dt /= int_stat % dt_next .or. t /= int_stat % t_new) then
       int_stat % first = .true.
       int_stat % kopt = int_stat % kmax
    endif

    reduce = .false.

    ierr = IERR_NONE

    converged = .false.

    do while (.not. converged .and. ierr == IERR_NONE)

       do k = 1, int_stat % kmax

          int_stat % t_new = t + dt
          if (int_stat % t_new == t) then
             ierr = IERR_DT_UNDERFLOW
             exit
          endif

          call semi_implicit_extrap(y_save, dydt, dfdy, neq, t, dt, nseq(k), yseq, f_rhs)

          xest = (dt/nseq(k))**2

          call poly_extrap(k, xest, yseq, y, yerr, neq, t_extrap, qcol)

          if (k /= 1) then
             err_max = max(SMALL, maxval(abs(yerr(:)/yscal(:))))
             err_max = err_max / eps
             km = k - 1
             err(km) = (err_max/S1)**(1.0/(2*km+1))
          endif

          if (k /= 1 .and. (k >=  int_stat % kopt-1 .or. int_stat % first)) then
             if (err_max < 1) then
                converged = .true.
                exit
             endif

             ! reduce stepsize if necessary
             if (k == int_stat % kmax .or. k == int_stat % kopt+1) then
                red = S2/err(km)
                reduce = .true.
                exit
             else if (k == int_stat % kopt) then
                if (int_stat % alpha(int_stat % kopt-1, int_stat % kopt-1) < err(km)) then
                   red = ONE/err(km)
                   reduce = .true.
                   exit
                endif
             else if (int_stat % kopt == int_stat % kmax) then
                if (int_stat % alpha(km, int_stat % kmax-1) < err(km)) then
                   red = int_stat % alpha(km, int_stat % kmax-1)*S2/err(km)
                   reduce = .true.
                   exit
                endif
             else if (int_stat % alpha(km, int_stat % kopt) < err(km)) then
                red = int_stat % alpha(km, int_stat % kopt - 1)/err(km)
                reduce = .true.
                exit
             endif

          endif

       enddo

       if (.not. converged .and. IERR == IERR_NONE) then
          red = max(min(red, REDMIN), REDMAX)
          dt = dt*red
       else
          ! we've either converged or hit and error
          exit
       endif
    enddo   ! while loop

    t = int_stat % t_new
    int_stat % dt_did = dt
    int_stat % first = .false.

    ! optimal convergence properties
    work_min = 1.e35
    do kk = 1, km
       fac = max(err(kk), SCALMX)
       work = fac*int_stat % a(kk+1)

       if (work < work_min) then
          scale = fac
          work_min = work
          int_stat % kopt = kk+1
       endif
    enddo

    ! increase in order
    int_stat % dt_next = dt / scale

    if (int_stat % kopt >= k .and. int_stat % kopt /= int_stat % kmax .and. .not. reduce) then
       fac = max(scale/int_stat % alpha(int_stat % kopt-1, int_stat % kopt), SCALMX)
       if (int_stat % a(int_stat % kopt+1)*fac <= work_min) then
          int_stat % dt_next = dt/fac
          int_stat % kopt = int_stat % kopt + 1
       endif
    endif

  end subroutine single_step


  subroutine poly_extrap(iest, test, yest, yz, dy, neq, t, qcol)

    integer, intent(in) :: iest, neq
    real(kind=dp_t), intent(in) :: test, yest(neq)
    real(kind=dp_t), intent(inout) :: yz(neq), dy(neq)

    ! these are for internal storage to save the state between calls
    real(kind=dp_t), intent(inout) :: t(KMAXX+1), qcol(neq, KMAXX+1)

    integer :: j, k
    real(kind=dp_t) :: delta, f1, f2, q, d(neq)

    t(iest) = test

    dy(:) = yest(:)
    yz(:) = yest(:)

    if (iest == 1) then
       qcol(:,1) = yest(:)
    else
       d(:) = yest(:)
       do k = 1, iest-1
          delta = ONE/(t(iest-k)-test)
          f1 = test*delta
          f2 = t(iest-k)*delta
          do j = 1, neq
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

