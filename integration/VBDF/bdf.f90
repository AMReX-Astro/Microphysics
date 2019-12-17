! BDF (backward differentiation formula) time-stepping routines.
!
! See
!
!   1. VODE: A variable-coefficient ODE solver; Brown, Byrne, and
!      Hindmarsh; SIAM J. Sci. Stat. Comput., vol. 10, no. 5, pp.
!      1035-1051, 1989.
!
!   2. An alternative implementation of variable step-size multistep
!      formulas for stiff ODES; Jackson and Sacks-Davis; ACM
!      Trans. Math. Soft., vol. 6, no. 3, pp. 295-318, 1980.
!
!   3. A polyalgorithm for the numerical solution of ODEs; Byrne and
!      Hindmarsh; ACM Trans. Math. Soft., vol. 1, no. 1, pp. 71-96,
!      1975.
!



module bdf

  use amrex_constants_module
  use amrex_fort_module, only : rt => amrex_real
  use bdf_type_module

  implicit none

  !TODO: Lowered iters for dev, change back
  integer, parameter :: bdf_max_iters = 666666666

  integer, parameter :: BDF_ERR_SUCCESS  = 0
  integer, parameter :: BDF_ERR_SOLVER   = 1
  integer, parameter :: BDF_ERR_MAXSTEPS = 2
  integer, parameter :: BDF_ERR_DTMIN    = 3

  character(len=64), parameter :: errors(0:3) = [ &
       'Success.                                                ', &
       'Newton solver failed to converge several times in a row.', &
       'Too many steps were taken.                              ', &
       'Minimum time-step reached several times in a row.       ' ]

  private :: &
       rescale_timestep, decrease_order, increase_order, &
       alpha0, alphahat0, xi_j, xi_star_inv, ewts, norm, eye_r, eye_i, &
       factorial, eoshift_local
  !public subroutines: bdf_advance, bdf_update, bdf_predict, bdf_solve, bdf_check
  !                    bdf_correct, bdf_dump, bdf_adjust, bdf_reset, print_y
  !                    bdf_ts_build, bdf_ts_destroy, bdf_wrap

contains

  !
  ! Advance system from t0 to t1.
  !
  subroutine bdf_advance(ts, y0, t0, y1, t1, dt0, reset, reuse, ierr, initial_call)
    !$acc routine seq
    type(bdf_ts), intent(inout) :: ts
    real(rt),   intent(in   ) :: y0(neqs,bdf_npt), t0, t1, dt0
    real(rt),   intent(  out) :: y1(neqs,bdf_npt)
    logical,      intent(in   ) :: reset, reuse
    integer,      intent(  out) :: ierr
    logical,      intent(in   ) :: initial_call
    integer  :: k, p, m, n
    logical  :: retry, linitial

    !TODO: We no longer have this argument as optional, so rewrite to get rid of linitial,
    !or maybe just get rid of it.  Commented out for now.  I prefer to use this,
    !but for GPU dev I'm trying to simplify.
    !linitial = initial_call

    if (reset) then
       ts%t = t0
       call bdf_reset(ts, y0, dt0, reuse)
    endif

    ierr = BDF_ERR_SUCCESS

    !print *, 'y0: ', y0
    ts%t1 = t1; ts%t = t0; ts%ncse = 0; ts%ncdtmin = 0;
    do k = 1, bdf_max_iters + 1
       if (ts%n > ts%max_steps .or. k > bdf_max_iters) then
          !ierr = BDF_ERR_MAXSTEPS; return
          ierr = BDF_ERR_MAXSTEPS; exit
       end if

       !TODO: Debug I/O not cool on GPUs. If we want to keep it, need to rewrite
       !if (k == 1) &
       !     call bdf_dump(ts)

       call bdf_update(ts)                ! update various coeffs (l, tq) based on time-step history
       call bdf_predict(ts)               ! predict nordsieck array using pascal matrix
       !if(linitial .and. k == 1) then
       !   !This is the initial solve, so use the user's initial value, 
       !   !not the predicted value.
       !   do p = 1, ts%npt
       !      do m = 1, ts%neq
       !         !Overwrite the predicted z0 with the user's y0
       !         ts%z0(m,p,0) = ts%y(m,p)
       !      end do
       !   end do
       !endif
       !print *, 'b4: ', ts%z0(1,1,0)
       call bdf_solve(ts)         ! solve for y_n based on predicted y and yd
       !print *, 'af: ', ts%z0(1,1,0)
       call bdf_check(ts, retry, ierr)    ! check for solver errors and test error estimate

       !if (ierr /= BDF_ERR_SUCCESS) return
       if (ierr /= BDF_ERR_SUCCESS) exit
       !TODO: cycle statements may lead to bad use of coalesced memory in OpenACC (or busy waiting),
       !look into this when tuning
       if (retry) cycle

       call bdf_correct(ts)               ! new solution looks good, correct history and advance

       !call bdf_dump(ts)
       !TODO: exit statements may lead to bad use of coalesced memory in OpenACC (or busy waiting),
       !look into this when tuning
       if (ts%t >= t1) exit

       call bdf_adjust(ts)                ! adjust step-size/order
    end do

    !TODO: GPUs don't like print statements.  Either delete this or work up alternative implementations
    !if (ts%verbose > 0) &
    !     print '("BDF: n:",i6,", fe:",i6,", je: ",i3,", lu: ",i3,", it: ",i3,", se: ",i3,", dt: ",e15.8,", k: ",i2)', &
    !     ts%n, ts%nfe, ts%nje, ts%nlu, ts%nit, ts%nse, ts%dt, ts%k

    !y1 = ts%z(:,:,0)
    do p = 1, ts%npt
       do m = 1, ts%neq
          y1(m,p) = ts%z(m,p,0)
       end do
    end do
    
    !print *, 'y1: ', y1
  end subroutine bdf_advance

  !
  ! Compute Nordsieck update coefficients l and error coefficients tq.
  !
  ! Regarding the l coefficients, see section 5, and in particular
  ! eqn. 5.2, of Jackson and Sacks-Davis (1980).
  !
  ! Regarding the error coefficients tq, these have been adapted from
  ! cvode.  The tq array is indexed as:
  !
  !  tq(-1) coeff. for order k-1 error est.
  !  tq(0)  coeff. for order k error est.
  !  tq(1)  coeff. for order k+1 error est.
  !  tq(2)  coeff. for order k+1 error est. (used for e_{n-1})
  !
  ! Note:
  !
  !   1. The input vector t = [ t_n, t_{n-1}, ... t_{n-k} ] where we
  !      are advancing from step n-1 to step n.
  !
  !   2. The step size h_n = t_n - t_{n-1}.
  !
  subroutine bdf_update(ts)
    !$acc routine seq
    type(bdf_ts), intent(inout) :: ts

    integer  :: j, o
    real(rt) :: a0, a0hat, a1, a2, a3, a4, a5, a6, xistar_inv, xi_inv, c

    !ts%l  = 0
    do o = 0, ts%max_order
       ts%l(o) = 0
    end do
    ts%tq = 0

    ! compute l vector
    ts%l(0) = 1
    ts%l(1) = xi_j(ts%h, 1)
    if (ts%k > 1) then
       do j = 2, ts%k-1
          !NOTE: this is causing a conformable error, had to replace with
          !explicit loop
          !  ts%l = ts%l + eoshift_local(ts%l, -1) / xi_j(ts%h, j)
          !l_shift = eoshift_local(ts%l, -1)
          call eoshift_local(ts%l, -1, ts%shift)
          do o = 0, ts%max_order
             ts%l(o) = ts%l(o) + ts%shift(o) / xi_j(ts%h, j)
          end do
       end do
       !l_shift = eoshift_local(ts%l, -1)
       call eoshift_local(ts%l, -1, ts%shift)
       do o = 0, ts%max_order
          ts%l(o) = ts%l(o) + ts%shift(o) * xi_star_inv(ts%k, ts%h)
       end do
    end if

    ! compute error coefficients (adapted from cvode)
    a0hat = alphahat0(ts%k, ts%h)
    a0    = alpha0(ts%k)

    xi_inv     = ONE
    xistar_inv = ONE
    if (ts%k > 1) then
       xi_inv     = ONE / xi_j(ts%h, ts%k)
       xistar_inv = xi_star_inv(ts%k, ts%h)
    end if

    a1 = ONE - a0hat + a0
    a2 = ONE + ts%k * a1
    ts%tq(0) = abs(a1 / (a0 * a2))
    ts%tq(2) = abs(a2 * xistar_inv / (ts%l(ts%k) * xi_inv))
    if (ts%k > 1) then
       c  = xistar_inv / ts%l(ts%k)
       a3 = a0 + ONE / ts%k
       a4 = a0hat + xi_inv
       ts%tq(-1) = abs(c * (ONE - a4 + a3) / a3)
    else
       ts%tq(-1) = ONE
    end if

    xi_inv = ts%h(0) / sum(ts%h(0:ts%k))
    a5 = a0 - ONE / (ts%k+1)
    a6 = a0hat - xi_inv
    ts%tq(1) = abs((ONE - a6 + a5) / a2 / (xi_inv * (ts%k+2) * a5))

    call ewts(ts)
  end subroutine bdf_update

  !
  ! Predict (apply Pascal matrix).
  !
  subroutine bdf_predict(ts)
    !$acc routine seq
    type(bdf_ts), intent(inout) :: ts
    integer :: i, j, m, p
    do i = 0, ts%k
       do p = 1, ts%npt
          !ts%z0(:,p,i) = 0
          do m = 1, ts%neq
             ts%z0(m,p,i) = 0
          end do
          do j = i, ts%k
             do m = 1, ts%neq
                ts%z0(m,p,i) = ts%z0(m,p,i) + A(i,j) * ts%z(m,p,j)
             end do
          end do
       end do
    end do
  end subroutine bdf_predict

  !
  ! Solve "y_n - dt f(y_n,t) = y - dt yd" for y_n where y and yd are
  ! predictors from the Nordsieck form.
  !
  ! Newton iteration is:
  !   solve:   P x = -c G(y(k)) for x
  !   update:  y(k+1) = y(k) + x
  ! where
  !   G(y) = y - dt * f(y,t) - rhs
  !
  subroutine bdf_solve(ts)
    !$acc routine seq
    !$acc routine(dgefa) seq
    !$acc routine(dgesl) seq

    use rhs_module, only: rhs, jac

    type(bdf_ts), intent(inout) :: ts

    integer  :: k, m, n, p, info
    real(rt) :: c, dt_adj, dt_rat, inv_l1
    logical  :: rebuild, iterating(ts%npt)

    inv_l1 = 1.0_rt / ts%l(1)
    do p = 1, ts%npt
       do m = 1, ts%neq
          ts%e(m,p)   = 0
          ts%rhs(m,p) = ts%z0(m,p,0) - ts%z0(m,p,1) * inv_l1
          ts%y(m,p)   = ts%z0(m,p,0)
       end do
    end do
    dt_adj = ts%dt / ts%l(1)

    dt_rat = dt_adj / ts%dt_nwt
    if (ts%p_age > ts%max_p_age) ts%refactor = .true.
    if (dt_rat < 0.7e0_rt .or. dt_rat > 1.429e0_rt) ts%refactor = .true.

    iterating = .true.

    do k = 1, ts%max_iters

       ! build iteration matrix and factor
       if (ts%refactor) then
          rebuild = .true.
          if (ts%ncse == 0 .and. ts%j_age < ts%max_j_age) rebuild = .false.
          if (ts%ncse > 0  .and. (dt_rat < 0.2e0_rt .or. dt_rat > 5.e0_rt)) rebuild = .false.

          if (rebuild) then
             call jac(ts)
             ts%nje   = ts%nje + 1*ts%npt
             ts%j_age = 0
          end if

          call eye_r(ts%P)

          do p =1, ts%npt
             do m = 1, ts%neq
                do n = 1, ts%neq
                   ts%P(n,m,p) = ts%P(n,m,p) - dt_adj * ts%J(n,m,p)
                enddo
             enddo

             call dgefa(ts%P(:,:,p), ts%neq, ts%neq, ts%ipvt(:,p), info)
             if (info /= 0) then
                print *, "error: singular matrix"
             endif

             ! lapack      call dgetrf(neq, neq, ts%P, neq, ts%ipvt, info)
             ts%nlu    = ts%nlu + 1
          end do

          ts%dt_nwt = dt_adj
          ts%p_age  = 0
          ts%refactor  = .false.
       end if

       c = 2 * ts%dt_nwt / (dt_adj + ts%dt_nwt)

       call rhs(ts)
       ts%nfe = ts%nfe + 1

       do p = 1, ts%npt
          if (.not. iterating(p)) cycle

          ! solve using factorized iteration matrix
          do m = 1, ts%neq
             ts%b(m,p) = c * (ts%rhs(m,p) - ts%y(m,p) + dt_adj * ts%yd(m,p))
          end do
          call dgesl(ts%P(:,:,p), ts%neq, ts%neq, ts%ipvt(:,p), ts%b(:,p), 0)
          ! lapack   call dgetrs ('N', neq, 1, ts%P, neq, ts%ipvt, ts%b, neq, info)
          ts%nit = ts%nit + 1

          do m = 1, ts%neq
             ts%e(m,p) = ts%e(m,p) + ts%b(m,p)
             ts%y(m,p) = ts%z0(m,p,0) + ts%e(m,p)
          end do
          if (norm(ts%b(:,p), ts%ewt(:,p)) < ONE) iterating(p) = .false.
       end do

       if (.not. any(iterating)) exit

    end do

    ts%ncit = k; ts%p_age = ts%p_age + 1; ts%j_age = ts%j_age + 1
  end subroutine bdf_solve

  !
  ! Check error estimates.
  !
  subroutine bdf_check(ts, retry, err)
    !$acc routine seq
    type(bdf_ts), intent(inout) :: ts
    logical,      intent(out)   :: retry
    integer,      intent(out)   :: err

    real(rt) :: error, eta
    integer  :: p

    retry = .false.; err = BDF_ERR_SUCCESS

    ! if solver failed many times, bail
    if (ts%ncit >= ts%max_iters .and. ts%ncse > 7) then
       err = BDF_ERR_SOLVER
       return
    end if

    ! if solver failed to converge, shrink dt and try again
    if (ts%ncit >= ts%max_iters) then
       ts%refactor = .true.; ts%nse = ts%nse + 1; ts%ncse = ts%ncse + 1
       call rescale_timestep(ts, 0.25e0_rt, .false.)
       retry = .true.
       return
    end if
    ts%ncse = 0

    ! if local error is too large, shrink dt and try again
    do p = 1, ts%npt
       error = ts%tq(0) * norm(ts%e(:,p), ts%ewt(:,p))
       if (error > ONE) then
          eta = ONE / ( (6.e0_rt * error) ** (ONE / ts%k) + 1.e-6_rt )
          call rescale_timestep(ts, eta, .false.)
          retry = .true.
          !if (ts%dt < ts%dt_min + epsilon(ts%dt_min)) ts%ncdtmin = ts%ncdtmin + 1
          if (ts%dt < ts%dt_min) ts%ncdtmin = ts%ncdtmin + 1
          if (ts%ncdtmin > 7) err = BDF_ERR_DTMIN
          return
       end if
    end do
    ts%ncdtmin = 0

  end subroutine bdf_check

  !
  ! Correct (apply l coeffs) and advance step.
  !
  subroutine bdf_correct(ts)
    !$acc routine seq
    type(bdf_ts), intent(inout) :: ts
    integer :: i, m, p, o

    do i = 0, ts%k
       do p = 1, ts%npt
          do m = 1, ts%neq
             ts%z(m,p,i) = ts%z0(m,p,i) + ts%e(m,p) * ts%l(i)
          end do
       end do
    end do

    !ts%h     = eoshift_local(ts%h, -1)
    !h_shift = eoshift_local(ts%h, -1)
    call eoshift_local(ts%h, -1, ts%shift)
    do o = 0, ts%max_order
       ts%h(o) = ts%shift(o)
    end do
    ts%h(0)  = ts%dt
    ts%t     = ts%t + ts%dt
    ts%n     = ts%n + 1
    ts%k_age = ts%k_age + 1
  end subroutine bdf_correct


  !
  ! Dump (for debugging)...
  !
  subroutine bdf_dump(ts)
    type(bdf_ts), intent(inout) :: ts
    integer :: i, m, p

    if (.not. ts%debug) return
    write(ts%dump_unit,*) ts%t, ts%z(:,:,0)
  end subroutine bdf_dump

  !
  ! Adjust step-size/order to maximize step-size.
  !
  subroutine bdf_adjust(ts)
    !$acc routine seq
    type(bdf_ts), intent(inout) :: ts

    real(rt) :: c, error, eta(-1:1), rescale, etamax(ts%npt), etaminmax
    real(rt) :: cxe1(ts%neq, ts%npt), emcxe1(ts%neq, ts%npt), delta(ts%npt)
    integer  :: p, m

    rescale = 0

    do p = 1, ts%npt
       ! compute eta(k-1), eta(k), eta(k+1)
       eta = 0
       error  = ts%tq(0) * norm(ts%e(:,p), ts%ewt(:,p))
       eta(0) = ONE / ( (6.e0_rt * error) ** (ONE / ts%k) + 1.e-6_rt )
       if (ts%k_age > ts%k) then
          if (ts%k > 1) then
             error     = ts%tq(-1) * norm(ts%z(:,p,ts%k), ts%ewt(:,p))
             eta(-1) = ONE / ( (6.e0_rt * error) ** (ONE / ts%k) + 1.e-6_rt )
          end if
          if (ts%k < ts%max_order) then
             c = (ts%tq(2) / ts%tq2save) * (ts%h(0) / ts%h(2)) ** (ts%k+1)
             !error  = ts%tq(1) * norm(ts%e(:,p) - c * ts%e1(:,p), ts%ewt(:,p))
             do m = 1, ts%neq
                !NOTE: we have to calculate these temporary arrays because
                !   the original code required an implicit allocation which is not
                !   allowed on GPUs
                cxe1(m,p) = c * ts%e1(m,p)
                emcxe1(m,p) = ts%e(m,p) - cxe1(m,p)
             end do
             error  = ts%tq(1) * norm(emcxe1(:,p), ts%ewt(:,p))
             eta(1) = ONE / ( (10.e0_rt * error) ** (ONE / (ts%k+2)) + 1.e-6_rt )
          end if
          ts%k_age = 0
       end if

       ! choose which eta will maximize the time step
       etamax(p) = 0
       if (eta(-1) > etamax(p)) then
          etamax(p) = eta(-1)
          delta(p)  = -1
       end if
       if (eta(1) > etamax(p)) then
          etamax(p) = eta(1)
          delta(p)  = 1
       end if
       if (eta(0) > etamax(p)) then
          etamax(p) = eta(0)
          delta(p)  = 0
       end if
    end do

    p = minloc(etamax)
    rescale = 0
    etaminmax = etamax(p)
    if (etaminmax > ts%eta_thresh) then
       if (delta(p) == -1) then
          call decrease_order(ts)
       else if (delta(p) == 1) then
          call increase_order(ts)
       end if
       rescale = etaminmax
    end if

    if (ts%t + ts%dt > ts%t1) then
       rescale = (ts%t1 - ts%t) / ts%dt
       call rescale_timestep(ts, rescale, .true.)
    else if (rescale /= 0) then
       call rescale_timestep(ts, rescale, .false.)
    end if

    ! save for next step (needed to compute eta(1))
    !ts%e1 = ts%e
    do p = 1, ts%npt
       do m = 1, ts%neq
          ts%e1(m,p) = ts%e(m,p)
       end do
    end do
    ts%tq2save = ts%tq(2)

  end subroutine bdf_adjust

  !
  ! Reset counters, set order to one, init Nordsieck history array.
  !
  subroutine bdf_reset(ts, y0, dt, reuse)
    !$acc routine seq

    use rhs_module, only: rhs

    type(bdf_ts), intent(inout) :: ts
    real(rt),   intent(in   ) :: y0(ts%neq, ts%npt), dt
    logical,      intent(in   ) :: reuse
    
    integer :: p,m,o

    ts%nfe = 0
    ts%nje = 0
    ts%nlu = 0
    ts%nit = 0
    ts%nse = 0

    !ts%y  = y0
    do p = 1, ts%npt
       do m = 1, ts%neq
          ts%y(m,p) = y0(m,p)
       end do
    end do
    ts%dt = dt
    ts%n  = 1
    ts%k  = 1

    !ts%h = ts%dt
    do o = 0, ts%max_order
       ts%h(o) = ts%dt
    enddo
    ts%dt_nwt   = ts%dt
    ts%refactor = .true.

    !print *, 'yd b4: ', ts%yd
    call rhs(ts)
    !print *, 'yd af: ', ts%yd
    ts%nfe = ts%nfe + 1

    !ts%z(:,:,0) = ts%y
    !ts%z(:,:,1) = ts%dt * ts%yd
    do p = 1, ts%npt
       do m = 1, ts%neq
          ts%z(m,p,0) = ts%y(m,p)
          ts%z(m,p,1) = ts%dt * ts%yd(m,p)
       end do
    end do

    ts%k_age = 0
    if (.not. reuse) then
       ts%j_age = ts%max_j_age + 1
       ts%p_age = ts%max_p_age + 1
    else
       ts%j_age = 0
       ts%p_age = 0
    end if

  end subroutine bdf_reset

  !
  ! Rescale time-step.
  !
  ! This consists of:
  !   1. bound eta to honor eta_min, eta_max, and dt_min
  !   2. scale dt and adjust time array t accordingly
  !   3. rescale Nordsieck history array
  !
  subroutine rescale_timestep(ts, eta_in, force)
    !$acc routine seq
    type(bdf_ts), intent(inout) :: ts
    real(rt),   intent(in   ) :: eta_in
    logical,      intent(in   ) :: force

    real(rt) :: eta
    integer  :: i

    if (force) then
       eta = eta_in
    else
       eta = max(eta_in, ts%dt_min / ts%dt, ts%eta_min)
       eta = min(eta, ts%eta_max)

       if (ts%t + eta*ts%dt > ts%t1) then
          eta = (ts%t1 - ts%t) / ts%dt
       end if
    end if

    ts%dt   = eta * ts%dt
    ts%h(0) = ts%dt

    do i = 1, ts%k
       ts%z(:,:,i) = eta**i * ts%z(:,:,i)
    end do
  end subroutine rescale_timestep

  !
  ! Decrease order.
  !
  subroutine decrease_order(ts)
    !$acc routine seq
    type(bdf_ts), intent(inout) :: ts
    integer  :: j, o, p, m
    real(rt) :: c(0:6), c_shift(0:6)

    if (ts%k > 2) then
       do o = 0, 6
          c(o) = 0
       end do
       c(2) = 1
       do j = 1, ts%k-2
          !c = eoshift_local(c, -1) + c * xi_j(ts%h, j)
          !c_shift = eoshift_local(c, -1)
          call eoshift_local(c, -1, c_shift)
          do o = 0, 6
             c(o) = c_shift(o) + c(o) *  xi_j(ts%h, j)
          end do
       end do

       do j = 2, ts%k-1
          !NOTE: We have to explicitly loop here because otherwise this breaks
          !   on the GPUs.  Assignments of form array = array + expr often
          !   require an implicit temporary array that Fortran or compilers generate in the
          !   background.  Such an array must be allocated and is thus not OK
          !   for GPUs.
          do p = 1, ts%npt
             do m = 1, ts%neq
                ts%z(m,p,j) = ts%z(m,p,j) - c(j) * ts%z(m,p,ts%k)
             end do
          end do
       end do
    end if

    ts%z(:,:,ts%k) = 0
    ts%k = ts%k - 1
  end subroutine decrease_order

  !
  ! Increase order.
  !
  subroutine increase_order(ts)
    !$acc routine seq
    type(bdf_ts), intent(inout) :: ts
    integer  :: j, o
    real(rt) :: c(0:6), c_shift(0:6)

    c = 0
    c(2) = 1
    do j = 1, ts%k-2
       !c = eoshift_local(c, -1) + c * xi_j(ts%h, j)
       !c_shift = eoshift_local(c, -1)
       call eoshift_local(c, -1, c_shift)
       do o = 0, 6
          c(o) = c_shift(o) + c(o) * xi_j(ts%h, j)
       end do
    end do

    ts%z(:,:,ts%k+1) = 0
    do j = 2, ts%k+1
       ts%z(:,:,j) = ts%z(:,:,j) + c(j) * ts%e
    end do

    ts%k = ts%k + 1
  end subroutine increase_order

  !
  ! Return $\alpha_0$.
  !
  function alpha0(k) result(a0)
    !$acc routine seq
    integer,  intent(in) :: k
    real(rt) :: a0
    integer  :: j
    a0 = -1
    do j = 2, k
       a0 = a0 - ONE / j
    end do
  end function alpha0

  !
  ! Return $\hat{\alpha}_{n,0}$.
  !
  function alphahat0(k, h) result(a0)
    !$acc routine seq
    integer,  intent(in) :: k
    real(rt), intent(in) :: h(0:k)
    real(rt) :: a0
    integer  :: j
    a0 = -1
    do j = 2, k
       a0 = a0 - h(0) / sum(h(0:j-1))
    end do
  end function alphahat0

  !
  ! Return 1 / $\xi^*_k$.
  !
  ! Note that a lot of simplifications can be made to the formula for
  ! $\xi^*_k$ that appears in Jackson and Sacks-Davis.
  !
  function xi_star_inv(k, h) result(xii)
    !$acc routine seq
    integer,  intent(in) :: k
    real(rt), intent(in) :: h(0:)
    real(rt) :: xii, hs
    integer  :: j
    hs = 0.0_rt
    xii = -alpha0(k)
    do j = 0, k-2
       hs  = hs + h(j)
       xii = xii - h(0) / hs
    end do
  end function xi_star_inv

  !
  ! Return $\xi_j$.
  !
  function xi_j(h, j) result(xi)
    !$acc routine seq
    integer,  intent(in) :: j
    real(rt), intent(in) :: h(0:)
    real(rt) :: xi
    xi = sum(h(0:j-1)) / h(0)
  end function xi_j

  !
  ! Pre-compute error weights.
  !
  subroutine ewts(ts)
    !$acc routine seq
    type(bdf_ts), intent(inout) :: ts
    integer :: m, p
    do p = 1, ts%npt
       do m = 1, ts%neq
          ts%ewt(m,p) = ONE / (ts%rtol(m) * abs(ts%y(m,p)) + ts%atol(m))
       end do
    end do
  end subroutine ewts

  subroutine print_y(ts)
    type(bdf_ts), intent(in) :: ts
    integer :: p
    do p = 1, ts%npt
       print *, ts%y(:,p)
    end do
  end subroutine print_y

  !
  ! Compute weighted norm of y.
  !
  function norm(y, ewt) result(r)
    !$acc routine seq
    real(rt), intent(in) :: y(1:), ewt(1:)
    real(rt) :: r
    integer :: m, n
    n = size(y)
    r = 0.0_rt
    do m = 1, n
       r = r + (y(m)*ewt(m))**2
    end do
    r = sqrt(r/n)
  end function norm

  !
  ! Build/destroy BDF time-stepper.
  ! Note that many operations are explicitly written out
  ! in this subroutine due to limitations with the PGI
  ! compiler intrinsics using the accelerator. Make sure to
  ! try compiling with OpenACC support enabled before
  ! committing any changes to this subroutine.
  !
  subroutine bdf_ts_build(ts)
    !$acc routine seq
    !$acc routine(dgemm) seq

    use extern_probin_module, only : dt_min, jac_age, p_age

    type(bdf_ts),   intent(inout) :: ts

    integer :: i, j, k, n

    ! these are set at build time
    ts%npt = bdf_npt
    ts%neq = neqs
    ts%max_order = bdf_max_order

    ! these are user-controllable
    ts%max_steps  = 1000000
    ts%max_iters  = 10
    ts%verbose    = 0
    ts%dt_min     = dt_min   !epsilon(ts%dt_min)
    ts%eta_min    = 0.2_rt
    ts%eta_max    = 10.0_rt
    ts%eta_thresh = 1.50_rt
    ts%max_j_age  = jac_age
    ts%max_p_age  = p_age

    ts%k = -1

    do n = 1, ts % npt
       do k = 1, neqs
          ts % yd(k,n) = 0.0e0_rt
          do j = 1, neqs
             ts % J(j,k,n) = 0.0e0_rt
             ts % P(j,k,n) = 0.0e0_rt
          enddo
       enddo
    enddo

    ! force a rebuild at the start
    ts%j_age = 666666666
    ts%p_age = 666666666

    ts%debug = .false.

  end subroutine bdf_ts_build

  subroutine bdf_ts_destroy(ts)
    type(bdf_ts), intent(inout) :: ts
  end subroutine bdf_ts_destroy

  !
  ! Various misc. helper functions
  !
  subroutine eye_r(A)
    !$acc routine seq
    real(rt), intent(inout) :: A(:,:,:)
    integer :: i
    A = 0
    do i = 1, size(A, 1)
       A(i,i,:) = 1
    end do
  end subroutine eye_r
  subroutine eye_i(A)
    integer, intent(inout) :: A(:,:)
    integer :: i
    A = 0
    do i = 1, size(A, 1)
       A(i,i) = 1
    end do
  end subroutine eye_i
  recursive function factorial(n) result(r)
    !$acc routine seq
    integer, intent(in) :: n
    integer :: r
    if (n == 1) then
       r = 1
    else
       r = n * factorial(n-1)
    end if
  end function factorial

  !
  ! A local, GPU-compiled version of intrinsic eoshift 
  ! Only what's needed for VBDF is implemented, also no
  ! error-checking.  And we assume 0-based indexing for the arrays as all uses
  ! of eoshift in bdf are with 0-based arrays.
  !
  ! NOTE: Array-valued functions are NOT allowed on the GPU (in PGI at least), had to rewrite this
  ! as a subroutine
  !
  subroutine eoshift_local(arr, sh, shifted_arr)
    !$acc routine seq
    real(rt), intent(in   ) :: arr(0:)
    integer,         intent(in   ) :: sh
    real(rt), intent(  out) :: shifted_arr(0:)
    
    integer :: i, hi_arr, hi_shift

    !TODO: These should be the same size, maybe do a consistency check here
    hi_arr = size(arr) - 1
    hi_shift = size(shifted_arr) - 1

    !shifted_arr = 0.0
    do i = 0, hi_shift
       shifted_arr(i) = 0.0
    enddo

    if(sh > 0) then
       do i = 0, hi_arr - sh
          shifted_arr(i) = arr(i+sh)
       enddo
    else if(sh < 0) then
       do i = hi_arr, abs(sh), -1
          shifted_arr(i) = arr(i+sh)
       enddo
    end if

    ! for debugging
    !shifted_arr(:) = eoshift(arr, sh)

  end subroutine eoshift_local

  !
  ! A local, GPU-compiled version of intrinsic minloc
  ! Only what's needed for VBDF is implemented, also no
  ! error-checking.
  ! TODO: Check if this is implemented on GPU, if so delete all this
  !
  function minloc(arr) result(ret)
    !$acc routine seq
    real(rt), intent(in   ) :: arr(:)
    
    integer :: ret
    integer :: i
    real(rt) :: cur_min

    cur_min = arr(1)
    ret = 1
    do i = 1, size(arr)
      if(arr(i) < cur_min) then
        cur_min = arr(i)
        ret = i
      endif
    enddo
  end function minloc

  !
  ! Initialize the pascal matrix.  We do this here because we should only do
  ! it once and use it for all bdf_ts types.
  !
  subroutine init_pascal()
     ! NOTE: bdf_max_order comes in from bdf_type_module
     integer :: U(bdf_max_order+1, bdf_max_order+1), Uk(bdf_max_order+1, bdf_max_order+1)
     integer :: k, n, r, c, sum_element

     ! build pascal matrix A using A = exp(U)
     !U = 0
     do r = 1, bdf_max_order+1
        do c = 1, bdf_max_order+1
           U(r,c) = 0
        enddo
     enddo
     do k = 1, bdf_max_order
        U(k,k+1) = k
     end do
     !Uk = U
     do r = 1, bdf_max_order+1
        do c = 1, bdf_max_order+1
           Uk(r,c) = U(r,c)
        enddo
     enddo
     ! NOTE: A comes in from bdf_type_module
     call eye_i(A)
     do k = 1, bdf_max_order+1
        A  = A + Uk / factorial(k)
        !TODO: This is an unoptimized, naive matrix multiply, might consider
        !using optimized.  Can't use Fortran intrinsic matmul() on GPU
        !Uk = matmul(U, Uk)
        do r = 1, bdf_max_order+1
           do c = 1, bdf_max_order+1
              sum_element = 0
              do n = 1, bdf_max_order+1
                 sum_element = sum_element + U(r,n) * Uk(n,c)
              enddo
              Uk(r,c) = sum_element
           enddo
        enddo
     end do
     !$acc update device(A)
  end subroutine init_pascal

end module bdf
