module sdc_quadrature_module

  use sdc_sizes_module
  use amrex_fort_module, only : rt => amrex_real
  use sdc_parameters_module
  use amrex_constants_module

  implicit none

contains

  subroutine sdc_C_source(m_start, dt, dt_m, f_old, C)

    use amrex_error_module, only : amrex_error

    implicit none

    integer, intent(in) :: m_start
    real(rt), intent(in) :: dt, dt_m
    real(rt), intent(in) :: f_old(0:SDC_NODES-1, SDC_NEQS)
    real(rt), intent(inout) :: C(SDC_NEQS)

    real(rt) :: integral(SDC_NEQS)

    if (m_start == 0) then
       ! compute the integral from [t_m, t_{m+1}], normalized by dt_m
       integral(:) = (dt/dt_m) * (ONE/1800.0_rt) * &
            ((-35.0_rt*sqrt(6.0_rt) + 440.0_rt)*f_old(1, :) + &
            (-169.0_rt*sqrt(6.0_rt) + 296.0_rt)*f_old(2, :) + &
            (-16.0_rt + 24.0_rt*sqrt(6.0_rt))*f_old(3, :))

       C(:) = - f_old(1, :) + integral

    else if (m_start == 1) then
       ! compute the integral from [t_m, t_{m+1}], normalized by dt_m
       integral(:) = (dt/dt_m) * (ONE/150.0_rt) * &
            ((-12.0_rt + 17.0_rt*sqrt(6.0_rt))*f_old(1, :) + &
            (12.0_rt + 17.0_rt*sqrt(6.0_rt))*f_old(2, :) + &
            (-4.0_rt*sqrt(6.0_rt))*f_old(3, :))

       C(:) = - f_old(2, :) + integral

    else if (m_start == 2) then
       ! compute the integral from [t_m, t_{m+1}], normalized by dt_m
       integral(:) = (dt/dt_m) * (ONE/600.0_rt) * &
            ((168.0_rt - 73.0_rt*sqrt(6.0_rt))*f_old(1, :) + &
            (120.0_rt + 5.0_rt*sqrt(6.0_rt))*f_old(2, :) + &
            (72.0_rt + 8.0_rt*sqrt(6.0_rt))*f_old(3, :))

       C(:) = - f_old(3, :) + integral

    else
       call amrex_error("error: invalid node in sdc_C_source")
    endif

  end subroutine sdc_C_source

  subroutine sdc_newton_solve(dt_m, sdc, y_new, f_source, sdc_iteration, ierr)
    ! the purpose of this function is to solve the system
    ! U - dt R(U) = U_old + dt C using a Newton solve.
    !
    ! here, U_new should come in as a guess for the new U and will be
    ! returned with the value that satisfies the nonlinear function

    use amrex_constants_module, only : ZERO, HALF, ONE
    use extern_probin_module, only : small_x
    use rhs_module, only : f_rhs, jac
    use sdc_parameters_module, only : SDC_NEQS

    implicit none

    real(rt), intent(in) :: dt_m
    type(sdc_t) :: sdc
    real(rt), intent(inout) :: y_new(SDC_NEQS)
    real(rt), intent(in) :: f_source(SDC_NEQS)
    integer, intent(in) :: sdc_iteration
    integer, intent(out) :: ierr

    real(rt) :: A(SDC_NEQS, SDC_NEQS)
    real(rt) :: rhs(SDC_NEQS)
    real(rt) :: dy(SDC_NEQS)

    real(rt) :: w(SDC_NEQS)

    integer :: ipvt(SDC_NEQS)
    integer :: info

    logical :: converged

    real(rt) :: eps_tot(SDC_NEQS)

    real(rt) :: err, eta

    integer, parameter :: MAX_ITER = 100
    integer :: iter

    integer :: max_newton_iter

    integer :: k, n

    ierr = NEWTON_SUCCESS

    ! do a simple Newton solve

    ! iterative loop
    iter = 0
    max_newton_iter = MAX_ITER

    err = 1.e30_rt
    converged = .false.

    do while (.not. converged .and. iter < max_newton_iter)

       ! get the Jacobian and the RHS
       call f_rhs(sdc)
       call jac(sdc)

       ! create the matrix for our linear system
       A(:,:) = 0.0_rt
       do n = 1, SDC_NEQS
          A(n, n) = 1.0_rt
       end do

       A(:,:) = A(:,:) - dt_m * sdc % burn_s % jac(:, :)

       rhs(:) = -sdc % y(:) + dt_m * sdc % burn_s % ydot(:) - f_source(:)

       ! solve the linear system: A dy_react = rhs
       call dgefa(A, SDC_NEQS, SDC_NEQS, ipvt, info)
       if (info /= 0) then
          ierr = SINGULAR_MATRIX
          return
       endif

       call dgesl(A, SDC_NEQS, SDC_NEQS, ipvt, rhs, 0)

       dy(:) = rhs(:)

       ! how much of dU_react should we apply?
       eta = ONE

       sdc % y(:) = sdc % y(:) + dy(:)

       ! we still need to normalize here
       call clean_state(sdc)

       eps_tot(:) = sdc % rtol(:) * abs(sdc % y(:)) + sdc % atol(:)

       ! compute the norm of the weighted error, where the weights are 1/eps_tot
       err = sqrt(sum((dy/eps_tot)**2)/SDC_NEQS)

       if (err < ONE) then
          converged = .true.
       endif

       iter = iter + 1
    enddo

    if (.not. converged) then
       ierr = CONVERGENCE_FAILURE
       return
    endif

  end subroutine sdc_newton_solve

end module sdc_quadrature_module
