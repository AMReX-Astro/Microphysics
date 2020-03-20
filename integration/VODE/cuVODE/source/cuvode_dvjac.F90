module cuvode_dvjac_module

  use cuvode_parameters_module, only: VODE_NEQS
  use cuvode_types_module, only: dvode_t, UROUND, SRUR, max_steps_between_jacobian_evals, &
                                 CCMXJ
  use amrex_fort_module, only: rt => amrex_real
  use linpack_module

  implicit none

contains

#if defined(AMREX_USE_CUDA) && !defined(AMREX_USE_GPU_PRAGMA)
  attributes(device) &
#endif
  subroutine dvjac(pivot, IERPJ, vstate)

    !$acc routine seq
    
    ! -----------------------------------------------------------------------
    !  DVJAC is called by DVNLSD to compute and process the matrix
    !  P = I - h*rl1*J , where J is an approximation to the Jacobian
    !  that we obtain either through direct evaluation or caching from
    !  a previous evaluation. P is then subjected to LU decomposition
    !  in preparation for later solution of linear systems with P as
    !  coefficient matrix. This is done by DGEFA.
    ! -----------------------------------------------------------------------

#ifdef TRUE_SDC
    use sdc_vode_rhs_module, only: f_rhs, jac
#else
    use vode_rhs_module, only: f_rhs, jac
#endif

    implicit none

    ! Declare arguments
    type(dvode_t), intent(inout) :: vstate
    integer,       intent(inout) :: pivot(VODE_NEQS)
    integer,       intent(  out) :: IERPJ

    ! Declare local variables
    real(rt) :: con, fac, hrl1, R, R0, yj
    integer  :: i, j, j1, IER, II, evaluate_jacobian

    !$gpu

    IERPJ = 0

    ! See whether the Jacobian should be evaluated. Start by basing
    ! the decision on whether we're caching the Jacobian.
    evaluate_jacobian = -vstate % JSV

    if (vstate % JSV == 1) then
       ! Now evaluate the cases where we're caching the Jacobian but aren't
       ! going to be using the cached Jacobian.

       ! On the first step we don't have a cached Jacobian. Also, after enough
       ! steps, we consider the cached Jacobian too old and will want to re-evaluate
       ! it, so we look at whether the step of the last Jacobian evaluation (NSLJ)
       ! is more than max_steps_between_jacobian_evals steps in the past.
       if (vstate % NST == 0 .or. vstate % NST > vstate % NSLJ + max_steps_between_jacobian_evals) then
          evaluate_jacobian = 1
       end if

       ! See the non-linear solver for details on these conditions.
       if (vstate % ICF == 1 .and. vstate % DRC .LT. CCMXJ) then
          evaluate_jacobian = 1
       end if

       if (vstate % ICF == 2) then
          evaluate_jacobian = 1
       end if
    end if

    if (evaluate_jacobian == 1) then

       ! We want to evaluate the Jacobian -- now the path depends on
       ! whether we're using the numerical or analytic Jacobian.

       if (vstate % jacobian == 1) then

          ! For the analytic Jacobian, call the user-supplied function.

          ! Increment the Jacobian evaluation counter.
          vstate % NJE = vstate % NJE + 1

          ! Refresh the timestep marker for the last Jacobian evaluation.
          vstate % NSLJ = vstate % NST

          ! Indicate that the Jacobian is current for this solve.
          vstate % JCUR = 1

          do i = 1, VODE_NEQS * VODE_NEQS
             vstate % JAC(i) = 0.0_rt
          end do

          call jac(vstate % tn, vstate, 0, 0, vstate % jac, VODE_NEQS)

          ! Store the Jacobian if we're caching.
          if (vstate % JSV == 1) then
             do i = 1, VODE_NEQS * VODE_NEQS
                vstate % jac_save(i) = vstate % jac(i)
             end do
          end if

       else

          ! For the numerical Jacobian, make N calls to the RHS to approximate it.

          ! Increment the Jacobian evaluation counter.
          vstate % NJE = vstate % NJE + 1

          ! Refresh the timestep marker for the last Jacobian evaluation.
          vstate % NSLJ = vstate % NST

          ! Indicate that the Jacobian is current for this solve.
          vstate % JCUR = 1

          fac = sqrt(sum((vstate % savf * vstate % ewt)**2) / VODE_NEQS)
          R0 = 1000.0_rt * abs(vstate % H) * UROUND * real(VODE_NEQS) * fac
          if (R0 == 0.0_rt) then
             R0 = 1.0_rt
          end if

          j1 = 0
          do j = 1, VODE_NEQS
             yj = vstate % y(j)

             R = max(SRUR * abs(yj), R0 / vstate % EWT(j))
             vstate % y(j) = vstate % y(j) + R
             fac = 1.0_rt/R

             call f_rhs(vstate % tn, vstate, vstate % acor)
             do i = 1, VODE_NEQS
                vstate % jac(i+j1) = (vstate % acor(i) - vstate % SAVF(i)) * fac
             end do

             vstate % y(j) = yj

             j1 = j1 + VODE_NEQS
          end do

          ! Increment the RHS evaluation counter by N.
          vstate % NFE = vstate % NFE + VODE_NEQS

          ! Store the Jacobian if we're caching.
          if (vstate % JSV == 1) then
             do i = 1, VODE_NEQS * VODE_NEQS
                vstate % jac_save(i) = vstate % jac(i)
             end do
          end if

       end if

    else

       ! Load the cached Jacobian.

       ! Indicate the Jacobian is not current for this step.
       vstate % JCUR = 0
       do i = 1, VODE_NEQS * VODE_NEQS
          vstate % jac(i) = vstate % jac_save(i)
       end do

    end if

    ! Multiply Jacobian by a scalar, add the identity matrix
    ! (along the diagonal), and do LU decomposition.

    hrl1 = vstate % H * vstate % RL1
    con = -hrl1
    vstate % jac(:) = vstate % jac(:) * con

    j = 1
    do i = 1, VODE_NEQS
       ! Add 1 from the identity matrix; then,
       ! advance to the next location on the diagonal,
       ! which is in the next row (advance by VODE_NEQS)
       ! and over by one column (advance by 1).
       vstate % jac(j) = vstate % jac(j) + 1.0_rt
       j = j + VODE_NEQS + 1
    end do

    call dgefa(vstate % jac, pivot, IER)

    if (IER /= 0) IERPJ = 1

  end subroutine dvjac

end module cuvode_dvjac_module
