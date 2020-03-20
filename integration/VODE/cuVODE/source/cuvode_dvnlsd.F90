module cuvode_dvnlsd_module

  use cuvode_parameters_module, only: VODE_LMAX, VODE_NEQS, VODE_MAXORD
  use cuvode_types_module, only: dvode_t
  use amrex_fort_module, only: rt => amrex_real
  use linpack_module
  use cuvode_dvjac_module

  implicit none

contains

#if defined(AMREX_USE_CUDA) && !defined(AMREX_USE_GPU_PRAGMA)
  attributes(device) &
#endif
  subroutine dvnlsd(pivot, NFLAG, vstate)

    !$acc routine seq

    ! -----------------------------------------------------------------------
    !  Subroutine DVNLSD is a nonlinear system solver, which uses functional
    !  iteration or a chord (modified Newton) method.  For the chord method
    !  direct linear algebraic system solvers are used.  Subroutine DVNLSD
    !  then handles the corrector phase of this integration package.
    ! -----------------------------------------------------------------------

#ifdef TRUE_SDC
    use sdc_vode_rhs_module, only: f_rhs, jac
#else
    use vode_rhs_module, only: f_rhs, jac
#endif
#ifdef CLEAN_INTEGRATOR_CORRECTION
    use vode_type_module, only: clean_state
#endif

    implicit none

    ! Declare arguments
    type(dvode_t), intent(inout) :: vstate
    integer,       intent(inout) :: NFLAG
    integer,       intent(inout) :: pivot(VODE_NEQS)

    ! Declare local variables
    real(rt) :: CSCALE, DCON, DEL, DELP
    integer  :: I, IERPJ, M

    ! Parameter declarations
    real(rt), parameter :: CCMAX = 0.3e0_rt
    real(rt), parameter :: CRDOWN = 0.3e0_rt
    real(rt), parameter :: RDIV  = 2.0e0_rt
    integer, parameter :: MAXCOR = 3
    integer, parameter :: MSBP = 20

    logical :: converged

    !$gpu

    ! -----------------------------------------------------------------------
    !  On the first step, on a change of method order, or after a
    !  nonlinear convergence failure with NFLAG = -2, set IPUP = 1
    !  to force a Jacobian update.
    ! -----------------------------------------------------------------------
    if (vstate % JSTART == 0) then
       vstate % NSLP = 0
    end if
    if (NFLAG == 0) then
       vstate % ICF = 0
    end if
    if (NFLAG == -2) then
       vstate % IPUP = 1
    end if
    if (vstate % JSTART == 0) then
       vstate % IPUP = 1
    end if

    ! -----------------------------------------------------------------------
    !  RC is the ratio of new to old values of the coefficient H/EL(2)=h/l1.
    !  When RC differs from 1 by more than CCMAX, IPUP is set to 1
    !  to force DVJAC to be called, if a Jacobian is involved.
    !  In any case, DVJAC is called at least every MSBP steps.
    ! -----------------------------------------------------------------------
    vstate % DRC = abs(vstate % RC - 1.0_rt)
    if (vstate % DRC > CCMAX .or. vstate % NST >= vstate % NSLP+MSBP) then
       vstate % IPUP = 1
    end if

    ! -----------------------------------------------------------------------
    !  Up to MAXCOR corrector iterations are taken.  A convergence test is
    !  made on the r.m.s. norm of each correction, weighted by the error
    !  weight array EWT.  The sum of the corrections is accumulated in the
    !  array ACOR(i).  The YH array is not altered in the corrector loop.
    ! -----------------------------------------------------------------------

    converged = .false.

    do while (.true.)

       M = 0
       DELP = 0.0_rt

       vstate % Y(1:VODE_NEQS) = vstate % yh(1:VODE_NEQS,1)
       call f_rhs(vstate % TN, vstate, vstate % savf)
       vstate % NFE = vstate % NFE + 1

       if (vstate % IPUP == 1) then
          ! -----------------------------------------------------------------------
          !  If indicated, the matrix P = I - h*rl1*J is reevaluated and
          !  preprocessed before starting the corrector iteration.  IPUP is set
          !  to 0 as an indicator that this has been done.
          ! -----------------------------------------------------------------------
          call dvjac(pivot, IERPJ, vstate)
          vstate % IPUP = 0
          vstate % RC = 1.0_rt
          vstate % DRC = 0.0_rt
          vstate % CRATE = 1.0_rt
          vstate % NSLP = vstate % NST

          ! If matrix is singular, take error return to force cut in step size. --
          if (IERPJ /= 0) then
             NFLAG = -1
             vstate % ICF = 2
             vstate % IPUP = 1
             return
          end if

       end if

       do I = 1, VODE_NEQS
          vstate % acor(I) = 0.0_rt
       end do

       ! This is a looping point for the corrector iteration.

       do while (.true.)

          ! -----------------------------------------------------------------------
          !  In the case of the chord method, compute the corrector error,
          !  and solve the linear system with that as right-hand side and
          !  P as coefficient matrix.  The correction is scaled by the factor
          !  2/(1+RC) to account for changes in h*rl1 since the last DVJAC call.
          ! -----------------------------------------------------------------------

          do I = 1, VODE_NEQS
             vstate % Y(I) = (vstate % RL1 * vstate % H) * vstate % SAVF(I) - &
                             (vstate % RL1 * vstate % YH(I,2) + vstate % ACOR(I))
          end do

          call dgesl(vstate % jac, pivot, vstate % Y(:))

          vstate % NNI = vstate % NNI + 1

          if (vstate % RC /= 1.0_rt) then
             CSCALE = 2.0_rt / (1.0_rt + vstate % RC)
             vstate % Y(:) = vstate % Y(:) * CSCALE
          end if

#ifdef CLEAN_INTEGRATOR_CORRECTION
          ! Clean the correction to Y. Use vstate % Y as scratch space.

          ! Find the corrected Y: Yc = Y_previous + Y_delta
          do I = 1, VODE_NEQS
             vstate % Y(I) = vstate % Y(I) + (vstate % YH(I,1) + vstate % ACOR(I))
          end do

          ! Clean the corrected Y: Yc' = clean(Yc)
          call clean_state(vstate % Y, vstate % RPAR)

          ! Find the cleaned correction: clean(Y_delta) = Yc' - Y_previous
          do I = 1,VODE_NEQS
             vstate % Y(I) = vstate % Y(I) - (vstate % YH(I,1) + vstate % ACOR(I))
          end do
#endif

          DEL = sqrt(sum((vstate % Y * vstate % EWT)**2) / VODE_NEQS)
          vstate % acor(:) = vstate % acor(:) + vstate % Y(:)

          do I = 1, VODE_NEQS
             vstate % Y(I) = vstate % YH(I,1) + vstate % ACOR(I)
          end do

          ! -----------------------------------------------------------------------
          !  Test for convergence.  If M .gt. 0, an estimate of the convergence
          !  rate constant is stored in CRATE, and this is used in the test.
          ! -----------------------------------------------------------------------

          if (M /= 0) vstate % CRATE = max(CRDOWN * vstate % CRATE, DEL / DELP)
          DCON = DEL * min(1.0_rt, vstate % CRATE) / vstate % TQ(4)
          if (DCON <= 1.0_rt) then
             ! we converged, exit the outer loop
             converged = .true.
             exit
          end if

          M = M + 1
          if (M == MAXCOR) then
             ! exit the inner correction iteration
             exit
          end if
          if (M >= 2 .and. DEL > RDIV*DELP) then
             ! exit the inner correction iteration
             exit
          end if
          DELP = DEL
          call f_rhs(vstate % TN, vstate, vstate % SAVF)
          vstate % NFE = vstate % NFE + 1
       end do

       if (converged) then
          exit
       end if

       if (vstate % JCUR == 1) then
          NFLAG = -1
          vstate % ICF = 2
          vstate % IPUP = 1
          return
       end if

       vstate % ICF = 1
       vstate % IPUP = 1

    end do

    ! Return for successful step.
    NFLAG = 0
    vstate % JCUR = 0
    vstate % ICF = 0
    if (M == 0) then
       vstate % ACNRM = DEL
    end if
    if (M > 0) then
       vstate % ACNRM = sqrt(sum((vstate % ACOR * vstate % EWT)**2) / VODE_NEQS)
    end if

  end subroutine dvnlsd

end module cuvode_dvnlsd_module
