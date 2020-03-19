module cuvode_dvstep_module

  use cuvode_parameters_module, only: VODE_LMAX, VODE_NEQS, VODE_LIW,   &
                                      VODE_LENWM, VODE_MAXORD
  use cuvode_types_module, only: dvode_t, rwork_t
  use amrex_fort_module, only: rt => amrex_real

  use cuvode_dvset_module
  use cuvode_dvjust_module
  use cuvode_dvnlsd_module
  use cuvode_nordsieck_module

  use cuvode_constants_module

  implicit none

contains

#if defined(AMREX_USE_CUDA) && !defined(AMREX_USE_GPU_PRAGMA)
  attributes(device) &
#endif
  subroutine dvstep(pivot, rwork, vstate)

    !$acc routine seq

    ! -----------------------------------------------------------------------
    !  DVSTEP performs one step of the integration of an initial value
    !  problem for a system of ordinary differential equations.
    !  DVSTEP calls subroutine VNLS for the solution of the nonlinear system
    !  arising in the time step.  Thus it is independent of the problem
    !  Jacobian structure and the type of nonlinear system solution method.
    !  DVSTEP returns a completion flag KFLAG (in COMMON).
    !  A return with KFLAG = -1 or -2 means either ABS(H) = HMIN or 10
    !  consecutive failures occurred.  On a return with KFLAG negative,
    !  the values of TN and the YH array are as of the beginning of the last
    !  step, and H is the last step size attempted.
    !
    !  Communication with DVSTEP is done with the following variables:
    !
    !  Y      = An array of length N used for the dependent variable array.
    !  YH     = An LDYH by LMAX array containing the dependent variables
    !           and their approximate scaled derivatives, where
    !           LMAX = MAXORD + 1.  YH(i,j+1) contains the approximate
    !           j-th derivative of y(i), scaled by H**j/factorial(j)
    !           (j = 0,1,...,NQ).  On entry for the first step, the first
    !           two columns of YH must be set from the initial values.
    !  LDYH   = A constant integer .ge. N, the first dimension of YH.
    !           N is the number of ODEs in the system.
    !  EWT    = An array of length N containing multiplicative weights
    !           for local error measurements.  Local errors in y(i) are
    !           compared to 1.0/EWT(i) in various error tests.
    !  SAVF   = An array of working storage, of length N.
    !           also used for input of YH(*,MAXORD+2) when JSTART = -1
    !           and MAXORD .lt. the current order NQ.
    !  VSAV   = A work array of length N passed to subroutine VNLS.
    !  ACOR   = A work array of length N, used for the accumulated
    !           corrections.  On a successful return, ACOR(i) contains
    !           the estimated one-step local error in y(i).
    !  F      = Dummy name for the user supplied subroutine for f.
    !  JAC    = Dummy name for the user supplied Jacobian subroutine.
    !  PSOL   = Dummy name for the subroutine passed to VNLS, for
    !           possible use there.
    !  VNLS   = Dummy name for the nonlinear system solving subroutine,
    !           whose real name is dependent on the method used.
    !  RPAR, IPAR = Dummy names for user's real and integer work arrays.
    ! -----------------------------------------------------------------------

    ! -----------------------------------------------------------------------
    !  On a successful return, ETAMAX is reset and ACOR is scaled.
    ! -----------------------------------------------------------------------

#ifdef TRUE_SDC
    use sdc_vode_rhs_module, only: f_rhs, jac
#else
    use vode_rhs_module, only: f_rhs, jac
#endif
    use cuvode_dvnorm_module, only: dvnorm ! function

    implicit none

    ! Declare arguments
    type(dvode_t), intent(inout) :: vstate
    type(rwork_t), intent(inout) :: rwork
    integer,       intent(inout) :: pivot(VODE_NEQS)

    ! Declare local variables
    real(rt) :: CNQUOT, DDN, DSM, DUP, TOLD
    real(rt) :: ETAQ, ETAQM1, ETAQP1, FLOTL, R
    integer    :: I, I1, I2, IBACK, J, JB, NCF, NFLAG

    ! Parameter declarations
    integer, parameter :: KFC = -3
    integer, parameter :: KFH = -7
    integer, parameter :: MXNCF = 10
    real(rt), parameter :: ADDON = 1.0e-6_rt
    real(rt), parameter :: BIAS1 = 6.0e0_rt
    real(rt), parameter :: BIAS2 = 6.0e0_rt
    real(rt), parameter :: BIAS3 = 10.0e0_rt
    real(rt), parameter :: ETACF = 0.25e0_rt
    real(rt), parameter :: ETAMIN = 0.1e0_rt
    real(rt), parameter :: ETAMXF = 0.2e0_rt
    real(rt), parameter :: ETAMX1 = 1.0e4_rt
    real(rt), parameter :: ETAMX2 = 10.0e0_rt
    real(rt), parameter :: ETAMX3 = 10.0e0_rt
    real(rt), parameter :: ONEPSM = 1.00001e0_rt
    real(rt), parameter :: THRESH = 1.5e0_rt

    logical :: do_initialization
    logical :: already_set_eta
    !$gpu

    ETAQ   = ONE
    ETAQM1 = ONE

    vstate % KFLAG = 0
    TOLD = vstate % TN
    NCF = 0
    vstate % JCUR = 0
    NFLAG = 0

    do_initialization = .false.

    if (vstate % jstart >= 0) then

       IF (vstate % JSTART == 0) then

          ! -----------------------------------------------------------------------
          !  On the first call, the order is set to 1, and other variables are
          !  initialized.  ETAMAX is the maximum ratio by which H can be increased
          !  in a single step.  It is normally 10, but is larger during the
          !  first step to compensate for the small initial H.  If a failure
          !  occurs (in corrector convergence or error test), ETAMAX is set to 1
          !  for the next increase.
          ! -----------------------------------------------------------------------
          vstate % NQ = 1
          vstate % L = 2
          vstate % NQNYH = vstate % NQ * VODE_NEQS
          vstate % TAU(1) = vstate % H
          vstate % PRL1 = ONE
          vstate % RC = ZERO
          vstate % ETAMAX = ETAMX1
          vstate % NQWAIT = 2
          vstate % HSCAL = vstate % H

       else
          ! -----------------------------------------------------------------------
          !  Take preliminary actions on a normal continuation step (JSTART.GT.0).
          !  If the driver changed H, then ETA must be reset and NEWH set to 1.
          !  If a change of order was dictated on the previous step, then
          !  it is done here and appropriate adjustments in the history are made.
          !  On an order decrease, the history array is adjusted by DVJUST.
          !  On an order increase, the history array is augmented by a column.
          !  On a change of step size H, the history array YH is rescaled.
          ! -----------------------------------------------------------------------

          IF (vstate % KUTH .EQ. 1) THEN
             vstate % ETA = MIN(vstate % ETA,vstate % H/vstate % HSCAL)
             vstate % NEWH = 1
          ENDIF
          IF (vstate % NEWH .EQ. 0) then
             do_initialization = .false.
          else

             IF (vstate % NEWQ .LT. vstate % NQ) THEN
                CALL DVJUST (-1, rwork, vstate)
                vstate % NQ = vstate % NEWQ
                vstate % L = vstate % NQ + 1
                vstate % NQWAIT = vstate % L
             else IF (vstate % NEWQ .GT. vstate % NQ) THEN
                CALL DVJUST (1, rwork, vstate)
                vstate % NQ = vstate % NEWQ
                vstate % L = vstate % NQ + 1
                vstate % NQWAIT = vstate % L
             ENDIF

             do_initialization = .true.
          end IF

       end IF

    else

       ! -----------------------------------------------------------------------
       !  The following block handles preliminaries needed when JSTART = -1.
       !  If N was reduced, zero out part of YH to avoid undefined references.
       !  If MAXORD was reduced to a value less than the tentative order NEWQ,
       !  then NQ is set to MAXORD, and a new H ratio ETA is chosen.
       !  Otherwise, we take the same preliminary actions as for JSTART .gt. 0.
       !  In any case, NQWAIT is reset to L = NQ + 1 to prevent further
       !  changes in order for that many steps.
       !  The new H ratio ETA is limited by the input H if KUTH = 1,
       !  by HMIN if KUTH = 0, and by HMXI in any case.
       !  Finally, the history array YH is rescaled.
       ! -----------------------------------------------------------------------

       IF (vstate % NEWQ > VODE_MAXORD) then
          FLOTL = REAL(VODE_LMAX)
          IF (VODE_MAXORD .LT. vstate % NQ-1) THEN
             DDN = DVNORM (rwork % SAVF, rwork % EWT)/vstate % TQ(1)
             vstate % ETA = ONE/((BIAS1*DDN)**(ONE/FLOTL) + ADDON)
          ENDIF
          IF (VODE_MAXORD .EQ. vstate % NQ .AND. vstate % NEWQ .EQ. vstate % NQ+1) vstate % ETA = ETAQ
          IF (VODE_MAXORD .EQ. vstate % NQ-1 .AND. vstate % NEWQ .EQ. vstate % NQ+1) THEN
             vstate % ETA = ETAQM1
             CALL DVJUST (-1, rwork, vstate)
          ENDIF
          IF (VODE_MAXORD .EQ. vstate % NQ-1 .AND. vstate % NEWQ .EQ. vstate % NQ) THEN
             DDN = DVNORM (rwork % SAVF, rwork % EWT)/vstate % TQ(1)
             vstate % ETA = ONE/((BIAS1*DDN)**(ONE/FLOTL) + ADDON)
             CALL DVJUST (-1, rwork, vstate)
          ENDIF
          vstate % ETA = MIN(vstate % ETA,ONE)
          vstate % NQ = VODE_MAXORD
          vstate % L = VODE_LMAX
       end IF
       IF (vstate % KUTH .EQ. 1) vstate % ETA = MIN(vstate % ETA,ABS(vstate % H/vstate % HSCAL))
       IF (vstate % KUTH .EQ. 0) vstate % ETA = MAX(vstate % ETA,vstate % HMIN/ABS(vstate % HSCAL))
       vstate % ETA = vstate % ETA/MAX(ONE,ABS(vstate % HSCAL)*vstate % HMXI*vstate % ETA)
       vstate % NEWH = 1
       vstate % NQWAIT = vstate % L
       IF (vstate % NEWQ .LE. VODE_MAXORD) then
          IF (vstate % NEWH .EQ. 0) then
             do_initialization = .false.
          else
             IF (vstate % NEWQ .LT. vstate % NQ) THEN
                CALL DVJUST (-1, rwork, vstate)
                vstate % NQ = vstate % NEWQ
                vstate % L = vstate % NQ + 1
                vstate % NQWAIT = vstate % L
             else if (vstate % NEWQ .GT. vstate % NQ) THEN
                CALL DVJUST (1, rwork, vstate)
                vstate % NQ = vstate % NEWQ
                vstate % L = vstate % NQ + 1
                vstate % NQWAIT = vstate % L
             ENDIF

             do_initialization = .true.
          end IF
       end IF
    end if

    if (do_initialization) then
       ! Rescale the history array for a change in H by a factor of ETA. ------

       R = ONE

       do J = 2, vstate % L
          R = R * vstate % ETA
          rwork % YH(:,J) = rwork % YH(:,J) * R
       end do
       vstate % H = vstate % HSCAL * vstate % ETA
       vstate % HSCAL = vstate % H
       vstate % RC = vstate % RC * vstate % ETA
       vstate % NQNYH = vstate % NQ*VODE_NEQS
    end if

    ! -----------------------------------------------------------------------
    !  This section computes the predicted values by effectively
    !  multiplying the YH array by the Pascal triangle matrix.
    !  DVSET is called to calculate all integration coefficients.
    !  RC is the ratio of new to old values of the coefficient H/EL(2)=h/l1.
    ! -----------------------------------------------------------------------

    do while (.true.)

       vstate % TN = vstate % TN + vstate % H

       call advance_nordsieck(rwork, vstate)

       CALL DVSET(vstate)
       vstate % RL1 = ONE/vstate % EL(2)
       vstate % RC = vstate % RC * (vstate % RL1/vstate % PRL1)
       vstate % PRL1 = vstate % RL1
       !
       !  Call the nonlinear system solver. ------------------------------------
       !
       call dvnlsd(pivot, NFLAG, rwork, vstate)

       IF (NFLAG /= 0) then
          ! -----------------------------------------------------------------------
          !  The VNLS routine failed to achieve convergence (NFLAG .NE. 0).
          !  The YH array is retracted to its values before prediction.
          !  The step size H is reduced and the step is retried, if possible.
          !  Otherwise, an error exit is taken.
          ! -----------------------------------------------------------------------
          NCF = NCF + 1
          vstate % NCFN = vstate % NCFN + 1
          vstate % ETAMAX = ONE
          vstate % TN = TOLD

          call retract_nordsieck(rwork, vstate)

          IF (NFLAG .LT. -1) then
             IF (NFLAG .EQ. -2) vstate % KFLAG = -3
             IF (NFLAG .EQ. -3) vstate % KFLAG = -4
             vstate % JSTART = 1
             RETURN
          end IF
          IF (ABS(vstate % H) .LE. vstate % HMIN*ONEPSM) then
             vstate % KFLAG = -2
             vstate % JSTART = 1
             RETURN
          end IF
          IF (NCF .EQ. MXNCF) then
             vstate % KFLAG = -2
             vstate % JSTART = 1
             RETURN
          end IF
          vstate % ETA = ETACF
          vstate % ETA = MAX(vstate % ETA,vstate % HMIN/ABS(vstate % H))
          NFLAG = -1

          ! Rescale the history array for a change in H by a factor of ETA. ------
          R = ONE

          do J = 2, vstate % L
             R = R * vstate % ETA
             rwork % YH(:,J) = rwork % YH(:,J) * R
          end do
          vstate % H = vstate % HSCAL * vstate % ETA
          vstate % HSCAL = vstate % H
          vstate % RC = vstate % RC * vstate % ETA
          vstate % NQNYH = vstate % NQ*VODE_NEQS

          cycle
       end IF

       ! -----------------------------------------------------------------------
       !  The corrector has converged (NFLAG = 0).  The local error test is
       !  made and control passes to statement 500 if it fails.
       ! -----------------------------------------------------------------------

       DSM = vstate % ACNRM/vstate % TQ(2)
       IF (DSM <= ONE) then
          ! -----------------------------------------------------------------------
          !  After a successful step, update the YH and TAU arrays and decrement
          !  NQWAIT.  If NQWAIT is then 1 and NQ .lt. MAXORD, then ACOR is saved
          !  for use in a possible order increase on the next step.
          !  If ETAMAX = 1 (a failure occurred this step), keep NQWAIT .ge. 2.
          ! -----------------------------------------------------------------------
          vstate % KFLAG = 0
          vstate % NST = vstate % NST + 1
          vstate % HU = vstate % H
          vstate % NQU = vstate % NQ
          do IBACK = 1, vstate % NQ
             I = vstate % L - IBACK
             vstate % TAU(I+1) = vstate % TAU(I)
          end do
          vstate % TAU(1) = vstate % H
          do J = 1, vstate % L
             rwork % yh(:,J) = rwork % yh(:,J) + vstate % EL(J) * rwork % acor(:)
          end do
          vstate % NQWAIT = vstate % NQWAIT - 1
          IF ((vstate % L /= VODE_LMAX) .and. (vstate % NQWAIT == 1)) then
             rwork % yh(1:VODE_NEQS,VODE_LMAX) = rwork % acor(1:VODE_NEQS)

             vstate % CONP = vstate % TQ(5)
          end IF
          IF (vstate % ETAMAX .NE. ONE) exit
          IF (vstate % NQWAIT .LT. 2) vstate % NQWAIT = 2
          vstate % NEWQ = vstate % NQ
          vstate % NEWH = 0
          vstate % ETA = ONE
          vstate % HNEW = vstate % H
          vstate % ETAMAX = ETAMX3
          IF (vstate % NST .LE. 10) vstate % ETAMAX = ETAMX2
          R = ONE/vstate % TQ(2)
          rwork % acor(:) = rwork % acor(:) * R
          vstate % JSTART = 1
          RETURN

       endif

       ! -----------------------------------------------------------------------
       !  The error test failed.  KFLAG keeps track of multiple failures.
       !  Restore TN and the YH array to their previous values, and prepare
       !  to try the step again.  Compute the optimum step size for the
       !  same order.  After repeated failures, H is forced to decrease
       !  more rapidly.
       ! -----------------------------------------------------------------------

       vstate % KFLAG = vstate % KFLAG - 1
       vstate % NETF = vstate % NETF + 1
       NFLAG = -2
       vstate % TN = TOLD

       call retract_nordsieck(rwork, vstate)

       IF (ABS(vstate % H) .LE. vstate % HMIN*ONEPSM) then
          vstate % KFLAG = -1
          vstate % JSTART = 1
          RETURN
       end IF
       vstate % ETAMAX = ONE
       IF (vstate % KFLAG > KFC) then
          ! Compute ratio of new H to current H at the current order. ------------
          FLOTL = REAL(vstate % L)
          vstate % ETA = ONE/((BIAS2*DSM)**(ONE/FLOTL) + ADDON)
          vstate % ETA = MAX(vstate % ETA,vstate % HMIN/ABS(vstate % H),ETAMIN)
          IF ((vstate % KFLAG .LE. -2) .AND. (vstate % ETA .GT. ETAMXF)) vstate % ETA = ETAMXF

          ! Rescale the history array for a change in H by a factor of ETA. ------
          R = ONE

          do J = 2, vstate % L
             R = R * vstate % ETA
             rwork % YH(:,J) = rwork % YH(:,J) * R
          end do
          vstate % H = vstate % HSCAL * vstate % ETA
          vstate % HSCAL = vstate % H
          vstate % RC = vstate % RC * vstate % ETA
          vstate % NQNYH = vstate % NQ*VODE_NEQS

          cycle
       end IF

       ! -----------------------------------------------------------------------
       !  Control reaches this section if 3 or more consecutive failures
       !  have occurred.  It is assumed that the elements of the YH array
       !  have accumulated errors of the wrong order.  The order is reduced
       !  by one, if possible.  Then H is reduced by a factor of 0.1 and
       !  the step is retried.  After a total of 7 consecutive failures,
       !  an exit is taken with KFLAG = -1.
       ! -----------------------------------------------------------------------
       IF (vstate % KFLAG .EQ. KFH) then
          vstate % KFLAG = -1
          vstate % JSTART = 1
          RETURN
       end IF
       IF (vstate % NQ /= 1) then
          vstate % ETA = MAX(ETAMIN,vstate % HMIN/ABS(vstate % H))
          CALL DVJUST (-1, rwork, vstate)
          vstate % L = vstate % NQ
          vstate % NQ = vstate % NQ - 1
          vstate % NQWAIT = vstate % L

          ! Rescale the history array for a change in H by a factor of ETA. ------
          R = ONE

          do J = 2, vstate % L
             R = R * vstate % ETA
             rwork % YH(:,J) = rwork % YH(:,J) * R
          end do
          vstate % H = vstate % HSCAL * vstate % ETA
          vstate % HSCAL = vstate % H
          vstate % RC = vstate % RC * vstate % ETA
          vstate % NQNYH = vstate % NQ*VODE_NEQS

          cycle

       end IF

       vstate % ETA = MAX(ETAMIN,vstate % HMIN/ABS(vstate % H))
       vstate % H = vstate % H * vstate % ETA
       vstate % HSCAL = vstate % H
       vstate % TAU(1) = vstate % H
       CALL f_rhs (vstate % TN, vstate, rwork % savf)
       vstate % NFE = vstate % NFE + 1
       do I = 1, VODE_NEQS
          rwork % yh(I,2) = vstate % H * rwork % savf(I)
       end do
       vstate % NQWAIT = 10

    end do


    ! -----------------------------------------------------------------------
    !  If NQWAIT = 0, an increase or decrease in order by one is considered.
    !  Factors ETAQ, ETAQM1, ETAQP1 are computed by which H could
    !  be multiplied at order q, q-1, or q+1, respectively.
    !  The largest of these is determined, and the new order and
    !  step size set accordingly.
    !  A change of H or NQ is made only if H increases by at least a
    !  factor of THRESH.  If an order change is considered and rejected,
    !  then NQWAIT is set to 2 (reconsider it after 2 steps).
    ! -----------------------------------------------------------------------

    already_set_eta = .false.

    !  Compute ratio of new H to current H at the current order. ------------
    FLOTL = REAL(vstate % L)
    ETAQ = ONE/((BIAS2*DSM)**(ONE/FLOTL) + ADDON)
    IF (vstate % NQWAIT == 0) then
       vstate % NQWAIT = 2
       ETAQM1 = ZERO
       IF (vstate % NQ /= 1) then
          ! Compute ratio of new H to current H at the current order less one. ---
          DDN = DVNORM (rwork % yh(:,vstate % L), rwork % ewt)/vstate % TQ(1)
          ETAQM1 = ONE/((BIAS1*DDN)**(ONE/(FLOTL - ONE)) + ADDON)
       end IF
       ETAQP1 = ZERO
       IF (vstate % L /= VODE_LMAX) then
          ! Compute ratio of new H to current H at current order plus one. -------
          CNQUOT = (vstate % TQ(5)/vstate % CONP)*(vstate % H/vstate % TAU(2))**vstate % L
          do I = 1, VODE_NEQS
             rwork % savf(I) = rwork % acor(I) - CNQUOT * rwork % yh(I,VODE_LMAX)
          end do
          DUP = DVNORM (rwork % savf, rwork % ewt)/vstate % TQ(3)
          ETAQP1 = ONE/((BIAS3*DUP)**(ONE/(FLOTL + ONE)) + ADDON)
       end IF
       IF (ETAQ < ETAQP1) then
          IF (ETAQP1 .GT. ETAQM1) then
             vstate % ETA = ETAQP1
             vstate % NEWQ = vstate % NQ + 1
             rwork % yh(1:VODE_NEQS,VODE_LMAX) = rwork % acor(1:VODE_NEQS)
          else
             vstate % ETA = ETAQM1
             vstate % NEWQ = vstate % NQ - 1
          end if
          already_set_eta = .true.
       end IF

       IF (ETAQ .LT. ETAQM1 .and. .not. already_set_eta) then
          vstate % ETA = ETAQM1
          vstate % NEWQ = vstate % NQ - 1
          already_set_eta = .true.
       end IF
    end IF

    if (.not. already_set_eta) then
       vstate % ETA = ETAQ
       vstate % NEWQ = vstate % NQ
    end if

    ! Test tentative new H against THRESH, ETAMAX, and HMXI, then exit. ----
    IF (vstate % ETA >= THRESH .and. vstate % ETAMAX /= ONE) then
       vstate % ETA = MIN(vstate % ETA,vstate % ETAMAX)
       vstate % ETA = vstate % ETA/MAX(ONE,ABS(vstate % H)*vstate % HMXI * vstate % ETA)
       vstate % NEWH = 1
       vstate % HNEW = vstate % H * vstate % ETA
       vstate % ETAMAX = ETAMX3
       IF (vstate % NST .LE. 10) vstate % ETAMAX = ETAMX2
       R = ONE/vstate % TQ(2)
       rwork % acor(:) = rwork % acor(:) * R
       vstate % JSTART = 1
       RETURN
    end IF

    vstate % NEWQ = vstate % NQ
    vstate % NEWH = 0
    vstate % ETA = ONE
    vstate % HNEW = vstate % H
    vstate % ETAMAX = ETAMX3
    IF (vstate % NST .LE. 10) vstate % ETAMAX = ETAMX2
    R = ONE/vstate % TQ(2)
    rwork % acor(:) = rwork % acor(:) * R
    vstate % JSTART = 1
    RETURN

  end subroutine dvstep

end module cuvode_dvstep_module
