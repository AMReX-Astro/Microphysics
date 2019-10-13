module cuvode_dvstep_module

  use cuvode_parameters_module, only: VODE_LMAX, VODE_NEQS, VODE_LIW,   &
                                      VODE_LENWM, VODE_MAXORD, VODE_ITOL
  use cuvode_types_module, only: dvode_t, rwork_t
  use amrex_fort_module, only: rt => amrex_real
  use blas_module

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
  subroutine dvstep(IWM, rwork, vstate)

    !$acc routine seq
    
    ! -----------------------------------------------------------------------
    !  Call sequence input -- Y, IWM, RPAR, IPAR, rwork, vstate
    !  Call sequence output -- YH, ACOR, WM, IWM
    !  COMMON block variables accessed:
    !      /DVOD01/  ACNRM, EL(13), H, HMIN, HMXI, HNEW, HSCAL, RC, TAU(13),
    !                TQ(5), TN, JCUR, JSTART, KFLAG, KUTH,
    !                L, LMAX, MAXORD, N, NEWQ, NQ, NQWAIT
    !      /DVOD02/  HU, NCFN, NETF, NFE, NQU, NST
    ! 
    !  Subroutines called by DVSTEP: F, DAXPY, DCOPY, DSCAL,
    !                                DVJUST, VNLS, DVSET
    !  Function routines called by DVSTEP: DVNORM
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
    !  WM,IWM = Real and integer work arrays associated with matrix
    !           operations in VNLS.
    !  F      = Dummy name for the user supplied subroutine for f.
    !  JAC    = Dummy name for the user supplied Jacobian subroutine.
    !  PSOL   = Dummy name for the subroutine passed to VNLS, for
    !           possible use there.
    !  VNLS   = Dummy name for the nonlinear system solving subroutine,
    !           whose real name is dependent on the method used.
    !  RPAR, IPAR = Dummy names for user's real and integer work arrays.
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
    integer,       intent(inout) :: IWM(VODE_LIW)

    ! Declare local variables
    real(rt) :: CNQUOT, DDN, DSM, DUP, TOLD
    real(rt) :: ETAQ, ETAQM1, ETAQP1, FLOTL, R
    integer    :: I, I1, I2, IBACK, J, JB, NCF, NFLAG

    ! Parameter declarations
    integer, parameter :: KFC = -3
    integer, parameter :: KFH = -7
    integer, parameter :: MXNCF = 10
    real(rt), parameter :: ADDON = 1.0D-6
    real(rt), parameter :: BIAS1 = 6.0D0
    real(rt), parameter :: BIAS2 = 6.0D0
    real(rt), parameter :: BIAS3 = 10.0D0
    real(rt), parameter :: ETACF = 0.25D0
    real(rt), parameter :: ETAMIN = 0.1D0
    real(rt), parameter :: ETAMXF = 0.2D0
    real(rt), parameter :: ETAMX1 = 1.0D4
    real(rt), parameter :: ETAMX2 = 10.0D0
    real(rt), parameter :: ETAMX3 = 10.0D0
    real(rt), parameter :: ONEPSM = 1.00001D0
    real(rt), parameter :: THRESH = 1.5D0

    !$gpu

    ETAQ   = ONE
    ETAQM1 = ONE

    vstate % KFLAG = 0
    TOLD = vstate % TN
    NCF = 0
    vstate % JCUR = 0
    NFLAG = 0
    IF (vstate % JSTART .GT. 0) GO TO 20
    IF (vstate % JSTART .EQ. -1) GO TO 100
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
    GO TO 200
    ! -----------------------------------------------------------------------
    !  Take preliminary actions on a normal continuation step (JSTART.GT.0).
    !  If the driver changed H, then ETA must be reset and NEWH set to 1.
    !  If a change of order was dictated on the previous step, then
    !  it is done here and appropriate adjustments in the history are made.
    !  On an order decrease, the history array is adjusted by DVJUST.
    !  On an order increase, the history array is augmented by a column.
    !  On a change of step size H, the history array YH is rescaled.
    ! -----------------------------------------------------------------------
20  CONTINUE
    IF (vstate % KUTH .EQ. 1) THEN
       vstate % ETA = MIN(vstate % ETA,vstate % H/vstate % HSCAL)
       vstate % NEWH = 1
    ENDIF
50  IF (vstate % NEWH .EQ. 0) GO TO 200
    IF (vstate % NEWQ .EQ. vstate % NQ) GO TO 150
    IF (vstate % NEWQ .LT. vstate % NQ) THEN
       CALL DVJUST (-1, rwork, vstate)
       vstate % NQ = vstate % NEWQ
       vstate % L = vstate % NQ + 1
       vstate % NQWAIT = vstate % L
       GO TO 150
    ENDIF
    IF (vstate % NEWQ .GT. vstate % NQ) THEN
       CALL DVJUST (1, rwork, vstate)
       vstate % NQ = vstate % NEWQ
       vstate % L = vstate % NQ + 1
       vstate % NQWAIT = vstate % L
       GO TO 150
    ENDIF
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
100 CONTINUE
    !don - remove the following logic we don't use?
    IF (VODE_NEQS .EQ. VODE_NEQS) GO TO 120
    I1 = 1 + (vstate % NEWQ + 1)*VODE_NEQS
    I2 = (VODE_MAXORD + 1)*VODE_NEQS
    IF (I1 .GT. I2) GO TO 120
    rwork % YH(:, vstate % NEWQ + 1:) = ZERO
120 IF (vstate % NEWQ .LE. VODE_MAXORD) GO TO 140
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
140 continue
    IF (vstate % KUTH .EQ. 1) vstate % ETA = MIN(vstate % ETA,ABS(vstate % H/vstate % HSCAL))
    IF (vstate % KUTH .EQ. 0) vstate % ETA = MAX(vstate % ETA,vstate % HMIN/ABS(vstate % HSCAL))
    vstate % ETA = vstate % ETA/MAX(ONE,ABS(vstate % HSCAL)*vstate % HMXI*vstate % ETA)
    vstate % NEWH = 1
    vstate % NQWAIT = vstate % L
    IF (vstate % NEWQ .LE. VODE_MAXORD) GO TO 50
    ! Rescale the history array for a change in H by a factor of ETA. ------
150 continue
    R = ONE

    do J = 2, vstate % L
       R = R * vstate % ETA
       CALL DSCALN (VODE_NEQS, R, rwork % YH(1:VODE_NEQS,J), 1)
    end do
    vstate % H = vstate % HSCAL * vstate % ETA
    vstate % HSCAL = vstate % H
    vstate % RC = vstate % RC * vstate % ETA
    vstate % NQNYH = vstate % NQ*VODE_NEQS
    ! -----------------------------------------------------------------------
    !  This section computes the predicted values by effectively
    !  multiplying the YH array by the Pascal triangle matrix.
    !  DVSET is called to calculate all integration coefficients.
    !  RC is the ratio of new to old values of the coefficient H/EL(2)=h/l1.
    ! -----------------------------------------------------------------------
200 continue
    vstate % TN = vstate % TN + vstate % H

    call advance_nordsieck(rwork, vstate)

    CALL DVSET(vstate)
    vstate % RL1 = ONE/vstate % EL(2)
    vstate % RC = vstate % RC * (vstate % RL1/vstate % PRL1)
    vstate % PRL1 = vstate % RL1
    ! 
    !  Call the nonlinear system solver. ------------------------------------
    !
    CALL dvnlsd (IWM, NFLAG, rwork, vstate)

    IF (NFLAG .EQ. 0) GO TO 450
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

    IF (NFLAG .LT. -1) GO TO 680
    IF (ABS(vstate % H) .LE. vstate % HMIN*ONEPSM) GO TO 670
    IF (NCF .EQ. MXNCF) GO TO 670
    vstate % ETA = ETACF
    vstate % ETA = MAX(vstate % ETA,vstate % HMIN/ABS(vstate % H))
    NFLAG = -1
    GO TO 150
    ! -----------------------------------------------------------------------
    !  The corrector has converged (NFLAG = 0).  The local error test is
    !  made and control passes to statement 500 if it fails.
    ! -----------------------------------------------------------------------
450 CONTINUE
    DSM = vstate % ACNRM/vstate % TQ(2)
    IF (DSM .GT. ONE) GO TO 500
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
       CALL DAXPYN(VODE_NEQS, vstate % EL(J), rwork % acor, 1, rwork % yh(:,J), 1)
    end do
    vstate % NQWAIT = vstate % NQWAIT - 1
    IF ((vstate % L .EQ. VODE_LMAX) .OR. (vstate % NQWAIT .NE. 1)) GO TO 490
    rwork % yh(1:VODE_NEQS,VODE_LMAX) = rwork % acor(1:VODE_NEQS)
    
    vstate % CONP = vstate % TQ(5)
490 IF (vstate % ETAMAX .NE. ONE) GO TO 560
    IF (vstate % NQWAIT .LT. 2) vstate % NQWAIT = 2
    vstate % NEWQ = vstate % NQ
    vstate % NEWH = 0
    vstate % ETA = ONE
    vstate % HNEW = vstate % H
    GO TO 690
    ! -----------------------------------------------------------------------
    !  The error test failed.  KFLAG keeps track of multiple failures.
    !  Restore TN and the YH array to their previous values, and prepare
    !  to try the step again.  Compute the optimum step size for the
    !  same order.  After repeated failures, H is forced to decrease
    !  more rapidly.
    ! -----------------------------------------------------------------------
500 continue
    vstate % KFLAG = vstate % KFLAG - 1
    vstate % NETF = vstate % NETF + 1
    NFLAG = -2
    vstate % TN = TOLD

    call retract_nordsieck(rwork, vstate)

    IF (ABS(vstate % H) .LE. vstate % HMIN*ONEPSM) GO TO 660
    vstate % ETAMAX = ONE
    IF (vstate % KFLAG .LE. KFC) GO TO 530
    ! Compute ratio of new H to current H at the current order. ------------
    FLOTL = REAL(vstate % L)
    vstate % ETA = ONE/((BIAS2*DSM)**(ONE/FLOTL) + ADDON)
    vstate % ETA = MAX(vstate % ETA,vstate % HMIN/ABS(vstate % H),ETAMIN)
    IF ((vstate % KFLAG .LE. -2) .AND. (vstate % ETA .GT. ETAMXF)) vstate % ETA = ETAMXF
    GO TO 150
    ! -----------------------------------------------------------------------
    !  Control reaches this section if 3 or more consecutive failures
    !  have occurred.  It is assumed that the elements of the YH array
    !  have accumulated errors of the wrong order.  The order is reduced
    !  by one, if possible.  Then H is reduced by a factor of 0.1 and
    !  the step is retried.  After a total of 7 consecutive failures,
    !  an exit is taken with KFLAG = -1.
    ! -----------------------------------------------------------------------
530 IF (vstate % KFLAG .EQ. KFH) GO TO 660
    IF (vstate % NQ .EQ. 1) GO TO 540
    vstate % ETA = MAX(ETAMIN,vstate % HMIN/ABS(vstate % H))
    CALL DVJUST (-1, rwork, vstate)
    vstate % L = vstate % NQ
    vstate % NQ = vstate % NQ - 1
    vstate % NQWAIT = vstate % L
    GO TO 150
540 continue
    vstate % ETA = MAX(ETAMIN,vstate % HMIN/ABS(vstate % H))
    vstate % H = vstate % H * vstate % ETA
    vstate % HSCAL = vstate % H
    vstate % TAU(1) = vstate % H
    CALL f_rhs (vstate % TN, vstate % Y, rwork % savf, vstate % RPAR)
    vstate % NFE = vstate % NFE + 1
    do I = 1, VODE_NEQS
       rwork % yh(I,2) = vstate % H * rwork % savf(I)
    end do
    vstate % NQWAIT = 10
    GO TO 200
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
    !  Compute ratio of new H to current H at the current order. ------------
560 continue
    FLOTL = REAL(vstate % L)
    ETAQ = ONE/((BIAS2*DSM)**(ONE/FLOTL) + ADDON)
    IF (vstate % NQWAIT .NE. 0) GO TO 600
    vstate % NQWAIT = 2
    ETAQM1 = ZERO
    IF (vstate % NQ .EQ. 1) GO TO 570
    ! Compute ratio of new H to current H at the current order less one. ---
    DDN = DVNORM (rwork % yh(:,vstate % L), rwork % ewt)/vstate % TQ(1)
    ETAQM1 = ONE/((BIAS1*DDN)**(ONE/(FLOTL - ONE)) + ADDON)
570 continue
    ETAQP1 = ZERO
    IF (vstate % L .EQ. VODE_LMAX) GO TO 580
    ! Compute ratio of new H to current H at current order plus one. -------
    CNQUOT = (vstate % TQ(5)/vstate % CONP)*(vstate % H/vstate % TAU(2))**vstate % L
    do I = 1, VODE_NEQS
       rwork % savf(I) = rwork % acor(I) - CNQUOT * rwork % yh(I,VODE_LMAX)
    end do
    DUP = DVNORM (rwork % savf, rwork % ewt)/vstate % TQ(3)
    ETAQP1 = ONE/((BIAS3*DUP)**(ONE/(FLOTL + ONE)) + ADDON)
580 IF (ETAQ .GE. ETAQP1) GO TO 590
    IF (ETAQP1 .GT. ETAQM1) GO TO 620
    GO TO 610
590 IF (ETAQ .LT. ETAQM1) GO TO 610
600 continue
    vstate % ETA = ETAQ
    vstate % NEWQ = vstate % NQ
    GO TO 630
610 continue
    vstate % ETA = ETAQM1
    vstate % NEWQ = vstate % NQ - 1
    GO TO 630
620 continue
    vstate % ETA = ETAQP1
    vstate % NEWQ = vstate % NQ + 1
    rwork % yh(1:VODE_NEQS,VODE_LMAX) = rwork % acor(1:VODE_NEQS)
    ! Test tentative new H against THRESH, ETAMAX, and HMXI, then exit. ----
630 IF (vstate % ETA .LT. THRESH .OR. vstate % ETAMAX .EQ. ONE) GO TO 640
    vstate % ETA = MIN(vstate % ETA,vstate % ETAMAX)
    vstate % ETA = vstate % ETA/MAX(ONE,ABS(vstate % H)*vstate % HMXI * vstate % ETA)
    vstate % NEWH = 1
    vstate % HNEW = vstate % H * vstate % ETA
    GO TO 690
640 continue
    vstate % NEWQ = vstate % NQ
    vstate % NEWH = 0
    vstate % ETA = ONE
    vstate % HNEW = vstate % H
    GO TO 690
    ! -----------------------------------------------------------------------
    !  All returns are made through this section.
    !  On a successful return, ETAMAX is reset and ACOR is scaled.
    ! -----------------------------------------------------------------------
660 continue
    vstate % KFLAG = -1
    GO TO 720
670 continue
    vstate % KFLAG = -2
    GO TO 720
680 continue
    IF (NFLAG .EQ. -2) vstate % KFLAG = -3
    IF (NFLAG .EQ. -3) vstate % KFLAG = -4
    GO TO 720
690 continue
    vstate % ETAMAX = ETAMX3
    IF (vstate % NST .LE. 10) vstate % ETAMAX = ETAMX2
    R = ONE/vstate % TQ(2)
    CALL DSCALN (VODE_NEQS, R, rwork % acor, 1)

720 continue
    vstate % JSTART = 1
    RETURN
  end subroutine dvstep

end module cuvode_dvstep_module
