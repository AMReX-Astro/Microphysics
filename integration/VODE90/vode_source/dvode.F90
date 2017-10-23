module dvode_module

  use dvode_constants_module
  use vode_rhs_module, only: f_rhs, jac
  use vode_type_module, only: rwork_t  
  use dvode_type_module, only: dvode_t
#ifndef CUDA  
  use dvode_output_module, only: xerrwd
#endif
  use rpar_indices
  use bl_types, only: dp_t
  use blas_module
  use linpack_module
#ifdef CUDA
  use cudafor
#endif
  
  implicit none

  public :: dvode

contains

#ifdef CUDA
  attributes(device) &
#endif
  subroutine dvode(NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK, &
       ISTATE, IOPT, RWORK, LRW, IWORK, LIW, MF, &
       RPAR, IPAR, vstate)       
    !$acc routine seq

    implicit none
       
    type(dvode_t), intent(inout) :: vstate
    
    integer    :: NEQ, ITOL, ITASK, ISTATE, IOPT, LRW, LIW, MF
    integer    :: IWORK(LIW)
    integer    :: IPAR(:)    
    real(dp_t) :: T, TOUT
    real(dp_t) :: Y(NEQ)
    real(dp_t) :: RTOL(NEQ), ATOL(NEQ)
    type(rwork_t) :: RWORK
    real(dp_t) :: RPAR(:)

    logical    :: IHIT
    real(dp_t) :: ATOLI, BIG, EWTI, H0, HMAX, HMX
    real(dp_t) :: RH, RTOLI, SIZE, TCRIT, TNEXT, TOLSF, TP
    integer    :: I, IER, IFLAG, IMXER, JCO, KGO, LENIW, LENJ, LENP, LENRW
    integer    :: LENWM, LF0, MBAND, MFA, ML, MU, NITER
    integer    :: NSLAST
    integer, dimension(2) :: MORD = [12, 5]
#ifndef CUDA    
    character (len=80) :: MSG
#endif
    
    ! Parameter declarations
    integer, parameter :: MXSTP0 = 500
    integer, parameter :: MXHNL0 = 10
    real(dp_t), parameter :: PT2 = 0.2D0
    real(dp_t), parameter :: HUN = 100.0D0

!     ! -----------------------------------------------------------------------
!     !  Block A.
!     !  This code block is executed on every call.
!     !  It tests ISTATE and ITASK for legality and branches appropriately.
!     !  If ISTATE .gt. 1 but the flag INIT shows that initialization has
!     !  not yet been done, an error return occurs.
!     !  If ISTATE = 1 and TOUT = T, return immediately.
!     ! -----------------------------------------------------------------------
!     IF (ISTATE .LT. 1 .OR. ISTATE .GT. 3) GO TO 601
!     IF (ITASK .LT. 1 .OR. ITASK .GT. 5) GO TO 602
!     IF (ISTATE .EQ. 1) GO TO 10
!     IF (vstate % INIT .NE. 1) GO TO 603
!     IF (ISTATE .EQ. 2) GO TO 200
!     GO TO 20
! 10  continue
!     vstate % INIT = 0
!     IF (TOUT .EQ. T) RETURN
    
!     ! -----------------------------------------------------------------------
!     !  Block B.
!     !  The next code block is executed for the initial call (ISTATE = 1),
!     !  or for a continuation call with parameter changes (ISTATE = 3).
!     !  It contains checking of all input and various initializations.
!     !
!     !  First check legality of the non-optional input NEQ, ITOL, IOPT,
!     !  MF, ML, and MU.
!     ! -----------------------------------------------------------------------
! 20  IF (NEQ .LE. 0) GO TO 604
!     IF (ISTATE .EQ. 1) GO TO 25
!     IF (NEQ .GT. vstate % N) GO TO 605
! 25  continue
!     vstate % N = NEQ
!     IF (ITOL .LT. 1 .OR. ITOL .GT. 4) GO TO 606
!     IF (IOPT .LT. 0 .OR. IOPT .GT. 1) GO TO 607
!     vstate % JSV = SIGN(1,MF)
!     MFA = ABS(MF)
!     vstate % METH = MFA/10
!     vstate % MITER = MFA - 10*vstate % METH
!     IF (vstate % METH .LT. 1 .OR. vstate % METH .GT. 2) GO TO 608
!     IF (vstate % MITER .LT. 0 .OR. vstate % MITER .GT. 5) GO TO 608
!     IF (vstate % MITER .LE. 3) GO TO 30
!     ML = IWORK(1)
!     MU = IWORK(2)
!     IF (ML .LT. 0 .OR. ML .GE. vstate % N) GO TO 609
!     IF (MU .LT. 0 .OR. MU .GE. vstate % N) GO TO 610
! 30  CONTINUE

!     ! Next process and check the optional input. ---------------------------
      
!     IF (IOPT .EQ. 1) GO TO 40
!     vstate % MAXORD = MORD(vstate % METH)
!     vstate % MXSTEP = MXSTP0
!     vstate % MXHNIL = MXHNL0
!     IF (ISTATE .EQ. 1) H0 = ZERO
!     vstate % HMXI = ZERO
!     vstate % HMIN = ZERO
!     GO TO 60
! 40  continue
!     vstate % MAXORD = IWORK(5)
!     IF (vstate % MAXORD .LT. 0) GO TO 611
!     IF (vstate % MAXORD .EQ. 0) vstate % MAXORD = 100
!     vstate % MAXORD = MIN(vstate % MAXORD,MORD(vstate % METH))
!     vstate % MXSTEP = IWORK(6)
!     IF (vstate % MXSTEP .LT. 0) GO TO 612
!     IF (vstate % MXSTEP .EQ. 0) vstate % MXSTEP = MXSTP0
!     vstate % MXHNIL = IWORK(7)
!     IF (vstate % MXHNIL .LT. 0) GO TO 613
!     !      EDIT 07/16/2016 -- see comments above about MXHNIL
!     !      IF (MXHNIL .EQ. 0) MXHNIL = MXHNL0
!     IF (ISTATE .NE. 1) GO TO 50
!     H0 = RWORK % CONDOPT(5)
!     IF ((TOUT - T)*H0 .LT. ZERO) GO TO 614
! 50  continue
!     HMAX = RWORK % CONDOPT(6)
!     IF (HMAX .LT. ZERO) GO TO 615
!     vstate % HMXI = ZERO
!     IF (HMAX .GT. ZERO) vstate % HMXI = ONE/HMAX
!     vstate % HMIN = RWORK % CONDOPT(7)
!     IF (vstate % HMIN .LT. ZERO) GO TO 616
    
!     ! -----------------------------------------------------------------------
!     !  Arrays stored in RWORK are denoted  CONDOPT, YH, WM, EWT, SAVF, ACOR.
!     !  Within WM, LOCJS is the location of the saved Jacobian (JSV .gt. 0).
!     ! -----------------------------------------------------------------------
    
! 60  continue
!     vstate % LYH = 21
!     IF (ISTATE .EQ. 1) vstate % NYH = vstate % N
!     vstate % LWM = vstate % LYH + (vstate % MAXORD + 1)*vstate % NYH
!     JCO = MAX(0,vstate % JSV)
!     IF (vstate % MITER .EQ. 0) LENWM = 0
!     IF (vstate % MITER .EQ. 1 .OR. vstate % MITER .EQ. 2) THEN
!        LENWM = 2 + (1 + JCO)*vstate % N*vstate % N
!        vstate % LOCJS = vstate % N*vstate % N + 3
!     ENDIF
!     IF (vstate % MITER .EQ. 3) LENWM = 2 + vstate % N
!     IF (vstate % MITER .EQ. 4 .OR. vstate % MITER .EQ. 5) THEN
!        MBAND = ML + MU + 1
!        LENP = (MBAND + ML)*vstate % N
!        LENJ = MBAND*vstate % N
!        LENWM = 2 + LENP + JCO*LENJ
!        vstate % LOCJS = LENP + 3
!     ENDIF
!     vstate % LMAX = vstate % MAXORD + 1
!     vstate % LEWT = vstate % LWM + LENWM
!     vstate % LSAVF = vstate % LEWT + vstate % N
!     vstate % LACOR = vstate % LSAVF + vstate % N
!     vstate % NEQ = NEQ
!     LENRW = vstate % LACOR + vstate % N - 1
!     IWORK(17) = LENRW
!     vstate % LIWM = 1
!     LENIW = 30 + vstate % N
!     IF (vstate % MITER .EQ. 0 .OR. vstate % MITER .EQ. 3) LENIW = 30
!     IWORK(18) = LENIW
!     IF (LENRW .GT. LRW) GO TO 617
!     IF (LENIW .GT. LIW) GO TO 618
!     ! Check RTOL and ATOL for legality. ------------------------------------
!     RTOLI = RTOL(1)
!     ATOLI = ATOL(1)
!     do I = 1,vstate % N
!        IF (ITOL .GE. 3) RTOLI = RTOL(I)
!        IF (ITOL .EQ. 2 .OR. ITOL .EQ. 4) ATOLI = ATOL(I)
!        IF (RTOLI .LT. ZERO) GO TO 619
!        IF (ATOLI .LT. ZERO) GO TO 620
!     end do
!     IF (ISTATE .EQ. 1) GO TO 100
!     ! If ISTATE = 3, set flag to signal parameter changes to DVSTEP. -------
!     vstate % JSTART = -1
!     IF (vstate % NQ .LE. vstate % MAXORD) GO TO 90
!     ! MAXORD was reduced below NQ.  Copy YH(*,MAXORD+2) into SAVF. ---------
!     CALL DCOPY(vstate % N, rwork % wm, 1, rwork % savf, 1)

!     ! Reload WM(1) = RWORK % wm(1), since LWM may have changed. ---------------
! 90  continue
!     IF (vstate % MITER .GT. 0) rwork % wm(1) = SQRT(vstate % UROUND)
!     GO TO 200

!     ! -----------------------------------------------------------------------
!     !  Block C.
!     !  The next block is for the initial call only (ISTATE = 1).
!     !  It contains all remaining initializations, the initial call to F,
!     !  and the calculation of the initial step size.
!     !  The error weights in EWT are inverted after being loaded.
!     ! -----------------------------------------------------------------------
        
! 100 continue
!     vstate % UROUND = DUMACH()
!     vstate % TN = T
!     IF (ITASK .NE. 4 .AND. ITASK .NE. 5) GO TO 110
!     TCRIT = RWORK % condopt(1)
!     IF ((TCRIT - TOUT)*(TOUT - T) .LT. ZERO) GO TO 625
!     if (H0 .NE. ZERO .AND. (T + H0 - TCRIT)*H0 .GT. ZERO) then
!        H0 = TCRIT - T
!     end if
! 110 continue
!     vstate % JSTART = 0
!     IF (vstate % MITER .GT. 0) RWORK % wm(1) = SQRT(vstate % UROUND)
!     vstate % CCMXJ = PT2
!     vstate % MSBJ = 50
!     vstate % NHNIL = 0
!     vstate % NST = 0
!     vstate % NJE = 0
!     vstate % NNI = 0
!     vstate % NCFN = 0
!     vstate % NETF = 0
!     vstate % NLU = 0
!     vstate % NSLJ = 0
!     NSLAST = 0
!     vstate % HU = ZERO
!     vstate % NQU = 0

!     ! Initial call to F.  (LF0 points to YH(*,2).) -------------------------
!     LF0 = vstate % LYH + vstate % NYH

!     CALL f_rhs (vstate % N, T, Y, rwork % yh(:,2), RPAR, IPAR)
!     vstate % NFE = 1
!     ! Load the initial value vector in YH. ---------------------------------
!     CALL DCOPY(vstate % N, Y, 1, rwork % YH(:,1), 1)

!     ! Load and invert the EWT array.  (H is temporarily set to 1.0.) -------
!     vstate % NQ = 1
!     vstate % H = ONE
!     CALL DEWSET (vstate % N, ITOL, RTOL, ATOL, rwork % YH(:,1), rwork % EWT)
!     do I = 1,vstate % N
!        IF (rwork % ewt(I) .LE. ZERO) GO TO 621
!        rwork % ewt(I) = ONE/rwork % ewt(I)
!     end do
!     IF (H0 .NE. ZERO) GO TO 180

!     ! Call DVHIN to set initial step size H0 to be attempted. --------------
!     CALL DVHIN (vstate % N, T, &
!          rwork % YH, &
!          RPAR, IPAR, TOUT, &
!          vstate % UROUND, &
!          rwork % EWT, &
!          ITOL, ATOL, Y, &
!          rwork % ACOR, &
!          H0, NITER, IER)
!     vstate % NFE = vstate % NFE + NITER
!     IF (IER .NE. 0) GO TO 622
!     ! Adjust H0 if necessary to meet HMAX bound. ---------------------------
! 180 continue
!     RH = ABS(H0)*vstate % HMXI
!     IF (RH .GT. ONE) H0 = H0/RH
!     ! Load H with H0 and scale YH(*,2) by H0. ------------------------------
!     vstate % H = H0
!     CALL DSCAL(vstate % N, H0, rwork % YH(:,2), 1)

!     GO TO 270
    
!     ! -----------------------------------------------------------------------
!     !  Block D.
!     !  The next code block is for continuation calls only (ISTATE = 2 or 3)
!     !  and is to check stop conditions before taking a step.
!     ! -----------------------------------------------------------------------
        
! 200 continue
!     NSLAST = vstate % NST
!     vstate % KUTH = 0

!     GO TO (210, 250, 220, 230, 240), ITASK
! 210 IF ((vstate % TN - TOUT) * vstate % H .LT. ZERO) GO TO 250
!     CALL DVINDY (rwork % yh, TOUT, 0, Y, IFLAG, vstate)
!     IF (IFLAG .NE. 0) GO TO 627
!     T = TOUT
!     GO TO 420
! 220 continue
!     TP = vstate % TN - vstate % HU * (ONE + HUN * vstate % UROUND)
!     IF ((TP - TOUT) * vstate % H .GT. ZERO) GO TO 623
!     IF ((vstate % TN - TOUT) * vstate % H .LT. ZERO) GO TO 250
!     GO TO 400
! 230 continue
!     TCRIT = RWORK % condopt(1)
!     IF ((vstate % TN - TCRIT) * vstate % H .GT. ZERO) GO TO 624
!     IF ((TCRIT - TOUT) * vstate % H .LT. ZERO) GO TO 625
!     IF ((vstate % TN - TOUT) * vstate % H .LT. ZERO) GO TO 245
!     CALL DVINDY (rwork % yh, TOUT, 0, Y, IFLAG, vstate)
!     IF (IFLAG .NE. 0) GO TO 627
!     T = TOUT
!     GO TO 420
! 240 continue
!     TCRIT = RWORK % condopt(1)
!     IF ((vstate % TN - TCRIT) * vstate % H .GT. ZERO) GO TO 624
! 245 continue
!     HMX = ABS(vstate % TN) + ABS(vstate % H)
!     IHIT = ABS(vstate % TN - TCRIT) .LE. HUN * vstate % UROUND * HMX
!     IF (IHIT) GO TO 400
!     TNEXT = vstate % TN + vstate % HNEW*(ONE + FOUR * vstate % UROUND)
!     IF ((TNEXT - TCRIT) * vstate % H .LE. ZERO) GO TO 250
!     vstate % H = (TCRIT - vstate % TN)*(ONE - FOUR * vstate % UROUND)
!     vstate % KUTH = 1
    
!     ! -----------------------------------------------------------------------
!     !  Block E.
!     !  The next block is normally executed for all calls and contains
!     !  the call to the one-step core integrator DVSTEP.
!     !
!     !  This is a looping point for the integration steps.
!     !
!     !  First check for too many steps being taken, update EWT (if not at
!     !  start of problem), check for too much accuracy being requested, and
!     !  check for H below the roundoff level in T.
!     ! -----------------------------------------------------------------------

! 250 CONTINUE
!     IF ((vstate % NST-NSLAST) .GE. vstate % MXSTEP) GO TO 500
!     CALL DEWSET (vstate % N, ITOL, RTOL, ATOL,  rwork % YH(:,1), rwork % EWT)
!     do I = 1,vstate % N
!        IF (rwork % ewt(I) .LE. ZERO) GO TO 510
!        rwork % ewt(I) = ONE/rwork % ewt(I)
!     end do
! 270 continue
!     TOLSF = vstate % UROUND * DVNORM (vstate % N, rwork % YH(:,1), rwork % EWT)
!     IF (TOLSF .LE. ONE) GO TO 280
!     TOLSF = TOLSF*TWO
!     IF (vstate % NST .EQ. 0) GO TO 626
!     GO TO 520
! 280 IF ((vstate % TN + vstate % H) .NE. vstate % TN) GO TO 290
!     vstate % NHNIL = vstate % NHNIL + 1
!     IF (vstate % NHNIL .GT. vstate % MXHNIL) GO TO 290
! #ifndef CUDA    
!     MSG = 'DVODE--  Warning: internal T (=R1) and H (=R2) are'
!     CALL XERRWD (MSG, 50, 101, 1, 0, 0, 0, 0, ZERO, ZERO)
!     MSG='      such that in the machine, T + H = T on the next step  '
!     CALL XERRWD (MSG, 60, 101, 1, 0, 0, 0, 0, ZERO, ZERO)
!     MSG = '      (H = step size). solver will continue anyway'
!     CALL XERRWD (MSG, 50, 101, 1, 0, 0, 0, 2, vstate % TN, vstate % H)
! #endif    
!     IF (vstate % NHNIL .LT. vstate % MXHNIL) GO TO 290
! #ifndef CUDA    
!     MSG = 'DVODE--  Above warning has been issued I1 times.  '
!     CALL XERRWD (MSG, 50, 102, 1, 0, 0, 0, 0, ZERO, ZERO)
!     MSG = '      it will not be issued again for this problem'
!     CALL XERRWD (MSG, 50, 102, 1, 1, vstate % MXHNIL, 0, 0, ZERO, ZERO)
! #endif
! 290 CONTINUE
    
!     ! -----------------------------------------------------------------------
!     !  CALL DVSTEP (Y, YH, NYH, YH, EWT, SAVF, VSAV, ACOR,
!     !               WM, IWM, F, JAC, F, DVNLSD, RPAR, IPAR)
!     ! -----------------------------------------------------------------------
!     CALL DVSTEP(Y, IWORK, RPAR, IPAR, rwork, vstate)
!     KGO = 1 - vstate % KFLAG
!     ! Branch on KFLAG.  Note: In this version, KFLAG can not be set to -3.
!     !  KFLAG .eq. 0,   -1,  -2
!     GO TO (300, 530, 540), KGO
    
!     ! -----------------------------------------------------------------------
!     !  Block F.
!     !  The following block handles the case of a successful return from the
!     !  core integrator (KFLAG = 0).  Test for stop conditions.
!     ! -----------------------------------------------------------------------
    
! 300 continue
!     vstate % INIT = 1
!     vstate % KUTH = 0
!     GO TO (310, 400, 330, 340, 350), ITASK
!     ! ITASK = 1.  If TOUT has been reached, interpolate. -------------------
! 310 IF ((vstate % TN - TOUT) * vstate % H .LT. ZERO) GO TO 250
!     CALL DVINDY (rwork % yh, TOUT, 0, Y, IFLAG, vstate)
!     T = TOUT
!     GO TO 420
!     ! ITASK = 3.  Jump to exit if TOUT was reached. ------------------------
! 330 IF ((vstate % TN - TOUT) * vstate % H .GE. ZERO) GO TO 400
!     GO TO 250
!     ! ITASK = 4.  See if TOUT or TCRIT was reached.  Adjust H if necessary.
! 340 IF ((vstate % TN - TOUT) * vstate % H .LT. ZERO) GO TO 345
!     CALL DVINDY (rwork % yh, TOUT, 0, Y, IFLAG, vstate)
!     T = TOUT
!     GO TO 420
! 345 continue
!     HMX = ABS(vstate % TN) + ABS(vstate % H)
!     IHIT = ABS(vstate % TN - TCRIT) .LE. HUN * vstate % UROUND * HMX
!     IF (IHIT) GO TO 400
!     TNEXT = vstate % TN + vstate % HNEW*(ONE + FOUR * vstate % UROUND)
!     IF ((TNEXT - TCRIT) * vstate % H .LE. ZERO) GO TO 250
!     vstate % H = (TCRIT - vstate % TN)*(ONE - FOUR * vstate % UROUND)
!     vstate % KUTH = 1
!     GO TO 250
!     ! ITASK = 5.  See if TCRIT was reached and jump to exit. ---------------
! 350 continue
!     HMX = ABS(vstate % TN) + ABS(vstate % H)
!     IHIT = ABS(vstate % TN - TCRIT) .LE. HUN * vstate % UROUND * HMX
    
!     ! -----------------------------------------------------------------------
!     !  Block G.
!     !  The following block handles all successful returns from DVODE.
!     !  If ITASK .ne. 1, Y is loaded from YH and T is set accordingly.
!     !  ISTATE is set to 2, and the optional output is loaded into the work
!     !  arrays before returning.
!     ! -----------------------------------------------------------------------
    
! 400 CONTINUE
!     CALL DCOPY(vstate % N, rwork % YH(:,1), 1, Y, 1)

!     T = vstate % TN
!     IF (ITASK .NE. 4 .AND. ITASK .NE. 5) GO TO 420
!     IF (IHIT) T = TCRIT
! 420 continue
!     ISTATE = 2
!     RWORK % condopt(11) = vstate % HU
!     RWORK % condopt(12) = vstate % HNEW
!     RWORK % condopt(13) = vstate % TN
!     IWORK(11) = vstate % NST
!     IWORK(12) = vstate % NFE
!     IWORK(13) = vstate % NJE
!     IWORK(14) = vstate % NQU
!     IWORK(15) = vstate % NEWQ
!     IWORK(19) = vstate % NLU
!     IWORK(20) = vstate % NNI
!     IWORK(21) = vstate % NCFN
!     IWORK(22) = vstate % NETF

!     return
    
!     ! -----------------------------------------------------------------------
!     !  Block H.
!     !  The following block handles all unsuccessful returns other than
!     !  those for illegal input.  First the error message routine is called.
!     !  if there was an error test or convergence test failure, IMXER is set.
!     !  Then Y is loaded from YH, and T is set to TN.
!     !  The optional output is loaded into the work arrays before returning.
!     ! -----------------------------------------------------------------------
    
!     ! The maximum number of steps was taken before reaching TOUT. ----------
! 500 continue
! #ifndef CUDA    
!     MSG = 'DVODE--  At current T (=R1), MXSTEP (=I1) steps   '
!     CALL XERRWD (MSG, 50, 201, 1, 0, 0, 0, 0, ZERO, ZERO)
!     MSG = '      taken on this call before reaching TOUT     '
!     CALL XERRWD (MSG, 50, 201, 1, 1, vstate % MXSTEP, 0, 1, vstate % TN, ZERO)
! #endif    
!     ISTATE = -1
!     GO TO 580
!     ! EWT(i) .le. 0.0 for some i (not at start of problem). ----------------
! 510 continue
! #ifndef CUDA        
!     EWTI = rwork % ewt(I)
!     MSG = 'DVODE--  At T (=R1), EWT(I1) has become R2 .le. 0.'
!     CALL XERRWD (MSG, 50, 202, 1, 1, I, 0, 2, vstate % TN, EWTI)
! #endif
!     ISTATE = -6
!     GO TO 580
!     ! Too much accuracy requested for machine precision. -------------------
! 520 continue
! #ifndef CUDA
!     MSG = 'DVODE--  At T (=R1), too much accuracy requested  '
!     CALL XERRWD (MSG, 50, 203, 1, 0, 0, 0, 0, ZERO, ZERO)
!     MSG = '      for precision of machine:   see TOLSF (=R2) '
!     CALL XERRWD (MSG, 50, 203, 1, 0, 0, 0, 2, vstate % TN, TOLSF)
! #endif
!     RWORK % condopt(14) = TOLSF
!     ISTATE = -2
!     GO TO 580
!     ! KFLAG = -1.  Error test failed repeatedly or with ABS(H) = HMIN. -----
! 530 continue
! #ifndef CUDA        
!     MSG = 'DVODE--  At T(=R1) and step size H(=R2), the error'
!     CALL XERRWD (MSG, 50, 204, 1, 0, 0, 0, 0, ZERO, ZERO)
!     MSG = '      test failed repeatedly or with abs(H) = HMIN'
!     CALL XERRWD (MSG, 50, 204, 1, 0, 0, 0, 2, vstate % TN, vstate % H)
! #endif    
!     ISTATE = -4
!     GO TO 560
!     ! KFLAG = -2.  Convergence failed repeatedly or with ABS(H) = HMIN. ----
! 540 continue
! #ifndef CUDA        
!     MSG = 'DVODE--  At T (=R1) and step size H (=R2), the    '
!     CALL XERRWD (MSG, 50, 205, 1, 0, 0, 0, 0, ZERO, ZERO)
!     MSG = '      corrector convergence failed repeatedly     '
!     CALL XERRWD (MSG, 50, 205, 1, 0, 0, 0, 0, ZERO, ZERO)
!     MSG = '      or with abs(H) = HMIN   '
!     CALL XERRWD (MSG, 30, 205, 1, 0, 0, 0, 2, vstate % TN, vstate % H)
! #endif    
!     ISTATE = -5
!     ! Compute IMXER if relevant. -------------------------------------------
! 560 continue
!     BIG = ZERO
!     IMXER = 1
!     do I = 1,vstate % N
!        SIZE = ABS(rwork % acor(I) * rwork % ewt(I))
!        IF (BIG .GE. SIZE) exit
!        BIG = SIZE
!        IMXER = I
!     end do
!     IWORK(16) = IMXER
!     ! Set Y vector, T, and optional output. --------------------------------
! 580 CONTINUE
!     CALL DCOPY(vstate % N, rwork % YH(:,1), 1, Y, 1)

!     T = vstate % TN
!     RWORK % condopt(11) = vstate % HU
!     RWORK % condopt(12) = vstate % H
!     RWORK % condopt(13) = vstate % TN
!     IWORK(11) = vstate % NST
!     IWORK(12) = vstate % NFE
!     IWORK(13) = vstate % NJE
!     IWORK(14) = vstate % NQU
!     IWORK(15) = vstate % NQ
!     IWORK(19) = vstate % NLU
!     IWORK(20) = vstate % NNI
!     IWORK(21) = vstate % NCFN
!     IWORK(22) = vstate % NETF

!     return
    
!     ! -----------------------------------------------------------------------
!     !  Block I.
!     !  The following block handles all error returns due to illegal input
!     !  (ISTATE = -3), as detected before calling the core integrator.
!     !  First the error message routine is called.   If the illegal input
!     !  is a negative ISTATE, the run is aborted (apparent infinite loop).
!     ! -----------------------------------------------------------------------
    
! 601 continue
! #ifndef CUDA    
!     MSG = 'DVODE--  ISTATE (=I1) illegal '
!     CALL XERRWD (MSG, 30, 1, 1, 1, ISTATE, 0, 0, ZERO, ZERO)
! #endif
!     IF (ISTATE .LT. 0) GO TO 800
!     GO TO 700
! 602 continue
! #ifndef CUDA    
!     MSG = 'DVODE--  ITASK (=I1) illegal  '
!     CALL XERRWD (MSG, 30, 2, 1, 1, ITASK, 0, 0, ZERO, ZERO)
! #endif
!     GO TO 700
! 603 continue
! #ifndef CUDA    
!     MSG='DVODE--  ISTATE (=I1) .gt. 1 but DVODE not initialized      '
!     CALL XERRWD (MSG, 60, 3, 1, 1, ISTATE, 0, 0, ZERO, ZERO)
! #endif    
!     GO TO 700
! 604 continue
! #ifndef CUDA    
!     MSG = 'DVODE--  NEQ (=I1) .lt. 1     '
!     CALL XERRWD (MSG, 30, 4, 1, 1, NEQ, 0, 0, ZERO, ZERO)
! #endif
!     GO TO 700
! 605 continue
! #ifndef CUDA
!     MSG = 'DVODE--  ISTATE = 3 and NEQ increased (I1 to I2)  '
!     CALL XERRWD (MSG, 50, 5, 1, 2, vstate % N, NEQ, 0, ZERO, ZERO)
! #endif
!     GO TO 700
! 606 continue
! #ifndef CUDA    
!     MSG = 'DVODE--  ITOL (=I1) illegal   '
!     CALL XERRWD (MSG, 30, 6, 1, 1, ITOL, 0, 0, ZERO, ZERO)
! #endif
!     GO TO 700
! 607 continue
! #ifndef CUDA    
!     MSG = 'DVODE--  IOPT (=I1) illegal   '
!     CALL XERRWD (MSG, 30, 7, 1, 1, IOPT, 0, 0, ZERO, ZERO)
! #endif
!     GO TO 700
! 608 continue
! #ifndef CUDA    
!     MSG = 'DVODE--  MF (=I1) illegal     '
!     CALL XERRWD (MSG, 30, 8, 1, 1, MF, 0, 0, ZERO, ZERO)
! #endif
!     GO TO 700
! 609 continue
! #ifndef CUDA    
!     MSG = 'DVODE--  ML (=I1) illegal:  .lt.0 or .ge.NEQ (=I2)'
!     CALL XERRWD (MSG, 50, 9, 1, 2, ML, NEQ, 0, ZERO, ZERO)
! #endif
!     GO TO 700
! 610 continue
! #ifndef CUDA    
!     MSG = 'DVODE--  MU (=I1) illegal:  .lt.0 or .ge.NEQ (=I2)'
!     CALL XERRWD (MSG, 50, 10, 1, 2, MU, NEQ, 0, ZERO, ZERO)
! #endif    
!     GO TO 700
! 611 continue
! #ifndef CUDA    
!     MSG = 'DVODE--  MAXORD (=I1) .lt. 0  '
!     CALL XERRWD (MSG, 30, 11, 1, 1, vstate % MAXORD, 0, 0, ZERO, ZERO)
! #endif    
!     GO TO 700
! 612 continue
! #ifndef CUDA    
!     MSG = 'DVODE--  MXSTEP (=I1) .lt. 0  '
!     CALL XERRWD (MSG, 30, 12, 1, 1, vstate % MXSTEP, 0, 0, ZERO, ZERO)
! #endif
!     GO TO 700
! 613 continue
! #ifndef CUDA    
!     MSG = 'DVODE--  MXHNIL (=I1) .lt. 0  '
!     CALL XERRWD (MSG, 30, 13, 1, 1, vstate % MXHNIL, 0, 0, ZERO, ZERO)
! #endif
!     GO TO 700
! 614 continue
! #ifndef CUDA    
!     MSG = 'DVODE--  TOUT (=R1) behind T (=R2)      '
!     CALL XERRWD (MSG, 40, 14, 1, 0, 0, 0, 2, TOUT, T)    
!     MSG = '      integration direction is given by H0 (=R1)  '
!     CALL XERRWD (MSG, 50, 14, 1, 0, 0, 0, 1, H0, ZERO)
! #endif
!     GO TO 700
! 615 continue
! #ifndef CUDA    
!     MSG = 'DVODE--  HMAX (=R1) .lt. 0.0  '
!     CALL XERRWD (MSG, 30, 15, 1, 0, 0, 0, 1, HMAX, ZERO)
! #endif
!     GO TO 700
! 616 continue
! #ifndef CUDA    
!     MSG = 'DVODE--  HMIN (=R1) .lt. 0.0  '
!     CALL XERRWD (MSG, 30, 16, 1, 0, 0, 0, 1, vstate % HMIN, ZERO)
! #endif
!     GO TO 700
! 617 CONTINUE
! #ifndef CUDA    
!     MSG='DVODE--  RWORK length needed, LENRW (=I1), exceeds LRW (=I2)'
!     CALL XERRWD (MSG, 60, 17, 1, 2, LENRW, LRW, 0, ZERO, ZERO)
! #endif
!     GO TO 700
! 618 CONTINUE
! #ifndef CUDA    
!     MSG='DVODE--  IWORK length needed, LENIW (=I1), exceeds LIW (=I2)'
!     CALL XERRWD (MSG, 60, 18, 1, 2, LENIW, LIW, 0, ZERO, ZERO)
! #endif
!     GO TO 700
! 619 continue
! #ifndef CUDA    
!     MSG = 'DVODE--  RTOL(I1) is R1 .lt. 0.0        '
!     CALL XERRWD (MSG, 40, 19, 1, 1, I, 0, 1, RTOLI, ZERO)
! #endif
!     GO TO 700
! 620 continue
! #ifndef CUDA    
!     MSG = 'DVODE--  ATOL(I1) is R1 .lt. 0.0        '
!     CALL XERRWD (MSG, 40, 20, 1, 1, I, 0, 1, ATOLI, ZERO)
! #endif
!     GO TO 700
! 621 continue
! #ifndef CUDA    
!     EWTI = rwork % ewt(I)
!     MSG = 'DVODE--  EWT(I1) is R1 .le. 0.0         '
!     CALL XERRWD (MSG, 40, 21, 1, 1, I, 0, 1, EWTI, ZERO)
! #endif
!     GO TO 700
! 622 CONTINUE
! #ifndef CUDA    
!     MSG='DVODE--  TOUT (=R1) too close to T(=R2) to start integration'
!     CALL XERRWD (MSG, 60, 22, 1, 0, 0, 0, 2, TOUT, T)
! #endif
!     GO TO 700
! 623 CONTINUE
! #ifndef CUDA    
!     MSG='DVODE--  ITASK = I1 and TOUT (=R1) behind TCUR - HU (= R2)  '
!     CALL XERRWD (MSG, 60, 23, 1, 1, ITASK, 0, 2, TOUT, TP)
! #endif
!     GO TO 700
! 624 CONTINUE
! #ifndef CUDA    
!     MSG='DVODE--  ITASK = 4 or 5 and TCRIT (=R1) behind TCUR (=R2)   '
!     CALL XERRWD (MSG, 60, 24, 1, 0, 0, 0, 2, TCRIT, vstate % TN)
! #endif
!     GO TO 700
! 625 CONTINUE
! #ifndef CUDA    
!     MSG='DVODE--  ITASK = 4 or 5 and TCRIT (=R1) behind TOUT (=R2)   '
!     CALL XERRWD (MSG, 60, 25, 1, 0, 0, 0, 2, TCRIT, TOUT)
! #endif
!     GO TO 700
! 626 continue
! #ifndef CUDA    
!     MSG = 'DVODE--  At start of problem, too much accuracy   '
!     CALL XERRWD (MSG, 50, 26, 1, 0, 0, 0, 0, ZERO, ZERO)
!     MSG='      requested for precision of machine:   see TOLSF (=R1) '
!     CALL XERRWD (MSG, 60, 26, 1, 0, 0, 0, 1, TOLSF, ZERO)
! #endif
!     RWORK % condopt(14) = TOLSF
!     GO TO 700
! 627 continue
! #ifndef CUDA
!     MSG='DVODE--  Trouble from DVINDY.  ITASK = I1, TOUT = R1.       '
!     CALL XERRWD (MSG, 60, 27, 1, 1, ITASK, 0, 1, TOUT, ZERO)
! #endif
    
! 700 CONTINUE
!     ISTATE = -3

!     return
    
! 800 continue
! #ifndef CUDA    
!     MSG = 'DVODE--  Run aborted:  apparent infinite loop     '
!     CALL XERRWD (MSG, 50, 303, 2, 0, 0, 0, 0, ZERO, ZERO)
! #endif
!     return
  end subroutine dvode
      
end module dvode_module
