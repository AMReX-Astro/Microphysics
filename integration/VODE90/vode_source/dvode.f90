module dvode_module

  use dvode_type_module, only: dvode_t
  
  implicit none

contains

  subroutine dvode(F, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK, &
       ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, MF, &
       RPAR, IPAR)

    use rpar_indices, only: n_rpar_comps, n_ipar_comps

    external F, JAC

    integer    :: NEQ, ITOL, ITASK, ISTATE, IOPT, LRW, LIW, MF
    integer    :: IWORK(LIW), IPAR(n_ipar_comps)
    real(dp_t) :: T, TOUT
    real(dp_t) :: Y(NEQ), RTOL(NEQ), ATOL(NEQ), RWORK(LRW), RPAR(n_rpar_comps)

    external DVNLSD
    logical    :: IHIT
    real(dp_t) :: ATOLI, BIG, EWTI, H0, HMAX, HMX
    real(dp_t) :: RH, RTOLI, SIZE, TCRIT, TNEXT, TOLSF, TP
    integer    :: I, IER, IFLAG, IMXER, JCO, KGO, LENIW, LENJ, LENP, LENRW
    integer    :: LENWM, LF0, MBAND, MFA, ML, MU, NITER
    integer    :: NSLAST
    integer, dimension(2) :: MORD = [12, 5]
    character (len=80) :: MSG

    ! Function declarations
    real(dp_t) :: DUMACH, DVNORM

    ! Parameter declarations
    integer, parameter :: MXSTP0 = 500
    integer, parameter :: MXHNL0 = 10
    real(dp_t), parameter :: ZERO = 0.0D0
    real(dp_t), parameter :: ONE = 1.0D0
    real(dp_t), parameter :: TWO = 2.0D0
    real(dp_t), parameter :: FOUR = 4.0D0
    real(dp_t), parameter :: PT2 = 0.2D0
    real(dp_t), parameter :: HUN = 100.0D0

    ! -----------------------------------------------------------------------
    !  Block A.
    !  This code block is executed on every call.
    !  It tests ISTATE and ITASK for legality and branches appropriately.
    !  If ISTATE .gt. 1 but the flag INIT shows that initialization has
    !  not yet been done, an error return occurs.
    !  If ISTATE = 1 and TOUT = T, return immediately.
    ! -----------------------------------------------------------------------
    IF (ISTATE .LT. 1 .OR. ISTATE .GT. 3) GO TO 601
    IF (ITASK .LT. 1 .OR. ITASK .GT. 5) GO TO 602
    IF (INIT .NE. 1) GO TO 603

    select case (ISTATE)
    case (1)
       INIT = 0
       IF (TOUT .EQ. T) RETURN
    case (2)
       GO TO 200
    case default
       GO TO 20
    end select
    
    ! -----------------------------------------------------------------------
    !  Block B.
    !  The next code block is executed for the initial call (ISTATE = 1),
    !  or for a continuation call with parameter changes (ISTATE = 3).
    !  It contains checking of all input and various initializations.
    !
    !  First check legality of the non-optional input NEQ, ITOL, IOPT,
    !  MF, ML, and MU.
    ! -----------------------------------------------------------------------
    
20  IF (NEQ .LE. 0) GO TO 604
    IF (ISTATE .EQ. 1) GO TO 25
    IF (NEQ .GT. N) GO TO 605
25  N = NEQ
    IF (ITOL .LT. 1 .OR. ITOL .GT. 4) GO TO 606
    IF (IOPT .LT. 0 .OR. IOPT .GT. 1) GO TO 607
    JSV = SIGN(1,MF)
    MFA = ABS(MF)
    METH = MFA/10
    MITER = MFA - 10*METH
    IF (METH .LT. 1 .OR. METH .GT. 2) GO TO 608
    IF (MITER .LT. 0 .OR. MITER .GT. 5) GO TO 608
    IF (MITER .LE. 3) GO TO 30
    ML = IWORK(1)
    MU = IWORK(2)
    IF (ML .LT. 0 .OR. ML .GE. N) GO TO 609
    IF (MU .LT. 0 .OR. MU .GE. N) GO TO 610
30  CONTINUE

    ! Next process and check the optional input. ---------------------------
      
    IF (IOPT .EQ. 1) GO TO 40
    MAXORD = MORD(METH)
    MXSTEP = MXSTP0
    MXHNIL = MXHNL0
    IF (ISTATE .EQ. 1) H0 = ZERO
    HMXI = ZERO
    HMIN = ZERO
    GO TO 60
40  MAXORD = IWORK(5)
    IF (MAXORD .LT. 0) GO TO 611
    IF (MAXORD .EQ. 0) MAXORD = 100
    MAXORD = MIN(MAXORD,MORD(METH))
    MXSTEP = IWORK(6)
    IF (MXSTEP .LT. 0) GO TO 612
    IF (MXSTEP .EQ. 0) MXSTEP = MXSTP0
    MXHNIL = IWORK(7)
    IF (MXHNIL .LT. 0) GO TO 613
    !      EDIT 07/16/2016 -- see comments above about MXHNIL
    !      IF (MXHNIL .EQ. 0) MXHNIL = MXHNL0
    IF (ISTATE .NE. 1) GO TO 50
    H0 = RWORK(5)
    IF ((TOUT - T)*H0 .LT. ZERO) GO TO 614
50  HMAX = RWORK(6)
    IF (HMAX .LT. ZERO) GO TO 615
    HMXI = ZERO
    IF (HMAX .GT. ZERO) HMXI = ONE/HMAX
    HMIN = RWORK(7)
    IF (HMIN .LT. ZERO) GO TO 616
    
    ! -----------------------------------------------------------------------
    !  Set work array pointers and check lengths LRW and LIW.
    !  Pointers to segments of RWORK and IWORK are named by prefixing L to
    !  the name of the segment.  E.g., the segment YH starts at RWORK(LYH).
    !  Segments of RWORK (in order) are denoted  YH, WM, EWT, SAVF, ACOR.
    !  Within WM, LOCJS is the location of the saved Jacobian (JSV .gt. 0).
    ! -----------------------------------------------------------------------
    
60  LYH = 21
    IF (ISTATE .EQ. 1) NYH = N
    LWM = LYH + (MAXORD + 1)*NYH
    JCO = MAX(0,JSV)
    IF (MITER .EQ. 0) LENWM = 0
    IF (MITER .EQ. 1 .OR. MITER .EQ. 2) THEN
       LENWM = 2 + (1 + JCO)*N*N
       LOCJS = N*N + 3
    ENDIF
    IF (MITER .EQ. 3) LENWM = 2 + N
    IF (MITER .EQ. 4 .OR. MITER .EQ. 5) THEN
       MBAND = ML + MU + 1
       LENP = (MBAND + ML)*N
       LENJ = MBAND*N
       LENWM = 2 + LENP + JCO*LENJ
       LOCJS = LENP + 3
    ENDIF
    LEWT = LWM + LENWM
    LSAVF = LEWT + N
    LACOR = LSAVF + N
    LENRW = LACOR + N - 1
    IWORK(17) = LENRW
    LIWM = 1
    LENIW = 30 + N
    IF (MITER .EQ. 0 .OR. MITER .EQ. 3) LENIW = 30
    IWORK(18) = LENIW
    IF (LENRW .GT. LRW) GO TO 617
    IF (LENIW .GT. LIW) GO TO 618
    ! Check RTOL and ATOL for legality. ------------------------------------
    RTOLI = RTOL(1)
    ATOLI = ATOL(1)
    do I = 1,N
       IF (ITOL .GE. 3) RTOLI = RTOL(I)
       IF (ITOL .EQ. 2 .OR. ITOL .EQ. 4) ATOLI = ATOL(I)
       IF (RTOLI .LT. ZERO) GO TO 619
       IF (ATOLI .LT. ZERO) GO TO 620
    end do
    IF (ISTATE .EQ. 1) GO TO 100
    ! If ISTATE = 3, set flag to signal parameter changes to DVSTEP. -------
    JSTART = -1
    IF (NQ .LE. MAXORD) GO TO 90
    ! MAXORD was reduced below NQ.  Copy YH(*,MAXORD+2) into SAVF. ---------
    CALL DCOPY (N, RWORK(LWM), 1, RWORK(LSAVF), 1)
    ! Reload WM(1) = RWORK(LWM), since LWM may have changed. ---------------
90  IF (MITER .GT. 0) RWORK(LWM) = SQRT(UROUND)
    GO TO 200

    ! -----------------------------------------------------------------------
    !  Block C.
    !  The next block is for the initial call only (ISTATE = 1).
    !  It contains all remaining initializations, the initial call to F,
    !  and the calculation of the initial step size.
    !  The error weights in EWT are inverted after being loaded.
    ! -----------------------------------------------------------------------
       
100 UROUND = DUMACH()
    TN = T
    IF (ITASK .NE. 4 .AND. ITASK .NE. 5) GO TO 110
    TCRIT = RWORK(1)
    IF ((TCRIT - TOUT)*(TOUT - T) .LT. ZERO) GO TO 625
    if (H0 .NE. ZERO .AND. (T + H0 - TCRIT)*H0 .GT. ZERO) then
       H0 = TCRIT - T
    end if
110 JSTART = 0
    IF (MITER .GT. 0) RWORK(LWM) = SQRT(UROUND)
    CCMXJ = PT2
    MSBJ = 50
    NHNIL = 0
    NST = 0
    NJE = 0
    NNI = 0
    NCFN = 0
    NETF = 0
    NLU = 0
    NSLJ = 0
    NSLAST = 0
    HU = ZERO
    NQU = 0
    ! Initial call to F.  (LF0 points to YH(*,2).) -------------------------
    LF0 = LYH + NYH
    CALL F (N, T, Y, RWORK(LF0), RPAR, IPAR)
    NFE = 1
    ! Load the initial value vector in YH. ---------------------------------
    CALL DCOPY (N, Y, 1, RWORK(LYH), 1)
    ! Load and invert the EWT array.  (H is temporarily set to 1.0.) -------
    NQ = 1
    H = ONE
    CALL DEWSET (N, ITOL, RTOL, ATOL, RWORK(LYH), RWORK(LEWT))
    do I = 1,N
       IF (RWORK(I+LEWT-1) .LE. ZERO) GO TO 621
       RWORK(I+LEWT-1) = ONE/RWORK(I+LEWT-1)
    end do
    IF (H0 .NE. ZERO) GO TO 180
    ! Call DVHIN to set initial step size H0 to be attempted. --------------
    CALL DVHIN (N, T, RWORK(LYH), RWORK(LF0), F, RPAR, IPAR, TOUT, &
         UROUND, RWORK(LEWT), ITOL, ATOL, Y, RWORK(LACOR), H0, &
         NITER, IER)
    NFE = NFE + NITER
    IF (IER .NE. 0) GO TO 622
    ! Adjust H0 if necessary to meet HMAX bound. ---------------------------
180 RH = ABS(H0)*HMXI
    IF (RH .GT. ONE) H0 = H0/RH
    ! Load H with H0 and scale YH(*,2) by H0. ------------------------------
    H = H0
    CALL DSCAL (N, H0, RWORK(LF0), 1)
    GO TO 270
    
    ! -----------------------------------------------------------------------
    !  Block D.
    !  The next code block is for continuation calls only (ISTATE = 2 or 3)
    !  and is to check stop conditions before taking a step.
    ! -----------------------------------------------------------------------
    
200 NSLAST = NST
    KUTH = 0
    GO TO (210, 250, 220, 230, 240), ITASK
210 IF ((TN - TOUT)*H .LT. ZERO) GO TO 250
    CALL DVINDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
    IF (IFLAG .NE. 0) GO TO 627
    T = TOUT
    GO TO 420
220 TP = TN - HU*(ONE + HUN*UROUND)
    IF ((TP - TOUT)*H .GT. ZERO) GO TO 623
    IF ((TN - TOUT)*H .LT. ZERO) GO TO 250
    GO TO 400
230 TCRIT = RWORK(1)
    IF ((TN - TCRIT)*H .GT. ZERO) GO TO 624
    IF ((TCRIT - TOUT)*H .LT. ZERO) GO TO 625
    IF ((TN - TOUT)*H .LT. ZERO) GO TO 245
    CALL DVINDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
    IF (IFLAG .NE. 0) GO TO 627
    T = TOUT
    GO TO 420
240 TCRIT = RWORK(1)
    IF ((TN - TCRIT)*H .GT. ZERO) GO TO 624
245 HMX = ABS(TN) + ABS(H)
    IHIT = ABS(TN - TCRIT) .LE. HUN*UROUND*HMX
    IF (IHIT) GO TO 400
    TNEXT = TN + HNEW*(ONE + FOUR*UROUND)
    IF ((TNEXT - TCRIT)*H .LE. ZERO) GO TO 250
    H = (TCRIT - TN)*(ONE - FOUR*UROUND)
    KUTH = 1
    
! -----------------------------------------------------------------------
!  Block E.
!  The next block is normally executed for all calls and contains
!  the call to the one-step core integrator DVSTEP.
!
!  This is a looping point for the integration steps.
!
!  First check for too many steps being taken, update EWT (if not at
!  start of problem), check for too much accuracy being requested, and
!  check for H below the roundoff level in T.
! -----------------------------------------------------------------------

250 CONTINUE
    IF ((NST-NSLAST) .GE. MXSTEP) GO TO 500
    CALL DEWSET (N, ITOL, RTOL, ATOL, RWORK(LYH), RWORK(LEWT))
    do I = 1,N
       IF (RWORK(I+LEWT-1) .LE. ZERO) GO TO 510
       RWORK(I+LEWT-1) = ONE/RWORK(I+LEWT-1)
    end do
270 TOLSF = UROUND*DVNORM (N, RWORK(LYH), RWORK(LEWT))
    IF (TOLSF .LE. ONE) GO TO 280
    TOLSF = TOLSF*TWO
    IF (NST .EQ. 0) GO TO 626
    GO TO 520
280 IF ((TN + H) .NE. TN) GO TO 290
    NHNIL = NHNIL + 1
    IF (NHNIL .GT. MXHNIL) GO TO 290
    MSG = 'DVODE--  Warning: internal T (=R1) and H (=R2) are'
    CALL XERRWD (MSG, 50, 101, 1, 0, 0, 0, 0, ZERO, ZERO)
    MSG='      such that in the machine, T + H = T on the next step  '
    CALL XERRWD (MSG, 60, 101, 1, 0, 0, 0, 0, ZERO, ZERO)
    MSG = '      (H = step size). solver will continue anyway'
    CALL XERRWD (MSG, 50, 101, 1, 0, 0, 0, 2, TN, H)
    IF (NHNIL .LT. MXHNIL) GO TO 290
    MSG = 'DVODE--  Above warning has been issued I1 times.  '
    CALL XERRWD (MSG, 50, 102, 1, 0, 0, 0, 0, ZERO, ZERO)
    MSG = '      it will not be issued again for this problem'
    CALL XERRWD (MSG, 50, 102, 1, 1, MXHNIL, 0, 0, ZERO, ZERO)
290 CONTINUE
    
    ! -----------------------------------------------------------------------
    !  CALL DVSTEP (Y, YH, NYH, YH, EWT, SAVF, VSAV, ACOR,
    !               WM, IWM, F, JAC, F, DVNLSD, RPAR, IPAR)
    ! -----------------------------------------------------------------------
    
    CALL DVSTEP (Y, RWORK(LYH), NYH, RWORK(LYH), RWORK(LEWT), &
         RWORK(LSAVF), Y, RWORK(LACOR), RWORK(LWM), IWORK(LIWM), &
         F, JAC, F, DVNLSD, RPAR, IPAR)
    KGO = 1 - KFLAG
    ! Branch on KFLAG.  Note: In this version, KFLAG can not be set to -3.
    !  KFLAG .eq. 0,   -1,  -2
    GO TO (300, 530, 540), KGO
    
    ! -----------------------------------------------------------------------
    !  Block F.
    !  The following block handles the case of a successful return from the
    !  core integrator (KFLAG = 0).  Test for stop conditions.
    ! -----------------------------------------------------------------------
    
300 INIT = 1
    KUTH = 0
    GO TO (310, 400, 330, 340, 350), ITASK
    ! ITASK = 1.  If TOUT has been reached, interpolate. -------------------
310 IF ((TN - TOUT)*H .LT. ZERO) GO TO 250
    CALL DVINDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
    T = TOUT
    GO TO 420
    ! ITASK = 3.  Jump to exit if TOUT was reached. ------------------------
330 IF ((TN - TOUT)*H .GE. ZERO) GO TO 400
    GO TO 250
    ! ITASK = 4.  See if TOUT or TCRIT was reached.  Adjust H if necessary.
340 IF ((TN - TOUT)*H .LT. ZERO) GO TO 345
    CALL DVINDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
    T = TOUT
    GO TO 420
345 HMX = ABS(TN) + ABS(H)
    IHIT = ABS(TN - TCRIT) .LE. HUN*UROUND*HMX
    IF (IHIT) GO TO 400
    TNEXT = TN + HNEW*(ONE + FOUR*UROUND)
    IF ((TNEXT - TCRIT)*H .LE. ZERO) GO TO 250
    H = (TCRIT - TN)*(ONE - FOUR*UROUND)
    KUTH = 1
    GO TO 250
    ! ITASK = 5.  See if TCRIT was reached and jump to exit. ---------------
350 HMX = ABS(TN) + ABS(H)
    IHIT = ABS(TN - TCRIT) .LE. HUN*UROUND*HMX
    
    ! -----------------------------------------------------------------------
    !  Block G.
    !  The following block handles all successful returns from DVODE.
    !  If ITASK .ne. 1, Y is loaded from YH and T is set accordingly.
    !  ISTATE is set to 2, and the optional output is loaded into the work
    !  arrays before returning.
    ! -----------------------------------------------------------------------
    
400 CONTINUE
    CALL DCOPY (N, RWORK(LYH), 1, Y, 1)
    T = TN
    IF (ITASK .NE. 4 .AND. ITASK .NE. 5) GO TO 420
    IF (IHIT) T = TCRIT
420 ISTATE = 2
    RWORK(11) = HU
    RWORK(12) = HNEW
    RWORK(13) = TN
    IWORK(11) = NST
    IWORK(12) = NFE
    IWORK(13) = NJE
    IWORK(14) = NQU
    IWORK(15) = NEWQ
    IWORK(19) = NLU
    IWORK(20) = NNI
    IWORK(21) = NCFN
    IWORK(22) = NETF
    RETURN
    
    ! -----------------------------------------------------------------------
    !  Block H.
    !  The following block handles all unsuccessful returns other than
    !  those for illegal input.  First the error message routine is called.
    !  if there was an error test or convergence test failure, IMXER is set.
    !  Then Y is loaded from YH, and T is set to TN.
    !  The optional output is loaded into the work arrays before returning.
    ! -----------------------------------------------------------------------
    
    ! The maximum number of steps was taken before reaching TOUT. ----------
500 MSG = 'DVODE--  At current T (=R1), MXSTEP (=I1) steps   '
    CALL XERRWD (MSG, 50, 201, 1, 0, 0, 0, 0, ZERO, ZERO)
    MSG = '      taken on this call before reaching TOUT     '
    CALL XERRWD (MSG, 50, 201, 1, 1, MXSTEP, 0, 1, TN, ZERO)
    ISTATE = -1
    GO TO 580
    ! EWT(i) .le. 0.0 for some i (not at start of problem). ----------------
510 EWTI = RWORK(LEWT+I-1)
    MSG = 'DVODE--  At T (=R1), EWT(I1) has become R2 .le. 0.'
    CALL XERRWD (MSG, 50, 202, 1, 1, I, 0, 2, TN, EWTI)
    ISTATE = -6
    GO TO 580
    ! Too much accuracy requested for machine precision. -------------------
520 MSG = 'DVODE--  At T (=R1), too much accuracy requested  '
    CALL XERRWD (MSG, 50, 203, 1, 0, 0, 0, 0, ZERO, ZERO)
    MSG = '      for precision of machine:   see TOLSF (=R2) '
    CALL XERRWD (MSG, 50, 203, 1, 0, 0, 0, 2, TN, TOLSF)
    RWORK(14) = TOLSF
    ISTATE = -2
    GO TO 580
    ! KFLAG = -1.  Error test failed repeatedly or with ABS(H) = HMIN. -----
530 MSG = 'DVODE--  At T(=R1) and step size H(=R2), the error'
    CALL XERRWD (MSG, 50, 204, 1, 0, 0, 0, 0, ZERO, ZERO)
    MSG = '      test failed repeatedly or with abs(H) = HMIN'
    CALL XERRWD (MSG, 50, 204, 1, 0, 0, 0, 2, TN, H)
    ISTATE = -4
    GO TO 560
    ! KFLAG = -2.  Convergence failed repeatedly or with ABS(H) = HMIN. ----
540 MSG = 'DVODE--  At T (=R1) and step size H (=R2), the    '
    CALL XERRWD (MSG, 50, 205, 1, 0, 0, 0, 0, ZERO, ZERO)
    MSG = '      corrector convergence failed repeatedly     '
    CALL XERRWD (MSG, 50, 205, 1, 0, 0, 0, 0, ZERO, ZERO)
    MSG = '      or with abs(H) = HMIN   '
    CALL XERRWD (MSG, 30, 205, 1, 0, 0, 0, 2, TN, H)
    ISTATE = -5
    ! Compute IMXER if relevant. -------------------------------------------
560 BIG = ZERO
    IMXER = 1
    do I = 1,N
       SIZE = ABS(RWORK(I+LACOR-1)*RWORK(I+LEWT-1))
       IF (BIG .GE. SIZE) exit
       BIG = SIZE
       IMXER = I
    end do
    IWORK(16) = IMXER
    ! Set Y vector, T, and optional output. --------------------------------
580 CONTINUE
    CALL DCOPY (N, RWORK(LYH), 1, Y, 1)
    T = TN
    RWORK(11) = HU
    RWORK(12) = H
    RWORK(13) = TN
    IWORK(11) = NST
    IWORK(12) = NFE
    IWORK(13) = NJE
    IWORK(14) = NQU
    IWORK(15) = NQ
    IWORK(19) = NLU
    IWORK(20) = NNI
    IWORK(21) = NCFN
    IWORK(22) = NETF
    RETURN
    
! -----------------------------------------------------------------------
!  Block I.
!  The following block handles all error returns due to illegal input
!  (ISTATE = -3), as detected before calling the core integrator.
!  First the error message routine is called.   If the illegal input
!  is a negative ISTATE, the run is aborted (apparent infinite loop).
! -----------------------------------------------------------------------

601 MSG = 'DVODE--  ISTATE (=I1) illegal '
    CALL XERRWD (MSG, 30, 1, 1, 1, ISTATE, 0, 0, ZERO, ZERO)
    IF (ISTATE .LT. 0) GO TO 800
    GO TO 700
602 MSG = 'DVODE--  ITASK (=I1) illegal  '
    CALL XERRWD (MSG, 30, 2, 1, 1, ITASK, 0, 0, ZERO, ZERO)
    GO TO 700
603 MSG='DVODE--  ISTATE (=I1) .gt. 1 but DVODE not initialized      '
    CALL XERRWD (MSG, 60, 3, 1, 1, ISTATE, 0, 0, ZERO, ZERO)
    GO TO 700
604 MSG = 'DVODE--  NEQ (=I1) .lt. 1     '
    CALL XERRWD (MSG, 30, 4, 1, 1, NEQ, 0, 0, ZERO, ZERO)
    GO TO 700
605 MSG = 'DVODE--  ISTATE = 3 and NEQ increased (I1 to I2)  '
    CALL XERRWD (MSG, 50, 5, 1, 2, N, NEQ, 0, ZERO, ZERO)
    GO TO 700
606 MSG = 'DVODE--  ITOL (=I1) illegal   '
    CALL XERRWD (MSG, 30, 6, 1, 1, ITOL, 0, 0, ZERO, ZERO)
    GO TO 700
607 MSG = 'DVODE--  IOPT (=I1) illegal   '
    CALL XERRWD (MSG, 30, 7, 1, 1, IOPT, 0, 0, ZERO, ZERO)
    GO TO 700
608 MSG = 'DVODE--  MF (=I1) illegal     '
    CALL XERRWD (MSG, 30, 8, 1, 1, MF, 0, 0, ZERO, ZERO)
    GO TO 700
609 MSG = 'DVODE--  ML (=I1) illegal:  .lt.0 or .ge.NEQ (=I2)'
    CALL XERRWD (MSG, 50, 9, 1, 2, ML, NEQ, 0, ZERO, ZERO)
    GO TO 700
610 MSG = 'DVODE--  MU (=I1) illegal:  .lt.0 or .ge.NEQ (=I2)'
    CALL XERRWD (MSG, 50, 10, 1, 2, MU, NEQ, 0, ZERO, ZERO)
    GO TO 700
611 MSG = 'DVODE--  MAXORD (=I1) .lt. 0  '
    CALL XERRWD (MSG, 30, 11, 1, 1, MAXORD, 0, 0, ZERO, ZERO)
    GO TO 700
612 MSG = 'DVODE--  MXSTEP (=I1) .lt. 0  '
    CALL XERRWD (MSG, 30, 12, 1, 1, MXSTEP, 0, 0, ZERO, ZERO)
    GO TO 700
613 MSG = 'DVODE--  MXHNIL (=I1) .lt. 0  '
    CALL XERRWD (MSG, 30, 13, 1, 1, MXHNIL, 0, 0, ZERO, ZERO)
    GO TO 700
614 MSG = 'DVODE--  TOUT (=R1) behind T (=R2)      '
    CALL XERRWD (MSG, 40, 14, 1, 0, 0, 0, 2, TOUT, T)
    MSG = '      integration direction is given by H0 (=R1)  '
    CALL XERRWD (MSG, 50, 14, 1, 0, 0, 0, 1, H0, ZERO)
    GO TO 700
615 MSG = 'DVODE--  HMAX (=R1) .lt. 0.0  '
    CALL XERRWD (MSG, 30, 15, 1, 0, 0, 0, 1, HMAX, ZERO)
    GO TO 700
616 MSG = 'DVODE--  HMIN (=R1) .lt. 0.0  '
    CALL XERRWD (MSG, 30, 16, 1, 0, 0, 0, 1, HMIN, ZERO)
    GO TO 700
617 CONTINUE
    MSG='DVODE--  RWORK length needed, LENRW (=I1), exceeds LRW (=I2)'
    CALL XERRWD (MSG, 60, 17, 1, 2, LENRW, LRW, 0, ZERO, ZERO)
    GO TO 700
618 CONTINUE
    MSG='DVODE--  IWORK length needed, LENIW (=I1), exceeds LIW (=I2)'
    CALL XERRWD (MSG, 60, 18, 1, 2, LENIW, LIW, 0, ZERO, ZERO)
    GO TO 700
619 MSG = 'DVODE--  RTOL(I1) is R1 .lt. 0.0        '
    CALL XERRWD (MSG, 40, 19, 1, 1, I, 0, 1, RTOLI, ZERO)
    GO TO 700
620 MSG = 'DVODE--  ATOL(I1) is R1 .lt. 0.0        '
    CALL XERRWD (MSG, 40, 20, 1, 1, I, 0, 1, ATOLI, ZERO)
    GO TO 700
621 EWTI = RWORK(LEWT+I-1)
    MSG = 'DVODE--  EWT(I1) is R1 .le. 0.0         '
    CALL XERRWD (MSG, 40, 21, 1, 1, I, 0, 1, EWTI, ZERO)
    GO TO 700
622 CONTINUE
    MSG='DVODE--  TOUT (=R1) too close to T(=R2) to start integration'
    CALL XERRWD (MSG, 60, 22, 1, 0, 0, 0, 2, TOUT, T)
    GO TO 700
623 CONTINUE
    MSG='DVODE--  ITASK = I1 and TOUT (=R1) behind TCUR - HU (= R2)  '
    CALL XERRWD (MSG, 60, 23, 1, 1, ITASK, 0, 2, TOUT, TP)
    GO TO 700
624 CONTINUE
    MSG='DVODE--  ITASK = 4 or 5 and TCRIT (=R1) behind TCUR (=R2)   '
    CALL XERRWD (MSG, 60, 24, 1, 0, 0, 0, 2, TCRIT, TN)
    GO TO 700
625 CONTINUE
    MSG='DVODE--  ITASK = 4 or 5 and TCRIT (=R1) behind TOUT (=R2)   '
    CALL XERRWD (MSG, 60, 25, 1, 0, 0, 0, 2, TCRIT, TOUT)
    GO TO 700
626 MSG = 'DVODE--  At start of problem, too much accuracy   '
    CALL XERRWD (MSG, 50, 26, 1, 0, 0, 0, 0, ZERO, ZERO)
    MSG='      requested for precision of machine:   see TOLSF (=R1) '
    CALL XERRWD (MSG, 60, 26, 1, 0, 0, 0, 1, TOLSF, ZERO)
    RWORK(14) = TOLSF
    GO TO 700
627 MSG='DVODE--  Trouble from DVINDY.  ITASK = I1, TOUT = R1.       '
    CALL XERRWD (MSG, 60, 27, 1, 1, ITASK, 0, 1, TOUT, ZERO)
    
700 CONTINUE
    ISTATE = -3
    RETURN
    
800 MSG = 'DVODE--  Run aborted:  apparent infinite loop     '
    CALL XERRWD (MSG, 50, 303, 2, 0, 0, 0, 0, ZERO, ZERO)
    RETURN
  end subroutine dvode
    
end module dvode_module
