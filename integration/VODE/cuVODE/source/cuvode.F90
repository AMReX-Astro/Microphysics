module cuvode_module

  use cuvode_parameters_module, only: VODE_LMAX, VODE_NEQS, VODE_LIW,   &
                                      VODE_LENWM, VODE_MAXORD, VODE_ITOL
  use cuvode_types_module, only: dvode_t, rwork_t
  use vode_rpar_indices
  use amrex_fort_module, only: rt => amrex_real
  use linpack_module
#ifdef AMREX_USE_CUDA
  use cudafor
#endif

  use cuvode_constants_module
  use cuvode_dewset_module
  use cuvode_dvhin_module
  use cuvode_dvindy_module
  use cuvode_dvstep_module
  
  implicit none

  public :: dvode
  
contains

#if defined(AMREX_USE_CUDA) && !defined(AMREX_USE_GPU_PRAGMA)
  attributes(device) &
#endif
  subroutine dvode(vstate, rwork, IWORK, ITASK, IOPT, MF)

    !$acc routine seq

#ifndef AMREX_USE_CUDA
    use cuvode_output_module, only: xerrwd
#endif
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
    integer,       intent(inout) :: IWORK(VODE_LIW)
    integer,       intent(in   ) :: ITASK, IOPT, MF

    ! Declare local variables
    logical    :: IHIT
    real(rt) :: ATOLI, BIG, EWTI, H0, HMAX, HMX
    real(rt) :: RH, RTOLI, SIZE, TCRIT, TNEXT, TOLSF, TP
    integer    :: I, IER, IFLAG, IMXER, JCO, KGO, LENJ, LENP
    integer    :: MBAND, MFA, ML, MU, NITER
    integer    :: NSLAST
#ifndef AMREX_USE_CUDA
    character (len=80) :: MSG
#endif

    ! Parameter declarations
    integer, parameter :: MXSTP0 = 500
    integer, parameter :: MXHNL0 = 10
    real(rt), parameter :: PT2 = 0.2e0_rt

    !$gpu

    ! -----------------------------------------------------------------------
    !  Block A.
    !  This code block is executed on every call.
    !  It tests vstate % ISTATE and ITASK for legality and branches appropriately.
    !  If vstate % ISTATE .gt. 1 but the flag INIT shows that initialization has
    !  not yet been done, an error return occurs.
    !  If vstate % ISTATE = 1 and TOUT = T, return immediately.
    ! -----------------------------------------------------------------------
    if (vstate % ISTATE .LT. 1 .OR. vstate % ISTATE .GT. 3) then
#ifndef AMREX_USE_CUDA
       MSG = 'DVODE--  vstate % ISTATE (=I1) illegal '
       CALL XERRWD (MSG, 30, 1, 1, 1, vstate % ISTATE, 0, 0, ZERO, ZERO)
#endif
       if (vstate % ISTATE .LT. 0) then
#ifndef AMREX_USE_CUDA
          MSG = 'DVODE--  Run aborted:  apparent infinite loop     '
          CALL XERRWD (MSG, 50, 303, 2, 0, 0, 0, 0, ZERO, ZERO)
#endif
       else
          vstate % ISTATE = -3
       end if

       return
    end if

    if (ITASK .LT. 1 .OR. ITASK .GT. 5) then
#ifndef AMREX_USE_CUDA
       MSG = 'DVODE--  ITASK (=I1) illegal  '
       CALL XERRWD (MSG, 30, 2, 1, 1, ITASK, 0, 0, ZERO, ZERO)
#endif
       vstate % ISTATE = -3
       return
    end if

    IF (vstate % ISTATE .EQ. 1) GO TO 10

    if (vstate % INIT .NE. 1) then
#ifndef AMREX_USE_CUDA
       MSG='DVODE--  vstate % ISTATE (=I1) .gt. 1 but DVODE not initialized      '
       CALL XERRWD (MSG, 60, 3, 1, 1, vstate % ISTATE, 0, 0, ZERO, ZERO)
#endif
       vstate % ISTATE = -3
       return
    end if

    IF (vstate % ISTATE .EQ. 2) GO TO 200
    GO TO 20
10  continue
    vstate % INIT = 0
    IF (vstate % TOUT .EQ. vstate % T) RETURN

    ! -----------------------------------------------------------------------
    !  Block B.
    !  The next code block is executed for the initial call (vstate % ISTATE = 1),
    !  or for a continuation call with parameter changes (vstate % ISTATE = 3).
    !  It contains checking of all input and various initializations.
    !
    !  First check legality of the non-optional input ITOL, IOPT,
    !  MF, ML, and MU.
    ! -----------------------------------------------------------------------
20  continue
    if (VODE_ITOL .LT. 1 .OR. VODE_ITOL .GT. 4) then
#ifndef AMREX_USE_CUDA
       MSG = 'DVODE--  ITOL (=I1) illegal   '
       CALL XERRWD (MSG, 30, 6, 1, 1, VODE_ITOL, 0, 0, ZERO, ZERO)
#endif
       vstate % ISTATE = -3
       return
    end if

    if (IOPT .LT. 0 .OR. IOPT .GT. 1) then
#ifndef AMREX_USE_CUDA
       MSG = 'DVODE--  IOPT (=I1) illegal   '
       CALL XERRWD (MSG, 30, 7, 1, 1, IOPT, 0, 0, ZERO, ZERO)
#endif
       vstate % ISTATE = -3
       return
    end if

    vstate % JSV = SIGN(1,MF)
    MFA = ABS(MF)
    vstate % METH = MFA/10
    vstate % MITER = MFA - 10*vstate % METH

    if (vstate % METH .LT. 1 .OR. vstate % METH .GT. 2 .or. &
        vstate % MITER .LT. 0 .OR. vstate % MITER .GT. 5) then
#ifndef AMREX_USE_CUDA
       MSG = 'DVODE--  MF (=I1) illegal     '
       CALL XERRWD (MSG, 30, 8, 1, 1, MF, 0, 0, ZERO, ZERO)
#endif
       vstate % ISTATE = -3
       return
    end if

    IF (vstate % MITER .LE. 3) GO TO 30
    ML = IWORK(1)
    MU = IWORK(2)

    if (ML .LT. 0 .OR. ML .GE. VODE_NEQS) then
#ifndef AMREX_USE_CUDA
       MSG = 'DVODE--  ML (=I1) illegal:  .lt.0 or .ge.NEQ (=I2)'
       CALL XERRWD (MSG, 50, 9, 1, 2, ML, VODE_NEQS, 0, ZERO, ZERO)
#endif
       vstate % ISTATE = -3
       return
    end if

    if (MU .LT. 0 .OR. MU .GE. VODE_NEQS) then
#ifndef AMREX_USE_CUDA
       MSG = 'DVODE--  MU (=I1) illegal:  .lt.0 or .ge.NEQ (=I2)'
       CALL XERRWD (MSG, 50, 10, 1, 2, MU, VODE_NEQS, 0, ZERO, ZERO)
#endif
       vstate % ISTATE = -3
       return
    end if

30  CONTINUE

    ! Next process and check the optional input. ---------------------------

    IF (IOPT .EQ. 1) GO TO 40
    vstate % MXSTEP = MXSTP0
    vstate % MXHNIL = MXHNL0
    IF (vstate % ISTATE .EQ. 1) H0 = ZERO
    vstate % HMXI = ZERO
    vstate % HMIN = ZERO
    GO TO 60
40  continue
    vstate % MXSTEP = IWORK(6)

    if (vstate % MXSTEP .LT. 0) then
#ifndef AMREX_USE_CUDA
       MSG = 'DVODE--  MXSTEP (=I1) .lt. 0  '
       CALL XERRWD (MSG, 30, 12, 1, 1, vstate % MXSTEP, 0, 0, ZERO, ZERO)
#endif
       vstate % ISTATE = -3
       return
    end if

    IF (vstate % MXSTEP .EQ. 0) vstate % MXSTEP = MXSTP0
    vstate % MXHNIL = IWORK(7)

    if (vstate % MXHNIL .LT. 0) then
#ifndef AMREX_USE_CUDA
       MSG = 'DVODE--  MXHNIL (=I1) .lt. 0  '
       CALL XERRWD (MSG, 30, 13, 1, 1, vstate % MXHNIL, 0, 0, ZERO, ZERO)
#endif
       vstate % ISTATE = -3
       return
    end if

    !      EDIT 07/16/2016 -- see comments above about MXHNIL
    !      IF (MXHNIL .EQ. 0) MXHNIL = MXHNL0

    IF (vstate % ISTATE .NE. 1) GO TO 50
    H0 = RWORK % CONDOPT(2)

    if ((vstate % TOUT - vstate % T)*H0 .LT. ZERO) then
#ifndef AMREX_USE_CUDA
       MSG = 'DVODE--  TOUT (=R1) behind T (=R2)      '
       CALL XERRWD (MSG, 40, 14, 1, 0, 0, 0, 2, vstate % TOUT, vstate % T)
       MSG = '      integration direction is given by H0 (=R1)  '
       CALL XERRWD (MSG, 50, 14, 1, 0, 0, 0, 1, H0, ZERO)
#endif
       vstate % ISTATE = -3
       return
    end if

50  continue
    HMAX = RWORK % CONDOPT(3)

    if (HMAX .LT. ZERO) then
#ifndef AMREX_USE_CUDA
       MSG = 'DVODE--  HMAX (=R1) .lt. 0.0_rt  '
       CALL XERRWD (MSG, 30, 15, 1, 0, 0, 0, 1, HMAX, ZERO)
#endif
       vstate % ISTATE = -3
       return
    end if

    vstate % HMXI = ZERO
    IF (HMAX .GT. ZERO) vstate % HMXI = ONE/HMAX
    vstate % HMIN = RWORK % CONDOPT(4)

    if (vstate % HMIN .LT. ZERO) then
#ifndef AMREX_USE_CUDA
       MSG = 'DVODE--  HMIN (=R1) .lt. 0.0_rt  '
       CALL XERRWD (MSG, 30, 16, 1, 0, 0, 0, 1, vstate % HMIN, ZERO)
#endif
       vstate % ISTATE = -3
       return
    end if


    ! -----------------------------------------------------------------------
    !  Arrays stored in RWORK are denoted  CONDOPT, YH, WM, EWT, SAVF, ACOR.
    !  Within WM, LOCJS is the location of the saved Jacobian (JSV .gt. 0).
    ! -----------------------------------------------------------------------

60  continue
    JCO = MAX(0,vstate % JSV)
    IF (vstate % MITER .EQ. 1 .OR. vstate % MITER .EQ. 2) THEN
       vstate % LOCJS = VODE_NEQS*VODE_NEQS + 3
    ENDIF
    IF (vstate % MITER .EQ. 4 .OR. vstate % MITER .EQ. 5) THEN
       MBAND = ML + MU + 1
       LENP = (MBAND + ML)*VODE_NEQS
       LENJ = MBAND*VODE_NEQS
       vstate % LOCJS = LENP + 3
    ENDIF

    ! Check RTOL and ATOL for legality. ------------------------------------
    RTOLI = vstate % RTOL(1)
    ATOLI = vstate % ATOL(1)
    do I = 1,VODE_NEQS

       IF (VODE_ITOL .GE. 3) RTOLI = vstate % RTOL(I)
       IF (VODE_ITOL .EQ. 2 .OR. VODE_ITOL .EQ. 4) ATOLI = vstate % ATOL(I)

       if (RTOLI .LT. ZERO) then
#ifndef AMREX_USE_CUDA
          MSG = 'DVODE--  RTOL(I1) is R1 .lt. 0.0_rt        '
          CALL XERRWD (MSG, 40, 19, 1, 1, I, 0, 1, RTOLI, ZERO)
#endif
          vstate % ISTATE = -3
          return
       end if

       if (ATOLI .LT. ZERO) then
#ifndef AMREX_USE_CUDA
          MSG = 'DVODE--  ATOL(I1) is R1 .lt. 0.0_rt        '
          CALL XERRWD (MSG, 40, 20, 1, 1, I, 0, 1, ATOLI, ZERO)
#endif
          vstate % ISTATE = -3
          return
       end if

    end do
    IF (vstate % ISTATE .EQ. 1) GO TO 100
    ! If vstate % ISTATE = 3, set flag to signal parameter changes to DVSTEP. -------
    vstate % JSTART = -1
    IF (vstate % NQ .LE. VODE_MAXORD) GO TO 90
    ! MAXORD was reduced below NQ.  Copy YH(*,MAXORD+2) into SAVF. ---------
    rwork % savf(1:VODE_NEQS) = rwork % wm(1:VODE_NEQS)

    ! Reload WM(1) = RWORK % wm(1), since LWM may have changed. ---------------
90  continue
    IF (vstate % MITER .GT. 0) rwork % wm(1) = SQRT(vstate % UROUND)
    GO TO 200

    ! -----------------------------------------------------------------------
    !  Block C.
    !  The next block is for the initial call only (vstate % ISTATE = 1).
    !  It contains all remaining initializations, the initial call to F,
    !  and the calculation of the initial step size.
    !  The error weights in EWT are inverted after being loaded.
    ! -----------------------------------------------------------------------

100 continue
    vstate % UROUND = epsilon(1.0_rt)
    vstate % TN = vstate % T
    IF (ITASK .NE. 4 .AND. ITASK .NE. 5) GO TO 110
    TCRIT = RWORK % condopt(1)
    if ((TCRIT - vstate % TOUT)*(vstate % TOUT - vstate % T) .LT. ZERO) then
#ifndef AMREX_USE_CUDA
       MSG='DVODE--  ITASK = 4 or 5 and TCRIT (=R1) behind TOUT (=R2)   '
       CALL XERRWD (MSG, 60, 25, 1, 0, 0, 0, 2, TCRIT, vstate % TOUT)
#endif
       vstate % ISTATE = -3
       return
    end if

    if (H0 .NE. ZERO .AND. (vstate % T + H0 - TCRIT)*H0 .GT. ZERO) then
       H0 = TCRIT - vstate % T
    end if
110 continue
    vstate % JSTART = 0
    IF (vstate % MITER .GT. 0) RWORK % wm(1) = SQRT(vstate % UROUND)
    vstate % CCMXJ = PT2
    vstate % MSBJ = 50
    vstate % NHNIL = 0
    vstate % NST = 0
    vstate % NJE = 0
    vstate % NNI = 0
    vstate % NCFN = 0
    vstate % NETF = 0
    vstate % NLU = 0
    vstate % NSLJ = 0
    NSLAST = 0
    vstate % HU = ZERO
    vstate % NQU = 0

    ! Initial call to F.  -------------------------

    CALL f_rhs (vstate % T, vstate % Y, rwork % yh(:,2), vstate % RPAR)
    vstate % NFE = 1
    ! Load the initial value array in YH. ---------------------------------
    rwork % YH(1:VODE_NEQS,1) = vstate % Y(1:VODE_NEQS)

    ! Load and invert the EWT array.  (H is temporarily set to 1.0.) -------
    vstate % NQ = 1
    vstate % H = ONE
    CALL DEWSET (vstate, rwork)
    do I = 1,VODE_NEQS
       if (rwork % ewt(I) .LE. ZERO) then
#ifndef AMREX_USE_CUDA
          EWTI = rwork % ewt(I)
          MSG = 'DVODE--  EWT(I1) is R1 .le. 0.0_rt         '
          CALL XERRWD (MSG, 40, 21, 1, 1, I, 0, 1, EWTI, ZERO)
#endif
          vstate % ISTATE = -3
          return
       end if

       rwork % ewt(I) = ONE/rwork % ewt(I)
    end do
    IF (H0 .NE. ZERO) GO TO 180

    ! Call DVHIN to set initial step size H0 to be attempted. --------------
    CALL DVHIN (vstate, rwork, H0, NITER, IER)
    vstate % NFE = vstate % NFE + NITER

    if (IER .NE. 0) then
#ifndef AMREX_USE_CUDA
       MSG='DVODE--  TOUT (=R1) too close to T(=R2) to start integration'
       CALL XERRWD (MSG, 60, 22, 1, 0, 0, 0, 2, vstate % TOUT, vstate % T)
#endif
       vstate % ISTATE = -3
       return
    end if

    ! Adjust H0 if necessary to meet HMAX bound. ---------------------------
180 continue
    RH = ABS(H0)*vstate % HMXI
    IF (RH .GT. ONE) H0 = H0/RH
    ! Load H with H0 and scale YH(*,2) by H0. ------------------------------
    vstate % H = H0
    rwork % YH(:,2) = rwork % YH(:,2) * H0

    GO TO 270

    ! -----------------------------------------------------------------------
    !  Block D.
    !  The next code block is for continuation calls only (vstate % ISTATE = 2 or 3)
    !  and is to check stop conditions before taking a step.
    ! -----------------------------------------------------------------------

200 continue
    NSLAST = vstate % NST
    vstate % KUTH = 0

    select case (ITASK)
    case (1)
       go to 210
    case (2)
       go to 250
    case (3)
       go to 220
    case (4)
       go to 230
    case (5)
       go to 240
    end select

210 IF ((vstate % TN - vstate % TOUT) * vstate % H .LT. ZERO) GO TO 250
    CALL DVINDY (vstate, rwork, IFLAG)

    if (IFLAG .NE. 0) then
#ifndef AMREX_USE_CUDA
       MSG='DVODE--  Trouble from DVINDY.  ITASK = I1, TOUT = R1.       '
       CALL XERRWD (MSG, 60, 27, 1, 1, ITASK, 0, 1, vstate % TOUT, ZERO)
#endif
       vstate % ISTATE = -3
       return
    end if

    vstate % T = vstate % TOUT
    GO TO 420
220 continue
    TP = vstate % TN - vstate % HU * (ONE + HUN * vstate % UROUND)

    if ((TP - vstate % TOUT) * vstate % H .GT. ZERO) then
#ifndef AMREX_USE_CUDA
       MSG='DVODE--  ITASK = I1 and TOUT (=R1) behind TCUR - HU (= R2)  '
       CALL XERRWD (MSG, 60, 23, 1, 1, ITASK, 0, 2, vstate % TOUT, TP)
#endif
       vstate % ISTATE = -3
       return
    end if

    IF ((vstate % TN - vstate % TOUT) * vstate % H .LT. ZERO) GO TO 250
    GO TO 400
230 continue
    TCRIT = RWORK % condopt(1)

    if ((vstate % TN - TCRIT) * vstate % H .GT. ZERO) then
#ifndef AMREX_USE_CUDA
       MSG='DVODE--  ITASK = 4 or 5 and TCRIT (=R1) behind TCUR (=R2)   '
       CALL XERRWD (MSG, 60, 24, 1, 0, 0, 0, 2, TCRIT, vstate % TN)
#endif
       vstate % ISTATE = -3
       return
    end if


    if ((TCRIT - vstate % TOUT) * vstate % H .LT. ZERO) then
#ifndef AMREX_USE_CUDA
       MSG='DVODE--  ITASK = 4 or 5 and TCRIT (=R1) behind TOUT (=R2)   '
       CALL XERRWD (MSG, 60, 25, 1, 0, 0, 0, 2, TCRIT, vstate % TOUT)
#endif
       vstate % ISTATE = -3
       return
    end if

    IF ((vstate % TN - vstate % TOUT) * vstate % H .LT. ZERO) GO TO 245
    CALL DVINDY (vstate, rwork, IFLAG)

    if (IFLAG .NE. 0) then
#ifndef AMREX_USE_CUDA
       MSG='DVODE--  Trouble from DVINDY.  ITASK = I1, TOUT = R1.       '
       CALL XERRWD (MSG, 60, 27, 1, 1, ITASK, 0, 1, vstate % TOUT, ZERO)
#endif
       vstate % ISTATE = -3
       return
    end if

    vstate % T = vstate % TOUT
    GO TO 420
240 continue
    TCRIT = RWORK % condopt(1)

    if ((vstate % TN - TCRIT) * vstate % H .GT. ZERO) then
#ifndef AMREX_USE_CUDA
       MSG='DVODE--  ITASK = 4 or 5 and TCRIT (=R1) behind TCUR (=R2)   '
       CALL XERRWD (MSG, 60, 24, 1, 0, 0, 0, 2, TCRIT, vstate % TN)
#endif
       vstate % ISTATE = -3
       return
    end if

245 continue
    HMX = ABS(vstate % TN) + ABS(vstate % H)
    IHIT = ABS(vstate % TN - TCRIT) .LE. HUN * vstate % UROUND * HMX
    IF (IHIT) GO TO 400
    TNEXT = vstate % TN + vstate % HNEW*(ONE + FOUR * vstate % UROUND)
    IF ((TNEXT - TCRIT) * vstate % H .LE. ZERO) GO TO 250
    vstate % H = (TCRIT - vstate % TN)*(ONE - FOUR * vstate % UROUND)
    vstate % KUTH = 1

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
    IF ((vstate % NST-NSLAST) .GE. vstate % MXSTEP) GO TO 500
    CALL DEWSET (vstate, rwork)
    do I = 1,VODE_NEQS
       IF (rwork % ewt(I) .LE. ZERO) GO TO 510
       rwork % ewt(I) = ONE/rwork % ewt(I)
    end do
270 continue
    TOLSF = vstate % UROUND * DVNORM (rwork % YH(:,1), rwork % EWT)
    IF (TOLSF .LE. ONE) GO TO 280
    TOLSF = TOLSF*TWO

    if (vstate % NST .EQ. 0) then
#ifndef AMREX_USE_CUDA
       MSG = 'DVODE--  At start of problem, too much accuracy   '
       CALL XERRWD (MSG, 50, 26, 1, 0, 0, 0, 0, ZERO, ZERO)
       MSG='      requested for precision of machine:   see TOLSF (=R1) '
       CALL XERRWD (MSG, 60, 26, 1, 0, 0, 0, 1, TOLSF, ZERO)
#endif
       vstate % ISTATE = -3
       return
    end if

    GO TO 520
280 IF ((vstate % TN + vstate % H) .NE. vstate % TN) GO TO 290
    vstate % NHNIL = vstate % NHNIL + 1
    IF (vstate % NHNIL .GT. vstate % MXHNIL) GO TO 290
#ifndef AMREX_USE_CUDA
    MSG = 'DVODE--  Warning: internal T (=R1) and H (=R2) are'
    CALL XERRWD (MSG, 50, 101, 1, 0, 0, 0, 0, ZERO, ZERO)
    MSG='      such that in the machine, T + H = T on the next step  '
    CALL XERRWD (MSG, 60, 101, 1, 0, 0, 0, 0, ZERO, ZERO)
    MSG = '      (H = step size). solver will continue anyway'
    CALL XERRWD (MSG, 50, 101, 1, 0, 0, 0, 2, vstate % TN, vstate % H)
#endif
    IF (vstate % NHNIL .LT. vstate % MXHNIL) GO TO 290
#ifndef AMREX_USE_CUDA
    MSG = 'DVODE--  Above warning has been issued I1 times.  '
    CALL XERRWD (MSG, 50, 102, 1, 0, 0, 0, 0, ZERO, ZERO)
    MSG = '      it will not be issued again for this problem'
    CALL XERRWD (MSG, 50, 102, 1, 1, vstate % MXHNIL, 0, 0, ZERO, ZERO)
#endif
290 CONTINUE

    ! -----------------------------------------------------------------------
    !  CALL DVSTEP (Y, YH, NYH, YH, EWT, SAVF, VSAV, ACOR,
    !               WM, IWM, F, JAC, F, DVNLSD, RPAR, IPAR)
    ! -----------------------------------------------------------------------
    CALL DVSTEP(IWORK, rwork, vstate)
    KGO = 1 - vstate % KFLAG
    ! Branch on KFLAG.  Note: In this version, KFLAG can not be set to -3.
    !  KFLAG .eq. 0,   -1,  -2

    select case (KGO)
    case (1)
       go to 300
    case (2)
       go to 530
    case (3)
       go to 540
    end select

    ! -----------------------------------------------------------------------
    !  Block F.
    !  The following block handles the case of a successful return from the
    !  core integrator (KFLAG = 0).  Test for stop conditions.
    ! -----------------------------------------------------------------------

300 continue
    vstate % INIT = 1
    vstate % KUTH = 0

    select case (ITASK)
    case (1)
       go to 310
    case (2)
       go to 400
    case (3)
       go to 330
    case (4)
       go to 340
    case (5)
       go to 350
    end select

    ! ITASK = 1.  If TOUT has been reached, interpolate. -------------------
310 IF ((vstate % TN - vstate % TOUT) * vstate % H .LT. ZERO) GO TO 250
    CALL DVINDY (vstate, rwork, IFLAG)
    vstate % T = vstate % TOUT
    GO TO 420
    ! ITASK = 3.  Jump to exit if TOUT was reached. ------------------------
330 IF ((vstate % TN - vstate % TOUT) * vstate % H .GE. ZERO) GO TO 400
    GO TO 250
    ! ITASK = 4.  See if TOUT or TCRIT was reached.  Adjust H if necessary.
340 IF ((vstate % TN - vstate % TOUT) * vstate % H .LT. ZERO) GO TO 345
    CALL DVINDY (vstate, rwork, IFLAG)
    vstate % T = vstate % TOUT
    GO TO 420
345 continue
    HMX = ABS(vstate % TN) + ABS(vstate % H)
    IHIT = ABS(vstate % TN - TCRIT) .LE. HUN * vstate % UROUND * HMX
    IF (IHIT) GO TO 400
    TNEXT = vstate % TN + vstate % HNEW*(ONE + FOUR * vstate % UROUND)
    IF ((TNEXT - TCRIT) * vstate % H .LE. ZERO) GO TO 250
    vstate % H = (TCRIT - vstate % TN)*(ONE - FOUR * vstate % UROUND)
    vstate % KUTH = 1
    GO TO 250
    ! ITASK = 5.  See if TCRIT was reached and jump to exit. ---------------
350 continue
    HMX = ABS(vstate % TN) + ABS(vstate % H)
    IHIT = ABS(vstate % TN - TCRIT) .LE. HUN * vstate % UROUND * HMX

    ! -----------------------------------------------------------------------
    !  Block G.
    !  The following block handles all successful returns from DVODE.
    !  If ITASK .ne. 1, Y is loaded from YH and T is set accordingly.
    !  vstate % ISTATE is set to 2, and the optional output is loaded into the work
    !  arrays before returning.
    ! -----------------------------------------------------------------------

400 CONTINUE
    vstate % Y(1:VODE_NEQS) = rwork % YH(1:VODE_NEQS,1)

    vstate % T = vstate % TN
    IF (ITASK .NE. 4 .AND. ITASK .NE. 5) GO TO 420
    IF (IHIT) vstate % T = TCRIT
420 continue
    vstate % ISTATE = 2
    IWORK(11) = vstate % NST
    IWORK(12) = vstate % NFE
    IWORK(13) = vstate % NJE
    IWORK(14) = vstate % NQU
    IWORK(15) = vstate % NEWQ
    IWORK(19) = vstate % NLU
    IWORK(20) = vstate % NNI
    IWORK(21) = vstate % NCFN
    IWORK(22) = vstate % NETF

    return

    ! -----------------------------------------------------------------------
    !  Block H.
    !  The following block handles all unsuccessful returns other than
    !  those for illegal input.  First the error message routine is called.
    !  if there was an error test or convergence test failure, IMXER is set.
    !  Then Y is loaded from YH, and T is set to TN.
    !  The optional output is loaded into the work arrays before returning.
    ! -----------------------------------------------------------------------

    ! The maximum number of steps was taken before reaching TOUT. ----------
500 continue
#ifndef AMREX_USE_CUDA
    MSG = 'DVODE--  At current T (=R1), MXSTEP (=I1) steps   '
    CALL XERRWD (MSG, 50, 201, 1, 0, 0, 0, 0, ZERO, ZERO)
    MSG = '      taken on this call before reaching TOUT     '
    CALL XERRWD (MSG, 50, 201, 1, 1, vstate % MXSTEP, 0, 1, vstate % TN, ZERO)
#endif
    vstate % ISTATE = -1
    GO TO 580
    ! EWT(i) .le. 0.0 for some i (not at start of problem). ----------------
510 continue
#ifndef AMREX_USE_CUDA
    EWTI = rwork % ewt(I)
    MSG = 'DVODE--  At T (=R1), EWT(I1) has become R2 .le. 0.'
    CALL XERRWD (MSG, 50, 202, 1, 1, I, 0, 2, vstate % TN, EWTI)
#endif
    vstate % ISTATE = -6
    GO TO 580
    ! Too much accuracy requested for machine precision. -------------------
520 continue
#ifndef AMREX_USE_CUDA
    MSG = 'DVODE--  At T (=R1), too much accuracy requested  '
    CALL XERRWD (MSG, 50, 203, 1, 0, 0, 0, 0, ZERO, ZERO)
    MSG = '      for precision of machine:   see TOLSF (=R2) '
    CALL XERRWD (MSG, 50, 203, 1, 0, 0, 0, 2, vstate % TN, TOLSF)
#endif
    vstate % ISTATE = -2
    GO TO 580
    ! KFLAG = -1.  Error test failed repeatedly or with ABS(H) = HMIN. -----
530 continue
#ifndef AMREX_USE_CUDA
    MSG = 'DVODE--  At T(=R1) and step size H(=R2), the error'
    CALL XERRWD (MSG, 50, 204, 1, 0, 0, 0, 0, ZERO, ZERO)
    MSG = '      test failed repeatedly or with abs(H) = HMIN'
    CALL XERRWD (MSG, 50, 204, 1, 0, 0, 0, 2, vstate % TN, vstate % H)
#endif
    vstate % ISTATE = -4
    GO TO 560
    ! KFLAG = -2.  Convergence failed repeatedly or with ABS(H) = HMIN. ----
540 continue
#ifndef AMREX_USE_CUDA
    MSG = 'DVODE--  At T (=R1) and step size H (=R2), the    '
    CALL XERRWD (MSG, 50, 205, 1, 0, 0, 0, 0, ZERO, ZERO)
    MSG = '      corrector convergence failed repeatedly     '
    CALL XERRWD (MSG, 50, 205, 1, 0, 0, 0, 0, ZERO, ZERO)
    MSG = '      or with abs(H) = HMIN   '
    CALL XERRWD (MSG, 30, 205, 1, 0, 0, 0, 2, vstate % TN, vstate % H)
#endif
    vstate % ISTATE = -5
    ! Compute IMXER if relevant. -------------------------------------------
560 continue
    BIG = ZERO
    IMXER = 1
    do I = 1,VODE_NEQS
       SIZE = ABS(rwork % acor(I) * rwork % ewt(I))
       IF (BIG .GE. SIZE) exit
       BIG = SIZE
       IMXER = I
    end do
    IWORK(16) = IMXER
    ! Set Y array, T, and optional output. --------------------------------
580 CONTINUE
    vstate % Y(1:VODE_NEQS) = rwork % YH(1:VODE_NEQS,1)

    vstate % T = vstate % TN
    IWORK(11) = vstate % NST
    IWORK(12) = vstate % NFE
    IWORK(13) = vstate % NJE
    IWORK(14) = vstate % NQU
    IWORK(15) = vstate % NQ
    IWORK(19) = vstate % NLU
    IWORK(20) = vstate % NNI
    IWORK(21) = vstate % NCFN
    IWORK(22) = vstate % NETF

    return

  end subroutine dvode
      
end module cuvode_module
