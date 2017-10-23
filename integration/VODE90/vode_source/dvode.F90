module dvode_module

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

  real(dp_t), parameter :: ZERO = 0.0D0
  real(dp_t), parameter :: ONE  = 1.0D0
  real(dp_t), parameter :: HALF = 0.5D0
  real(dp_t), parameter :: TWO  = 2.0D0
  real(dp_t), parameter :: FOUR = 4.0D0
  real(dp_t), parameter :: SIX  = 6.0D0
  real(dp_t), parameter :: HUN  = 100.0D0
  real(dp_t), parameter :: THOU = 1000.0D0

  public :: dvode
  
contains

#ifdef CUDA
  attributes(device) &
#endif  
  function dumach() result(dum)

    !$acc routine seq
    
    ! ***BEGIN PROLOGUE  DUMACH
    ! ***PURPOSE  Compute the unit roundoff of the machine.
    ! ***CATEGORY  R1
    ! ***TYPE      DOUBLE PRECISION (RUMACH-S, DUMACH-D)
    ! ***KEYWORDS  MACHINE CONSTANTS
    ! ***AUTHOR  Hindmarsh, Alan C., (LLNL)
    ! ***DESCRIPTION
    !  *Usage:
    !         DOUBLE PRECISION  A, DUMACH
    !         A = DUMACH()
    ! 
    !  *Function Return Values:
    !      A : the unit roundoff of the machine.
    ! 
    !  *Description:
    !      The unit roundoff is defined as the smallest positive machine
    !      number u such that  1.0 + u .ne. 1.0.  This is computed by DUMACH
    !      in a machine-independent manner.
    ! 
    ! ***REFERENCES  (NONE)
    ! ***ROUTINES CALLED  (NONE)
    ! ***REVISION HISTORY  (YYMMDD)
    !    930216  DATE WRITTEN
    !    930818  Added SLATEC-format prologue.  (FNF)
    ! ***END PROLOGUE  DUMACH
    ! 
    ! *Internal Notes:
    ! -----------------------------------------------------------------------
    !  Subroutines/functions called by DUMACH.. None
    ! -----------------------------------------------------------------------
    ! **End
    !

    implicit none
  
    real(dp_t) :: U, dum
    dum = EPSILON(U)
  end function dumach

#ifdef CUDA
  attributes(device) &
#endif  
  subroutine dewset(N, ITOL, RTOL, ATOL, YCUR, EWT)

    !$acc routine seq
    
    ! ***BEGIN PROLOGUE  DEWSET
    ! ***SUBSIDIARY
    ! ***PURPOSE  Set error weight vector.
    ! ***TYPE      DOUBLE PRECISION (SEWSET-S, DEWSET-D)
    ! ***AUTHOR  Hindmarsh, Alan C., (LLNL)
    ! ***DESCRIPTION
    ! 
    !   This subroutine sets the error weight vector EWT according to
    !       EWT(i) = RTOL(i)*ABS(YCUR(i)) + ATOL(i),  i = 1,...,N,
    !   with the subscript on RTOL and/or ATOL possibly replaced by 1 above,
    !   depending on the value of ITOL.
    ! 
    ! ***SEE ALSO  DLSODE
    ! ***ROUTINES CALLED  (NONE)
    ! ***REVISION HISTORY  (YYMMDD)
    !    791129  DATE WRITTEN
    !    890501  Modified prologue to SLATEC/LDOC format.  (FNF)
    !    890503  Minor cosmetic changes.  (FNF)
    !    930809  Renamed to allow single/double precision versions. (ACH)
    ! ***END PROLOGUE  DEWSET
    ! **End

    implicit none
  
    integer    :: I, N, ITOL
    real(dp_t) :: RTOL(:), ATOL(:)
    real(dp_t) :: YCUR(:), EWT(:)

    GO TO (10, 20, 30, 40), ITOL
10  CONTINUE
    do I = 1,N
       EWT(I) = RTOL(1)*ABS(YCUR(I)) + ATOL(1)
    end do
    RETURN
20  CONTINUE
    do I = 1,N
       EWT(I) = RTOL(1)*ABS(YCUR(I)) + ATOL(I)
    end do
    RETURN
30  CONTINUE
    do I = 1,N
       EWT(I) = RTOL(I)*ABS(YCUR(I)) + ATOL(1)
    end do
    RETURN
40  CONTINUE

    do I = 1,N
       EWT(I) = RTOL(I)*ABS(YCUR(I)) + ATOL(I)
    end do
    RETURN
  end subroutine dewset

#ifdef CUDA
  attributes(device) &
#endif  
  function dvnorm(N, V, W) result(dvn)

    !$acc routine seq
    
    ! ***BEGIN PROLOGUE  DVNORM
    ! ***SUBSIDIARY
    ! ***PURPOSE  Weighted root-mean-square vector norm.
    ! ***TYPE      DOUBLE PRECISION (SVNORM-S, DVNORM-D)
    ! ***AUTHOR  Hindmarsh, Alan C., (LLNL)
    ! ***DESCRIPTION
    ! 
    !   This function routine computes the weighted root-mean-square norm
    !   of the vector of length N contained in the array V, with weights
    !   contained in the array W of length N:
    !     DVNORM = SQRT( (1/N) * SUM( V(i)*W(i) )**2 )
    ! 
    ! ***SEE ALSO  DLSODE
    ! ***ROUTINES CALLED  (NONE)
    ! ***REVISION HISTORY  (YYMMDD)
    !    791129  DATE WRITTEN
    !    890501  Modified prologue to SLATEC/LDOC format.  (FNF)
    !    890503  Minor cosmetic changes.  (FNF)
    !    930809  Renamed to allow single/double precision versions. (ACH)
    ! ***END PROLOGUE  DVNORM
    ! **End

    implicit none

    integer    :: N, I
    real(dp_t) :: V(:), W(:)
    real(dp_t) :: SUM, dvn, dscratch

    SUM = 0.0D0
    do I = 1,N
       SUM = SUM + (V(I)*W(I))**2
    end do
    dscratch = SUM/N
    dvn = SQRT(dscratch)
    RETURN
  end function dvnorm
  
#ifdef CUDA
  attributes(device) &
#endif
  subroutine dvhin(N, T0, YH, RPAR, IPAR, TOUT, UROUND, &
       EWT, ITOL, ATOL, Y, TEMP, H0, NITER, IER)

    !$acc routine seq
    
    ! -----------------------------------------------------------------------
    !  Call sequence input -- N, T0, Y0, YDOT, F, RPAR, IPAR, TOUT, UROUND,
    !                         EWT, ITOL, ATOL, Y, TEMP
    !  Call sequence output -- H0, NITER, IER
    !  COMMON block variables accessed -- None
    ! 
    !  Subroutines called by DVHIN:  F
    !  Function routines called by DVHI: DVNORM
    ! -----------------------------------------------------------------------
    !  This routine computes the step size, H0, to be attempted on the
    !  first step, when the user has not supplied a value for this.
    ! 
    !  First we check that TOUT - T0 differs significantly from zero.  Then
    !  an iteration is done to approximate the initial second derivative
    !  and this is used to define h from w.r.m.s.norm(h**2 * yddot / 2) = 1.
    !  A bias factor of 1/2 is applied to the resulting h.
    !  The sign of H0 is inferred from the initial values of TOUT and T0.
    ! 
    !  Communication with DVHIN is done with the following variables:
    ! 
    !  N      = Size of ODE system, input.
    !  T0     = Initial value of independent variable, input.
    !
    !  Y0     = Vector of initial conditions, input.
    !  YDOT   = Vector of initial first derivatives, input.
    !  [NOTE: Y0 = YH(:,1) and YDOT = YH(:,2) in this subroutine now]
    !
    !  F      = Name of subroutine for right-hand side f(t,y), input.
    !  RPAR, IPAR = Dummy names for user's real and integer work arrays.
    !  TOUT   = First output value of independent variable
    !  UROUND = Machine unit roundoff
    !  EWT, ITOL, ATOL = Error weights and tolerance parameters
    !                    as described in the driver routine, input.
    !  Y, TEMP = Work arrays of length N.
    !  H0     = Step size to be attempted, output.
    !  NITER  = Number of iterations (and of f evaluations) to compute H0,
    !           output.
    !  IER    = The error flag, returned with the value
    !           IER = 0  if no trouble occurred, or
    !           IER = -1 if TOUT and T0 are considered too close to proceed.
    ! -----------------------------------------------------------------------

    use vode_rhs_module, only: f_rhs
  
    implicit none
  
    real(dp_t) :: T0, RPAR(:), TOUT, UROUND
    real(dp_t) :: ATOL(:), Y(:), H0
    integer    :: N, IPAR(:), ITOL, NITER, IER
    real(dp_t) :: YH(:,:)
    real(dp_t) :: TEMP(:), EWT(:)

    real(dp_t) :: AFI, ATOLI, DELYI, H, HG, HLB, HNEW, HRAT
    real(dp_t) :: HUB, T1, TDIST, TROUND, YDDNRM, dscratch
    integer    :: I, ITER

    real(dp_t), parameter :: PT1 = 0.1D0

    NITER = 0
    TDIST = ABS(TOUT - T0)
    TROUND = UROUND*MAX(ABS(T0),ABS(TOUT))
    IF (TDIST .LT. TWO*TROUND) GO TO 100

    ! Set a lower bound on h based on the roundoff level in T0 and TOUT. ---
    HLB = HUN*TROUND
    ! Set an upper bound on h based on TOUT-T0 and the initial Y and YDOT. -
    HUB = PT1*TDIST
    ATOLI = ATOL(1)
    
    do I = 1, N
       IF (ITOL .EQ. 2 .OR. ITOL .EQ. 4) ATOLI = ATOL(I)
       DELYI = PT1*ABS(YH(I,1)) + ATOLI
       AFI = ABS(YH(I,2))
       IF (AFI*HUB .GT. DELYI) HUB = DELYI/AFI
    end do

    ! Set initial guess for h as geometric mean of upper and lower bounds. -
    ITER = 0
    dscratch = HLB*HUB
    HG = SQRT(dscratch)
    ! If the bounds have crossed, exit with the mean value. ----------------
    IF (HUB .LT. HLB) THEN
       H0 = HG
       GO TO 90
    ENDIF

    ! Looping point for iteration. -----------------------------------------
50  CONTINUE
    ! Estimate the second derivative as a difference quotient in f. --------
    dscratch = TOUT - T0
    H = SIGN (HG, dscratch)
    T1 = T0 + H
    do I = 1, N
       Y(I) = YH(I,1) + H*YH(I,2)
    end do
    CALL f_rhs(N, T1, Y, TEMP, RPAR, IPAR)
    do I = 1, N
       TEMP(I) = (TEMP(I) - YH(I,2))/H
    end do
    YDDNRM = DVNORM (N, TEMP, EWT)
    ! Get the corresponding new value of h. --------------------------------
    IF (YDDNRM*HUB*HUB .GT. TWO) THEN
       dscratch = TWO/YDDNRM
       HNEW = SQRT(dscratch)
    ELSE
       dscratch = HG*HUB
       HNEW = SQRT(dscratch)
    ENDIF
    ITER = ITER + 1
    ! -----------------------------------------------------------------------
    !  Test the stopping conditions.
    !  Stop if the new and previous h values differ by a factor of .lt. 2.
    !  Stop if four iterations have been done.  Also, stop with previous h
    !  if HNEW/HG .gt. 2 after first iteration, as this probably means that
    !  the second derivative value is bad because of cancellation error.
    ! -----------------------------------------------------------------------
    IF (ITER .GE. 4) GO TO 80
    HRAT = HNEW/HG
    IF ( (HRAT .GT. HALF) .AND. (HRAT .LT. TWO) ) GO TO 80
    IF ( (ITER .GE. 2) .AND. (HNEW .GT. TWO*HG) ) THEN
       HNEW = HG
       GO TO 80
    ENDIF
    HG = HNEW
    GO TO 50

    ! Iteration done.  Apply bounds, bias factor, and sign.  Then exit. ----
80  continue
    H0 = HNEW*HALF
    IF (H0 .LT. HLB) H0 = HLB
    IF (H0 .GT. HUB) H0 = HUB
90  continue
    dscratch = TOUT - T0
    H0 = SIGN(H0, dscratch)
    NITER = ITER
    IER = 0
    RETURN
    ! Error return for TOUT - T0 too small. --------------------------------
100 continue
    IER = -1
    RETURN
  end subroutine dvhin
  
#ifdef CUDA
  attributes(device) &
#endif
  subroutine dvindy(YH, T, K, DKY, IFLAG, vstate)

    !$acc routine seq
    
    ! -----------------------------------------------------------------------
    !  Call sequence input -- T, K, YH, LDYH
    !  Call sequence output -- DKY, IFLAG
    !  COMMON block variables accessed:
    !      /DVOD01/ --  H, TN, UROUND, L, N, NQ
    !      /DVOD02/ --  HU
    ! 
    !  Subroutines called by DVINDY: DSCAL, XERRWD
    !  Function routines called by DVINDY: None
    ! -----------------------------------------------------------------------
    !  DVINDY computes interpolated values of the K-th derivative of the
    !  dependent variable vector y, and stores it in DKY.  This routine
    !  is called within the package with K = 0 and T = TOUT, but may
    !  also be called by the user for any K up to the current order.
    !  (See detailed instructions in the usage documentation.)
    ! -----------------------------------------------------------------------
    !  The computed values in DKY are gotten by interpolation using the
    !  Nordsieck history array YH.  This array corresponds uniquely to a
    !  vector-valued polynomial of degree NQCUR or less, and DKY is set
    !  to the K-th derivative of this polynomial at T.
    !  The formula for DKY is:
    !               q
    !   DKY(i)  =  sum  c(j,K) * (T - TN)**(j-K) * H**(-j) * YH(i,j+1)
    !              j=K
    !  where  c(j,K) = j*(j-1)*...*(j-K+1), q = NQCUR, TN = TCUR, H = HCUR.
    !  The quantities  NQ = NQCUR, L = NQ+1, N, TN, and H are
    !  communicated by COMMON.  The above sum is done in reverse order.
    !  IFLAG is returned negative if either K or T is out of bounds.
    ! 
    !  Discussion above and comments in driver explain all variables.
    ! -----------------------------------------------------------------------
    !
  
    implicit none
  
    type(dvode_t) :: vstate
    real(dp_t) :: T
    real(dp_t) :: YH(:,:)
    real(dp_t) :: DKY(:)
    integer    :: K, IFLAG

    real(dp_t) :: C, R, S, TFUZZ, TN1, TP
    integer    :: I, IC, J, JB, JB2, JJ, JJ1, JP1
#ifndef CUDA
    character (len=80) :: MSG
#endif

    IFLAG = 0
    IF (K .LT. 0 .OR. K .GT. vstate % NQ) GO TO 80
    TFUZZ = HUN * vstate % UROUND * (vstate % TN + vstate % HU)
    TP = vstate % TN - vstate % HU - TFUZZ
    TN1 = vstate % TN + TFUZZ
    IF ((T-TP)*(T-TN1) .GT. ZERO) GO TO 90

    S = (T - vstate % TN)/vstate % H
    IC = 1
    IF (K .EQ. 0) GO TO 15
    JJ1 = vstate % L - K
    do JJ = JJ1, vstate % NQ
       IC = IC*JJ
    end do
15  continue
    C = REAL(IC)
    do I = 1, vstate % N
       DKY(I) = C * YH(I,vstate % L)
    end do
    IF (K .EQ. vstate % NQ) GO TO 55
    JB2 = vstate % NQ - K
    do JB = 1, JB2
       J = vstate % NQ - JB
       JP1 = J + 1
       IC = 1
       IF (K .EQ. 0) GO TO 35
       JJ1 = JP1 - K
       do JJ = JJ1, J
          IC = IC*JJ
       end do
35     continue
       C = REAL(IC)
       do I = 1, vstate % N
          DKY(I) = C * YH(I,JP1) + S*DKY(I)
       end do
    end do
    IF (K .EQ. 0) RETURN
55  continue
    R = vstate % H**(-K)

    CALL DSCAL (vstate % N, R, DKY, 1)
    RETURN

80  continue
#ifndef CUDA    
    MSG = 'DVINDY-- K (=I1) illegal      '
    CALL XERRWD (MSG, 30, 51, 1, 1, K, 0, 0, ZERO, ZERO)
#endif    
    IFLAG = -1
    RETURN
90  continue
#ifndef CUDA    
    MSG = 'DVINDY-- T (=R1) illegal      '
    CALL XERRWD (MSG, 30, 52, 1, 0, 0, 0, 1, T, ZERO)
    MSG='      T not in interval TCUR - HU (= R1) to TCUR (=R2)      '
    CALL XERRWD (MSG, 60, 52, 1, 0, 0, 0, 2, TP, vstate % TN)
#endif    
     IFLAG = -2
    RETURN
  end subroutine dvindy
  
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

    ! Set LIW and LENWM in vstate to allow for explicit-shape array passing.
    vstate % LIW = LIW
    vstate % LENWM = LENWM

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
    IF (ISTATE .EQ. 1) GO TO 10
    IF (vstate % INIT .NE. 1) GO TO 603
    IF (ISTATE .EQ. 2) GO TO 200
    GO TO 20
10  continue
    vstate % INIT = 0
    IF (TOUT .EQ. T) RETURN
    
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
    IF (NEQ .GT. vstate % N) GO TO 605
25  continue
    vstate % N = NEQ
    IF (ITOL .LT. 1 .OR. ITOL .GT. 4) GO TO 606
    IF (IOPT .LT. 0 .OR. IOPT .GT. 1) GO TO 607
    vstate % JSV = SIGN(1,MF)
    MFA = ABS(MF)
    vstate % METH = MFA/10
    vstate % MITER = MFA - 10*vstate % METH
    IF (vstate % METH .LT. 1 .OR. vstate % METH .GT. 2) GO TO 608
    IF (vstate % MITER .LT. 0 .OR. vstate % MITER .GT. 5) GO TO 608
    IF (vstate % MITER .LE. 3) GO TO 30
    ML = IWORK(1)
    MU = IWORK(2)
    IF (ML .LT. 0 .OR. ML .GE. vstate % N) GO TO 609
    IF (MU .LT. 0 .OR. MU .GE. vstate % N) GO TO 610
30  CONTINUE

    ! Next process and check the optional input. ---------------------------
      
    IF (IOPT .EQ. 1) GO TO 40
    vstate % MAXORD = MORD(vstate % METH)
    vstate % MXSTEP = MXSTP0
    vstate % MXHNIL = MXHNL0
    IF (ISTATE .EQ. 1) H0 = ZERO
    vstate % HMXI = ZERO
    vstate % HMIN = ZERO
    GO TO 60
40  continue
    vstate % MAXORD = IWORK(5)
    IF (vstate % MAXORD .LT. 0) GO TO 611
    IF (vstate % MAXORD .EQ. 0) vstate % MAXORD = 100
    vstate % MAXORD = MIN(vstate % MAXORD,MORD(vstate % METH))
    vstate % MXSTEP = IWORK(6)
    IF (vstate % MXSTEP .LT. 0) GO TO 612
    IF (vstate % MXSTEP .EQ. 0) vstate % MXSTEP = MXSTP0
    vstate % MXHNIL = IWORK(7)
    IF (vstate % MXHNIL .LT. 0) GO TO 613
    !      EDIT 07/16/2016 -- see comments above about MXHNIL
    !      IF (MXHNIL .EQ. 0) MXHNIL = MXHNL0
    IF (ISTATE .NE. 1) GO TO 50
    H0 = RWORK % CONDOPT(5)
    IF ((TOUT - T)*H0 .LT. ZERO) GO TO 614
50  continue
    HMAX = RWORK % CONDOPT(6)
    IF (HMAX .LT. ZERO) GO TO 615
    vstate % HMXI = ZERO
    IF (HMAX .GT. ZERO) vstate % HMXI = ONE/HMAX
    vstate % HMIN = RWORK % CONDOPT(7)
    IF (vstate % HMIN .LT. ZERO) GO TO 616
    
    ! -----------------------------------------------------------------------
    !  Arrays stored in RWORK are denoted  CONDOPT, YH, WM, EWT, SAVF, ACOR.
    !  Within WM, LOCJS is the location of the saved Jacobian (JSV .gt. 0).
    ! -----------------------------------------------------------------------
    
60  continue
    vstate % LYH = 21
    IF (ISTATE .EQ. 1) vstate % NYH = vstate % N
    vstate % LWM = vstate % LYH + (vstate % MAXORD + 1)*vstate % NYH
    JCO = MAX(0,vstate % JSV)
    IF (vstate % MITER .EQ. 0) LENWM = 0
    IF (vstate % MITER .EQ. 1 .OR. vstate % MITER .EQ. 2) THEN
       LENWM = 2 + (1 + JCO)*vstate % N*vstate % N
       vstate % LOCJS = vstate % N*vstate % N + 3
    ENDIF
    IF (vstate % MITER .EQ. 3) LENWM = 2 + vstate % N
    IF (vstate % MITER .EQ. 4 .OR. vstate % MITER .EQ. 5) THEN
       MBAND = ML + MU + 1
       LENP = (MBAND + ML)*vstate % N
       LENJ = MBAND*vstate % N
       LENWM = 2 + LENP + JCO*LENJ
       vstate % LOCJS = LENP + 3
    ENDIF
    vstate % LMAX = vstate % MAXORD + 1
    vstate % LEWT = vstate % LWM + LENWM
    vstate % LSAVF = vstate % LEWT + vstate % N
    vstate % LACOR = vstate % LSAVF + vstate % N
    vstate % NEQ = NEQ
    LENRW = vstate % LACOR + vstate % N - 1
    IWORK(17) = LENRW
    vstate % LIWM = 1
    LENIW = 30 + vstate % N
    IF (vstate % MITER .EQ. 0 .OR. vstate % MITER .EQ. 3) LENIW = 30
    IWORK(18) = LENIW
    IF (LENRW .GT. LRW) GO TO 617
    IF (LENIW .GT. LIW) GO TO 618
    ! Check RTOL and ATOL for legality. ------------------------------------
    RTOLI = RTOL(1)
    ATOLI = ATOL(1)
    do I = 1,vstate % N
       IF (ITOL .GE. 3) RTOLI = RTOL(I)
       IF (ITOL .EQ. 2 .OR. ITOL .EQ. 4) ATOLI = ATOL(I)
       IF (RTOLI .LT. ZERO) GO TO 619
       IF (ATOLI .LT. ZERO) GO TO 620
    end do
    IF (ISTATE .EQ. 1) GO TO 100
    ! If ISTATE = 3, set flag to signal parameter changes to DVSTEP. -------
    vstate % JSTART = -1
    IF (vstate % NQ .LE. vstate % MAXORD) GO TO 90
    ! MAXORD was reduced below NQ.  Copy YH(*,MAXORD+2) into SAVF. ---------
    CALL DCOPY(vstate % N, rwork % wm, 1, rwork % savf, 1)

    ! Reload WM(1) = RWORK % wm(1), since LWM may have changed. ---------------
90  continue
    IF (vstate % MITER .GT. 0) rwork % wm(1) = SQRT(vstate % UROUND)
    GO TO 200

    ! -----------------------------------------------------------------------
    !  Block C.
    !  The next block is for the initial call only (ISTATE = 1).
    !  It contains all remaining initializations, the initial call to F,
    !  and the calculation of the initial step size.
    !  The error weights in EWT are inverted after being loaded.
    ! -----------------------------------------------------------------------
        
100 continue
    vstate % UROUND = DUMACH()
    vstate % TN = T
    IF (ITASK .NE. 4 .AND. ITASK .NE. 5) GO TO 110
    TCRIT = RWORK % condopt(1)
    IF ((TCRIT - TOUT)*(TOUT - T) .LT. ZERO) GO TO 625
    if (H0 .NE. ZERO .AND. (T + H0 - TCRIT)*H0 .GT. ZERO) then
       H0 = TCRIT - T
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

    ! Initial call to F.  (LF0 points to YH(*,2).) -------------------------
    LF0 = vstate % LYH + vstate % NYH

    CALL f_rhs (vstate % N, T, Y, rwork % yh(:,2), RPAR, IPAR)
    vstate % NFE = 1
    ! Load the initial value vector in YH. ---------------------------------
    CALL DCOPY(vstate % N, Y, 1, rwork % YH(:,1), 1)

    ! Load and invert the EWT array.  (H is temporarily set to 1.0.) -------
    vstate % NQ = 1
    vstate % H = ONE
    CALL DEWSET (vstate % N, ITOL, RTOL, ATOL, rwork % YH(:,1), rwork % EWT)
    do I = 1,vstate % N
       IF (rwork % ewt(I) .LE. ZERO) GO TO 621
       rwork % ewt(I) = ONE/rwork % ewt(I)
    end do
    IF (H0 .NE. ZERO) GO TO 180

    ! Call DVHIN to set initial step size H0 to be attempted. --------------
    CALL DVHIN (vstate % N, T, &
         rwork % YH, &
         RPAR, IPAR, TOUT, &
         vstate % UROUND, &
         rwork % EWT, &
         ITOL, ATOL, Y, &
         rwork % ACOR, &
         H0, NITER, IER)
    vstate % NFE = vstate % NFE + NITER
    IF (IER .NE. 0) GO TO 622
    ! Adjust H0 if necessary to meet HMAX bound. ---------------------------
180 continue
    RH = ABS(H0)*vstate % HMXI
    IF (RH .GT. ONE) H0 = H0/RH
    ! Load H with H0 and scale YH(*,2) by H0. ------------------------------
    vstate % H = H0
    CALL DSCAL(vstate % N, H0, rwork % YH(:,2), 1)

    GO TO 270
    
    ! -----------------------------------------------------------------------
    !  Block D.
    !  The next code block is for continuation calls only (ISTATE = 2 or 3)
    !  and is to check stop conditions before taking a step.
    ! -----------------------------------------------------------------------
        
200 continue
    NSLAST = vstate % NST
    vstate % KUTH = 0

    GO TO (210, 250, 220, 230, 240), ITASK
210 IF ((vstate % TN - TOUT) * vstate % H .LT. ZERO) GO TO 250
    CALL DVINDY (rwork % yh, TOUT, 0, Y, IFLAG, vstate)
    IF (IFLAG .NE. 0) GO TO 627
    T = TOUT
    GO TO 420
220 continue
    TP = vstate % TN - vstate % HU * (ONE + HUN * vstate % UROUND)
    IF ((TP - TOUT) * vstate % H .GT. ZERO) GO TO 623
    IF ((vstate % TN - TOUT) * vstate % H .LT. ZERO) GO TO 250
    GO TO 400
230 continue
    TCRIT = RWORK % condopt(1)
    IF ((vstate % TN - TCRIT) * vstate % H .GT. ZERO) GO TO 624
    IF ((TCRIT - TOUT) * vstate % H .LT. ZERO) GO TO 625
    IF ((vstate % TN - TOUT) * vstate % H .LT. ZERO) GO TO 245
    CALL DVINDY (rwork % yh, TOUT, 0, Y, IFLAG, vstate)
    IF (IFLAG .NE. 0) GO TO 627
    T = TOUT
    GO TO 420
240 continue
    TCRIT = RWORK % condopt(1)
    IF ((vstate % TN - TCRIT) * vstate % H .GT. ZERO) GO TO 624
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
    CALL DEWSET (vstate % N, ITOL, RTOL, ATOL,  rwork % YH(:,1), rwork % EWT)
    do I = 1,vstate % N
       IF (rwork % ewt(I) .LE. ZERO) GO TO 510
       rwork % ewt(I) = ONE/rwork % ewt(I)
    end do
270 continue
    TOLSF = vstate % UROUND * DVNORM (vstate % N, rwork % YH(:,1), rwork % EWT)
    IF (TOLSF .LE. ONE) GO TO 280
    TOLSF = TOLSF*TWO
    IF (vstate % NST .EQ. 0) GO TO 626
    GO TO 520
280 IF ((vstate % TN + vstate % H) .NE. vstate % TN) GO TO 290
    vstate % NHNIL = vstate % NHNIL + 1
    IF (vstate % NHNIL .GT. vstate % MXHNIL) GO TO 290
#ifndef CUDA    
    MSG = 'DVODE--  Warning: internal T (=R1) and H (=R2) are'
    CALL XERRWD (MSG, 50, 101, 1, 0, 0, 0, 0, ZERO, ZERO)
    MSG='      such that in the machine, T + H = T on the next step  '
    CALL XERRWD (MSG, 60, 101, 1, 0, 0, 0, 0, ZERO, ZERO)
    MSG = '      (H = step size). solver will continue anyway'
    CALL XERRWD (MSG, 50, 101, 1, 0, 0, 0, 2, vstate % TN, vstate % H)
#endif    
    IF (vstate % NHNIL .LT. vstate % MXHNIL) GO TO 290
#ifndef CUDA    
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
    CALL DVSTEP(Y, IWORK, RPAR, IPAR, rwork, vstate)
    KGO = 1 - vstate % KFLAG
    ! Branch on KFLAG.  Note: In this version, KFLAG can not be set to -3.
    !  KFLAG .eq. 0,   -1,  -2
    GO TO (300, 530, 540), KGO
    
    ! -----------------------------------------------------------------------
    !  Block F.
    !  The following block handles the case of a successful return from the
    !  core integrator (KFLAG = 0).  Test for stop conditions.
    ! -----------------------------------------------------------------------
    
300 continue
    vstate % INIT = 1
    vstate % KUTH = 0
    GO TO (310, 400, 330, 340, 350), ITASK
    ! ITASK = 1.  If TOUT has been reached, interpolate. -------------------
310 IF ((vstate % TN - TOUT) * vstate % H .LT. ZERO) GO TO 250
    CALL DVINDY (rwork % yh, TOUT, 0, Y, IFLAG, vstate)
    T = TOUT
    GO TO 420
    ! ITASK = 3.  Jump to exit if TOUT was reached. ------------------------
330 IF ((vstate % TN - TOUT) * vstate % H .GE. ZERO) GO TO 400
    GO TO 250
    ! ITASK = 4.  See if TOUT or TCRIT was reached.  Adjust H if necessary.
340 IF ((vstate % TN - TOUT) * vstate % H .LT. ZERO) GO TO 345
    CALL DVINDY (rwork % yh, TOUT, 0, Y, IFLAG, vstate)
    T = TOUT
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
    !  ISTATE is set to 2, and the optional output is loaded into the work
    !  arrays before returning.
    ! -----------------------------------------------------------------------
    
400 CONTINUE
    CALL DCOPY(vstate % N, rwork % YH(:,1), 1, Y, 1)

    T = vstate % TN
    IF (ITASK .NE. 4 .AND. ITASK .NE. 5) GO TO 420
    IF (IHIT) T = TCRIT
420 continue
    ISTATE = 2
    RWORK % condopt(11) = vstate % HU
    RWORK % condopt(12) = vstate % HNEW
    RWORK % condopt(13) = vstate % TN
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
#ifndef CUDA    
    MSG = 'DVODE--  At current T (=R1), MXSTEP (=I1) steps   '
    CALL XERRWD (MSG, 50, 201, 1, 0, 0, 0, 0, ZERO, ZERO)
    MSG = '      taken on this call before reaching TOUT     '
    CALL XERRWD (MSG, 50, 201, 1, 1, vstate % MXSTEP, 0, 1, vstate % TN, ZERO)
#endif    
    ISTATE = -1
    GO TO 580
    ! EWT(i) .le. 0.0 for some i (not at start of problem). ----------------
510 continue
#ifndef CUDA        
    EWTI = rwork % ewt(I)
    MSG = 'DVODE--  At T (=R1), EWT(I1) has become R2 .le. 0.'
    CALL XERRWD (MSG, 50, 202, 1, 1, I, 0, 2, vstate % TN, EWTI)
#endif
    ISTATE = -6
    GO TO 580
    ! Too much accuracy requested for machine precision. -------------------
520 continue
#ifndef CUDA
    MSG = 'DVODE--  At T (=R1), too much accuracy requested  '
    CALL XERRWD (MSG, 50, 203, 1, 0, 0, 0, 0, ZERO, ZERO)
    MSG = '      for precision of machine:   see TOLSF (=R2) '
    CALL XERRWD (MSG, 50, 203, 1, 0, 0, 0, 2, vstate % TN, TOLSF)
#endif
    RWORK % condopt(14) = TOLSF
    ISTATE = -2
    GO TO 580
    ! KFLAG = -1.  Error test failed repeatedly or with ABS(H) = HMIN. -----
530 continue
#ifndef CUDA        
    MSG = 'DVODE--  At T(=R1) and step size H(=R2), the error'
    CALL XERRWD (MSG, 50, 204, 1, 0, 0, 0, 0, ZERO, ZERO)
    MSG = '      test failed repeatedly or with abs(H) = HMIN'
    CALL XERRWD (MSG, 50, 204, 1, 0, 0, 0, 2, vstate % TN, vstate % H)
#endif    
    ISTATE = -4
    GO TO 560
    ! KFLAG = -2.  Convergence failed repeatedly or with ABS(H) = HMIN. ----
540 continue
#ifndef CUDA        
    MSG = 'DVODE--  At T (=R1) and step size H (=R2), the    '
    CALL XERRWD (MSG, 50, 205, 1, 0, 0, 0, 0, ZERO, ZERO)
    MSG = '      corrector convergence failed repeatedly     '
    CALL XERRWD (MSG, 50, 205, 1, 0, 0, 0, 0, ZERO, ZERO)
    MSG = '      or with abs(H) = HMIN   '
    CALL XERRWD (MSG, 30, 205, 1, 0, 0, 0, 2, vstate % TN, vstate % H)
#endif    
    ISTATE = -5
    ! Compute IMXER if relevant. -------------------------------------------
560 continue
    BIG = ZERO
    IMXER = 1
    do I = 1,vstate % N
       SIZE = ABS(rwork % acor(I) * rwork % ewt(I))
       IF (BIG .GE. SIZE) exit
       BIG = SIZE
       IMXER = I
    end do
    IWORK(16) = IMXER
    ! Set Y vector, T, and optional output. --------------------------------
580 CONTINUE
    CALL DCOPY(vstate % N, rwork % YH(:,1), 1, Y, 1)

    T = vstate % TN
    RWORK % condopt(11) = vstate % HU
    RWORK % condopt(12) = vstate % H
    RWORK % condopt(13) = vstate % TN
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
    
    ! -----------------------------------------------------------------------
    !  Block I.
    !  The following block handles all error returns due to illegal input
    !  (ISTATE = -3), as detected before calling the core integrator.
    !  First the error message routine is called.   If the illegal input
    !  is a negative ISTATE, the run is aborted (apparent infinite loop).
    ! -----------------------------------------------------------------------
    
601 continue
#ifndef CUDA    
    MSG = 'DVODE--  ISTATE (=I1) illegal '
    CALL XERRWD (MSG, 30, 1, 1, 1, ISTATE, 0, 0, ZERO, ZERO)
#endif
    IF (ISTATE .LT. 0) GO TO 800
    GO TO 700
602 continue
#ifndef CUDA    
    MSG = 'DVODE--  ITASK (=I1) illegal  '
    CALL XERRWD (MSG, 30, 2, 1, 1, ITASK, 0, 0, ZERO, ZERO)
#endif
    GO TO 700
603 continue
#ifndef CUDA    
    MSG='DVODE--  ISTATE (=I1) .gt. 1 but DVODE not initialized      '
    CALL XERRWD (MSG, 60, 3, 1, 1, ISTATE, 0, 0, ZERO, ZERO)
#endif    
    GO TO 700
604 continue
#ifndef CUDA    
    MSG = 'DVODE--  NEQ (=I1) .lt. 1     '
    CALL XERRWD (MSG, 30, 4, 1, 1, NEQ, 0, 0, ZERO, ZERO)
#endif
    GO TO 700
605 continue
#ifndef CUDA
    MSG = 'DVODE--  ISTATE = 3 and NEQ increased (I1 to I2)  '
    CALL XERRWD (MSG, 50, 5, 1, 2, vstate % N, NEQ, 0, ZERO, ZERO)
#endif
    GO TO 700
606 continue
#ifndef CUDA    
    MSG = 'DVODE--  ITOL (=I1) illegal   '
    CALL XERRWD (MSG, 30, 6, 1, 1, ITOL, 0, 0, ZERO, ZERO)
#endif
    GO TO 700
607 continue
#ifndef CUDA    
    MSG = 'DVODE--  IOPT (=I1) illegal   '
    CALL XERRWD (MSG, 30, 7, 1, 1, IOPT, 0, 0, ZERO, ZERO)
#endif
    GO TO 700
608 continue
#ifndef CUDA    
    MSG = 'DVODE--  MF (=I1) illegal     '
    CALL XERRWD (MSG, 30, 8, 1, 1, MF, 0, 0, ZERO, ZERO)
#endif
    GO TO 700
609 continue
#ifndef CUDA    
    MSG = 'DVODE--  ML (=I1) illegal:  .lt.0 or .ge.NEQ (=I2)'
    CALL XERRWD (MSG, 50, 9, 1, 2, ML, NEQ, 0, ZERO, ZERO)
#endif
    GO TO 700
610 continue
#ifndef CUDA    
    MSG = 'DVODE--  MU (=I1) illegal:  .lt.0 or .ge.NEQ (=I2)'
    CALL XERRWD (MSG, 50, 10, 1, 2, MU, NEQ, 0, ZERO, ZERO)
#endif    
    GO TO 700
611 continue
#ifndef CUDA    
    MSG = 'DVODE--  MAXORD (=I1) .lt. 0  '
    CALL XERRWD (MSG, 30, 11, 1, 1, vstate % MAXORD, 0, 0, ZERO, ZERO)
#endif    
    GO TO 700
612 continue
#ifndef CUDA    
    MSG = 'DVODE--  MXSTEP (=I1) .lt. 0  '
    CALL XERRWD (MSG, 30, 12, 1, 1, vstate % MXSTEP, 0, 0, ZERO, ZERO)
#endif
    GO TO 700
613 continue
#ifndef CUDA    
    MSG = 'DVODE--  MXHNIL (=I1) .lt. 0  '
    CALL XERRWD (MSG, 30, 13, 1, 1, vstate % MXHNIL, 0, 0, ZERO, ZERO)
#endif
    GO TO 700
614 continue
#ifndef CUDA    
    MSG = 'DVODE--  TOUT (=R1) behind T (=R2)      '
    CALL XERRWD (MSG, 40, 14, 1, 0, 0, 0, 2, TOUT, T)    
    MSG = '      integration direction is given by H0 (=R1)  '
    CALL XERRWD (MSG, 50, 14, 1, 0, 0, 0, 1, H0, ZERO)
#endif
    GO TO 700
615 continue
#ifndef CUDA    
    MSG = 'DVODE--  HMAX (=R1) .lt. 0.0  '
    CALL XERRWD (MSG, 30, 15, 1, 0, 0, 0, 1, HMAX, ZERO)
#endif
    GO TO 700
616 continue
#ifndef CUDA    
    MSG = 'DVODE--  HMIN (=R1) .lt. 0.0  '
    CALL XERRWD (MSG, 30, 16, 1, 0, 0, 0, 1, vstate % HMIN, ZERO)
#endif
    GO TO 700
617 CONTINUE
#ifndef CUDA    
    MSG='DVODE--  RWORK length needed, LENRW (=I1), exceeds LRW (=I2)'
    CALL XERRWD (MSG, 60, 17, 1, 2, LENRW, LRW, 0, ZERO, ZERO)
#endif
    GO TO 700
618 CONTINUE
#ifndef CUDA    
    MSG='DVODE--  IWORK length needed, LENIW (=I1), exceeds LIW (=I2)'
    CALL XERRWD (MSG, 60, 18, 1, 2, LENIW, LIW, 0, ZERO, ZERO)
#endif
    GO TO 700
619 continue
#ifndef CUDA    
    MSG = 'DVODE--  RTOL(I1) is R1 .lt. 0.0        '
    CALL XERRWD (MSG, 40, 19, 1, 1, I, 0, 1, RTOLI, ZERO)
#endif
    GO TO 700
620 continue
#ifndef CUDA    
    MSG = 'DVODE--  ATOL(I1) is R1 .lt. 0.0        '
    CALL XERRWD (MSG, 40, 20, 1, 1, I, 0, 1, ATOLI, ZERO)
#endif
    GO TO 700
621 continue
#ifndef CUDA    
    EWTI = rwork % ewt(I)
    MSG = 'DVODE--  EWT(I1) is R1 .le. 0.0         '
    CALL XERRWD (MSG, 40, 21, 1, 1, I, 0, 1, EWTI, ZERO)
#endif
    GO TO 700
622 CONTINUE
#ifndef CUDA    
    MSG='DVODE--  TOUT (=R1) too close to T(=R2) to start integration'
    CALL XERRWD (MSG, 60, 22, 1, 0, 0, 0, 2, TOUT, T)
#endif
    GO TO 700
623 CONTINUE
#ifndef CUDA    
    MSG='DVODE--  ITASK = I1 and TOUT (=R1) behind TCUR - HU (= R2)  '
    CALL XERRWD (MSG, 60, 23, 1, 1, ITASK, 0, 2, TOUT, TP)
#endif
    GO TO 700
624 CONTINUE
#ifndef CUDA    
    MSG='DVODE--  ITASK = 4 or 5 and TCRIT (=R1) behind TCUR (=R2)   '
    CALL XERRWD (MSG, 60, 24, 1, 0, 0, 0, 2, TCRIT, vstate % TN)
#endif
    GO TO 700
625 CONTINUE
#ifndef CUDA    
    MSG='DVODE--  ITASK = 4 or 5 and TCRIT (=R1) behind TOUT (=R2)   '
    CALL XERRWD (MSG, 60, 25, 1, 0, 0, 0, 2, TCRIT, TOUT)
#endif
    GO TO 700
626 continue
#ifndef CUDA    
    MSG = 'DVODE--  At start of problem, too much accuracy   '
    CALL XERRWD (MSG, 50, 26, 1, 0, 0, 0, 0, ZERO, ZERO)
    MSG='      requested for precision of machine:   see TOLSF (=R1) '
    CALL XERRWD (MSG, 60, 26, 1, 0, 0, 0, 1, TOLSF, ZERO)
#endif
    RWORK % condopt(14) = TOLSF
    GO TO 700
627 continue
#ifndef CUDA
    MSG='DVODE--  Trouble from DVINDY.  ITASK = I1, TOUT = R1.       '
    CALL XERRWD (MSG, 60, 27, 1, 1, ITASK, 0, 1, TOUT, ZERO)
#endif
    
700 CONTINUE
    ISTATE = -3

    return
    
800 continue
#ifndef CUDA    
    MSG = 'DVODE--  Run aborted:  apparent infinite loop     '
    CALL XERRWD (MSG, 50, 303, 2, 0, 0, 0, 0, ZERO, ZERO)
#endif
    return
  end subroutine dvode

  

#ifdef CUDA
  attributes(device) &
#endif  
  subroutine dvsol(WM, IWM, X, IERSL, vstate)

    !$acc routine seq
    
    ! -----------------------------------------------------------------------
    !  Call sequence input -- WM, IWM, X
    !  Call sequence output -- X, IERSL
    !  COMMON block variables accessed:
    !      /DVOD01/ -- H, RL1, MITER, N
    ! 
    !  Subroutines called by DVSOL: DGESL, DGBSL
    !  Function routines called by DVSOL: None
    ! -----------------------------------------------------------------------
    !  This routine manages the solution of the linear system arising from
    !  a chord iteration.  It is called if MITER .ne. 0.
    !  If MITER is 1 or 2, it calls DGESL to accomplish this.
    !  If MITER = 3 it updates the coefficient H*RL1 in the diagonal
    !  matrix, and then computes the solution.
    !  If MITER is 4 or 5, it calls DGBSL.
    !  Communication with DVSOL uses the following variables:
    !  WM    = Real work space containing the inverse diagonal matrix if
    !          MITER = 3 and the LU decomposition of the matrix otherwise.
    !          Storage of matrix elements starts at WM(3).
    !          WM also contains the following matrix-related data:
    !          WM(1) = SQRT(UROUND) (not used here),
    !          WM(2) = HRL1, the previous value of H*RL1, used if MITER = 3.
    !  IWM   = Integer work space containing pivot information, starting at
    !          IWM(31), if MITER is 1, 2, 4, or 5.  IWM also contains band
    !          parameters ML = IWM(1) and MU = IWM(2) if MITER is 4 or 5.
    !  X     = The right-hand side vector on input, and the solution vector
    !          on output, of length N.
    !  IERSL = Output flag.  IERSL = 0 if no trouble occurred.
    !          IERSL = 1 if a singular matrix arose with MITER = 3.
    ! -----------------------------------------------------------------------
    ! 

    implicit none

    real(dp_t) :: WM(:)
    real(dp_t) :: X(:)
    integer    :: IWM(:), IERSL
    type(dvode_t) :: vstate

    integer    :: I, MEBAND, ML, MU
    real(dp_t) :: DI, HRL1, PHRL1, R

    IERSL = 0
    GO TO (100, 100, 300, 400, 400), vstate % MITER
100 continue
    CALL DGESL (WM(3:3 + vstate % N**2 - 1), vstate % N, vstate % N, IWM(31:31 + vstate % N - 1), X, 0)
    RETURN

300 continue
    PHRL1 = WM(2)
    HRL1 = vstate % H*vstate % RL1
    WM(2) = HRL1
    IF (HRL1 .EQ. PHRL1) GO TO 330
    R = HRL1/PHRL1
    do I = 1,vstate % N
       DI = ONE - R*(ONE - ONE/WM(I+2))
       IF (ABS(DI) .EQ. ZERO) GO TO 390
       WM(I+2) = ONE/DI
    end do

330 continue
    do I = 1,vstate % N
       X(I) = WM(I+2)*X(I)
    end do
    RETURN
390 continue
    IERSL = 1
    RETURN

400 continue
    ML = IWM(1)
    MU = IWM(2)
    MEBAND = 2*ML + MU + 1
    CALL DGBSL (WM(3:3 + MEBAND * vstate % N - 1), MEBAND, vstate % N, &
         ML, MU, IWM(31:31 + vstate % N - 1), X, 0)
    RETURN
  end subroutine dvsol

#ifdef CUDA
  attributes(device) &
#endif  
  subroutine dacopy(NROW, NCOL, A, NROWA, B, NROWB)

    !$acc routine seq
    
    ! -----------------------------------------------------------------------
    !  Call sequence input -- NROW, NCOL, A, NROWA, NROWB
    !  Call sequence output -- B
    !  COMMON block variables accessed -- None
    ! 
    !  Subroutines called by DACOPY: DCOPY
    !  Function routines called by DACOPY: None
    ! -----------------------------------------------------------------------
    !  This routine copies one rectangular array, A, to another, B,
    !  where A and B may have different row dimensions, NROWA and NROWB.
    !  The data copied consists of NROW rows and NCOL columns.
    ! -----------------------------------------------------------------------

    implicit none

    integer    :: NROW, NCOL, NROWA, NROWB
    real(dp_t) :: A(NROWA,NCOL), B(NROWB,NCOL)
    integer    :: IC

    do IC = 1,NCOL
       CALL DCOPY (NROW, A(:,IC), 1, B(:,IC), 1)
    end do
    RETURN
  end subroutine dacopy

#ifdef CUDA
  attributes(device) &
#endif  
  subroutine dvjac(Y, IWM, IERPJ, RPAR, IPAR, rwork, vstate)

    !$acc routine seq
    
    ! -----------------------------------------------------------------------
    !  Call sequence input -- Y, YH, LDYH, EWT, FTEM, SAVF, WM, IWM,
    !                         F, JAC, RPAR, IPAR
    !  Call sequence output -- WM, IWM, IERPJ
    !  COMMON block variables accessed:
    !      /DVOD01/  CCMXJ, DRC, H, RL1, TN, UROUND, ICF, JCUR, LOCJS,
    !                MITER, MSBJ, N, NSLJ
    !      /DVOD02/  NFE, NST, NJE, NLU
    ! 
    !  Subroutines called by DVJAC: F, JAC, DACOPY, DCOPY, DGBFA, DGEFA,
    !                               DSCAL
    !  Function routines called by DVJAC: DVNORM
    ! -----------------------------------------------------------------------
    !  DVJAC is called by DVNLSD to compute and process the matrix
    !  P = I - h*rl1*J , where J is an approximation to the Jacobian.
    !  Here J is computed by the user-supplied routine JAC if
    !  MITER = 1 or 4, or by finite differencing if MITER = 2, 3, or 5.
    !  If MITER = 3, a diagonal approximation to J is used.
    !  If JSV = -1, J is computed from scratch in all cases.
    !  If JSV = 1 and MITER = 1, 2, 4, or 5, and if the saved value of J is
    !  considered acceptable, then P is constructed from the saved J.
    !  J is stored in wm and replaced by P.  If MITER .ne. 3, P is then
    !  subjected to LU decomposition in preparation for later solution
    !  of linear systems with P as coefficient matrix. This is done
    !  by DGEFA if MITER = 1 or 2, and by DGBFA if MITER = 4 or 5.
    ! 
    !  Communication with DVJAC is done with the following variables.  (For
    !  more details, please see the comments in the driver subroutine.)
    !  Y          = Vector containing predicted values on entry.
    !  YH         = The Nordsieck array, an LDYH by LMAX array, input.
    !  LDYH       = A constant .ge. N, the first dimension of YH, input.
    !  EWT        = An error weight vector of length N.
    !  SAVF       = Array containing f evaluated at predicted y, input.
    !  WM         = Real work space for matrices.  In the output, it containS
    !               the inverse diagonal matrix if MITER = 3 and the LU
    !               decomposition of P if MITER is 1, 2 , 4, or 5.
    !               Storage of matrix elements starts at WM(3).
    !               Storage of the saved Jacobian starts at WM(LOCJS).
    !               WM also contains the following matrix-related data:
    !               WM(1) = SQRT(UROUND), used in numerical Jacobian step.
    !               WM(2) = H*RL1, saved for later use if MITER = 3.
    !  IWM        = Integer work space containing pivot information,
    !               starting at IWM(31), if MITER is 1, 2, 4, or 5.
    !               IWM also contains band parameters ML = IWM(1) and
    !               MU = IWM(2) if MITER is 4 or 5.
    !  F          = Dummy name for the user supplied subroutine for f.
    !  JAC        = Dummy name for the user supplied Jacobian subroutine.
    !  RPAR, IPAR = Dummy names for user's real and integer work arrays.
    !  RL1        = 1/EL(2) (input).
    !  IERPJ      = Output error flag,  = 0 if no trouble, 1 if the P
    !               matrix is found to be singular.
    !  JCUR       = Output flag to indicate whether the Jacobian matrix
    !               (or approximation) is now current.
    !               JCUR = 0 means J is not current.
    !               JCUR = 1 means J is current.
    ! -----------------------------------------------------------------------
    ! 

    implicit none
  
    type(dvode_t) :: vstate
    type(rwork_t) :: rwork
    
    real(dp_t) :: Y(:)
    real(dp_t) :: RPAR(:)
    integer    :: IWM(:), IERPJ, IPAR(:)

    real(dp_t) :: CON, DI, FAC, HRL1, R, R0, SRUR, YI, YJ, YJJ
    integer    :: I, I1, I2, IER, II, J, J1, JJ, JOK, LENP, MBA, MBAND
    integer    :: MEB1, MEBAND, ML, ML3, MU, NP1

    ! Parameter declarations
    real(dp_t), parameter :: PT1 = 0.1D0

    IERPJ = 0
    HRL1 = vstate % H*vstate % RL1
    ! See whether J should be evaluated (JOK = -1) or not (JOK = 1). -------
    JOK = vstate % JSV
    IF (vstate % JSV .EQ. 1) THEN
       IF (vstate % NST .EQ. 0 .OR. vstate % NST .GT. vstate % NSLJ+vstate % MSBJ) JOK = -1
       IF (vstate % ICF .EQ. 1 .AND. vstate % DRC .LT. vstate % CCMXJ) JOK = -1
       IF (vstate % ICF .EQ. 2) JOK = -1
    ENDIF
    ! End of setting JOK. --------------------------------------------------

    IF (JOK .EQ. -1 .AND. vstate % MITER .EQ. 1) THEN
       ! If JOK = -1 and MITER = 1, call JAC to evaluate Jacobian. ------------
       vstate % NJE = vstate % NJE + 1
       vstate % NSLJ = vstate % NST
       vstate % JCUR = 1
       LENP = vstate % N * vstate % N
       do I = 1,LENP
          rwork % WM(I+2) = ZERO
       end do
       CALL JAC (vstate % N, vstate % TN, Y, 0, 0, &
            rwork % WM(3:3 + vstate % N**2 - 1), vstate % N, RPAR, IPAR)
       if (vstate % JSV .EQ. 1) then
          CALL DCOPY (LENP, rwork % WM(3:3 + LENP - 1), 1, &
               rwork % WM(vstate % LOCJS:vstate % LOCJS + LENP - 1), 1)
       endif
    ENDIF

    IF (JOK .EQ. -1 .AND. vstate % MITER .EQ. 2) THEN
       ! If MITER = 2, make N calls to F to approximate the Jacobian. ---------
       vstate % NJE = vstate % NJE + 1
       vstate % NSLJ = vstate % NST
       vstate % JCUR = 1
       FAC = DVNORM (vstate % N, rwork % SAVF, rwork % EWT)
       R0 = THOU*ABS(vstate % H) * vstate % UROUND * REAL(vstate % N)*FAC
       IF (R0 .EQ. ZERO) R0 = ONE
       SRUR = rwork % WM(1)
       J1 = 2
       do J = 1,vstate % N
          YJ = Y(J)
          R = MAX(SRUR*ABS(YJ),R0/rwork % EWT(J))
          Y(J) = Y(J) + R
          FAC = ONE/R
          CALL f_rhs (vstate % N, vstate % TN, Y, rwork % acor, RPAR, IPAR)
          do I = 1,vstate % N
             rwork % WM(I+J1) = (rwork % acor(I) - rwork % SAVF(I))*FAC
          end do
          Y(J) = YJ
          J1 = J1 + vstate % N
       end do
       vstate % NFE = vstate % NFE + vstate % N
       LENP = vstate % N * vstate % N
       if (vstate % JSV .EQ. 1) then
          CALL DCOPY (LENP, rwork % WM(3:3 + LENP - 1), 1, &
               rwork % WM(vstate % LOCJS:vstate % LOCJS + LENP - 1), 1)
       end if
    ENDIF

    IF (JOK .EQ. 1 .AND. (vstate % MITER .EQ. 1 .OR. vstate % MITER .EQ. 2)) THEN
       vstate % JCUR = 0
       LENP = vstate % N * vstate % N
       CALL DCOPY (LENP, rwork % WM(vstate % LOCJS:vstate % LOCJS + LENP - 1), 1, &
            rwork % WM(3:3 + LENP - 1), 1)
    ENDIF

    IF (vstate % MITER .EQ. 1 .OR. vstate % MITER .EQ. 2) THEN
       ! Multiply Jacobian by scalar, add identity, and do LU decomposition. --
       CON = -HRL1
       CALL DSCAL (LENP, CON, rwork % WM(3:3 + LENP - 1), 1)
       J = 3
       NP1 = vstate % N + 1
       do I = 1,vstate % N
          rwork % WM(J) = rwork % WM(J) + ONE
          J = J + NP1
       end do
       vstate % NLU = vstate % NLU + 1
       CALL DGEFA (rwork % WM(3:3 + vstate % N**2 - 1), vstate % N, &
            vstate % N, IWM(31:31 + vstate % N - 1), IER)
       IF (IER .NE. 0) IERPJ = 1
       RETURN
    ENDIF
    ! End of code block for MITER = 1 or 2. --------------------------------

    IF (vstate % MITER .EQ. 3) THEN
       ! If MITER = 3, construct a diagonal approximation to J and P. ---------
       vstate % NJE = vstate % NJE + 1
       vstate % JCUR = 1
       rwork % WM(2) = HRL1
       R = vstate % RL1*PT1
       do I = 1,vstate % N
          Y(I) = Y(I) + R*(vstate % H * rwork % SAVF(I) - rwork % YH(I,2))
       end do
       CALL f_rhs (vstate % N, vstate % TN, Y, &
            rwork % WM(3:3 + vstate % N - 1), RPAR, IPAR)
       vstate % NFE = vstate % NFE + 1
       do I = 1,vstate % N
          R0 = vstate % H * rwork % SAVF(I) - rwork % YH(I,2)
          DI = PT1*R0 - vstate % H*(rwork % WM(I+2) - rwork % SAVF(I))
          rwork % WM(I+2) = ONE
          IF (ABS(R0) .LT. vstate % UROUND/rwork % EWT(I)) cycle
          IF (ABS(DI) .EQ. ZERO) GO TO 330
          rwork % WM(I+2) = PT1*R0/DI
       end do
       RETURN
330    continue
       IERPJ = 1
       RETURN
    ENDIF
    ! End of code block for MITER = 3. -------------------------------------

    ! Set constants for MITER = 4 or 5. ------------------------------------
    ML = IWM(1)
    MU = IWM(2)
    ML3 = ML + 3
    MBAND = ML + MU + 1
    MEBAND = MBAND + ML
    LENP = MEBAND * vstate % N

    if (JOK .EQ. -1 .AND. vstate % MITER .EQ. 4) then
       ! If JOK = -1 and MITER = 4, call JAC to evaluate Jacobian. ------------
       vstate % NJE = vstate % NJE + 1
       vstate % NSLJ = vstate % NST
       vstate % JCUR = 1
       do I = 1,LENP
          rwork % WM(I+2) = ZERO
       end do
       CALL JAC (vstate % N, vstate % TN, Y, ML, MU, rwork % WM(ML3:ML3 + MEBAND * vstate % N - 1), MEBAND, RPAR, IPAR)
       if (vstate % JSV .EQ. 1) then
          CALL DACOPY(MBAND, vstate % N, &
               rwork % WM(ML3:ML3 + MEBAND * vstate % N - 1), MEBAND, &
               rwork % WM(vstate % LOCJS:vstate % LOCJS + MBAND * vstate % N - 1), MBAND)
       end if

    else if (JOK .EQ. -1 .AND. vstate % MITER .EQ. 5) then
       ! If MITER = 5, make ML+MU+1 calls to F to approximate the Jacobian. ---
       vstate % NJE = vstate % NJE + 1
       vstate % NSLJ = vstate % NST
       vstate % JCUR = 1
       MBA = MIN(MBAND,vstate % N)
       MEB1 = MEBAND - 1
       SRUR = rwork % WM(1)
       FAC = DVNORM (vstate % N, rwork % SAVF, rwork % EWT)
       R0 = THOU*ABS(vstate % H) * vstate % UROUND * REAL(vstate % N)*FAC
       IF (R0 .EQ. ZERO) R0 = ONE
       do J = 1,MBA
          do I = J,vstate % N,MBAND
             YI = Y(I)
             R = MAX(SRUR*ABS(YI),R0/rwork % EWT(I))
             Y(I) = Y(I) + R
          end do
          CALL f_rhs (vstate % N, vstate % TN, Y, rwork % acor, RPAR, IPAR)
          do JJ = J,vstate % N,MBAND
             Y(JJ) = rwork % YH(JJ,1)
             YJJ = Y(JJ)
             R = MAX(SRUR*ABS(YJJ),R0/rwork % EWT(JJ))
             FAC = ONE/R
             I1 = MAX(JJ-MU,1)
             I2 = MIN(JJ+ML,vstate % N)
             II = JJ*MEB1 - ML + 2
             do I = I1,I2
                rwork % WM(II+I) = (rwork % acor(I) - rwork % SAVF(I))*FAC
             end do
          end do
       end do
       vstate % NFE = vstate % NFE + MBA
       if (vstate % JSV .EQ. 1) then
          CALL DACOPY(MBAND, vstate % N, &
               rwork % WM(ML3:ML3 + MEBAND * vstate % N - 1), MEBAND, &
               rwork % WM(vstate % LOCJS:vstate % LOCJS + MBAND * vstate % N - 1), MBAND)
       end if
    end if

    IF (JOK .EQ. 1) THEN
       vstate % JCUR = 0
       CALL DACOPY(MBAND, vstate % N, &
            rwork % WM(vstate % LOCJS:vstate % LOCJS + MBAND * vstate % N - 1), MBAND, &
            rwork % WM(ML3:ML3 + MEBAND * vstate % N - 1), MEBAND)
    ENDIF

    ! Multiply Jacobian by scalar, add identity, and do LU decomposition.
    CON = -HRL1
    CALL DSCAL (LENP, CON, rwork % WM(3:3 + LENP - 1), 1 )
    II = MBAND + 2
    do I = 1,vstate % N
       rwork % WM(II) = rwork % WM(II) + ONE
       II = II + MEBAND
    end do
    vstate % NLU = vstate % NLU + 1
    CALL DGBFA (rwork % WM(3:3 + MEBAND * vstate % N - 1), MEBAND, vstate % N, ML, &
         MU, IWM(31:31 + vstate % N - 1), IER)
    if (IER .NE. 0) then
       IERPJ = 1
    end if
    RETURN
    ! End of code block for MITER = 4 or 5. --------------------------------
  end subroutine dvjac

#ifdef CUDA
  attributes(device) &
#endif  
  subroutine dvnlsd(Y, IWM, NFLAG, RPAR, IPAR, rwork, vstate)

    !$acc routine seq
    
    ! -----------------------------------------------------------------------
    !  Call sequence input -- Y, YH, LDYH, SAVF, EWT, ACOR, IWM, WM,
    !                         F, JAC, NFLAG, RPAR, IPAR
    !  Call sequence output -- YH, ACOR, WM, IWM, NFLAG
    !  COMMON block variables accessed:
    !      /DVOD01/ ACNRM, CRATE, DRC, H, RC, RL1, TQ(5), TN, ICF,
    !                 JCUR, METH, MITER, N, NSLP
    !      /DVOD02/ HU, NCFN, NETF, NFE, NJE, NLU, NNI, NQU, NST
    ! 
    !  Subroutines called by DVNLSD: F, DAXPY, DCOPY, DSCAL, DVJAC, DVSOL
    !  Function routines called by DVNLSD: DVNORM
    ! -----------------------------------------------------------------------
    !  Subroutine DVNLSD is a nonlinear system solver, which uses functional
    !  iteration or a chord (modified Newton) method.  For the chord method
    !  direct linear algebraic system solvers are used.  Subroutine DVNLSD
    !  then handles the corrector phase of this integration package.
    ! 
    !  Communication with DVNLSD is done with the following variables. (For
    !  more details, please see the comments in the driver subroutine.)
    ! 
    !  Y          = The dependent variable, a vector of length N, input.
    !  YH         = The Nordsieck (Taylor) array, LDYH by LMAX, input
    !               and output.  On input, it contains predicted values.
    !  LDYH       = A constant .ge. N, the first dimension of YH, input.
    !  VSAV       = Unused work array.
    !  SAVF       = A work array of length N.
    !  EWT        = An error weight vector of length N, input.
    !  ACOR       = A work array of length N, used for the accumulated
    !               corrections to the predicted y vector.
    !  WM,IWM     = Real and integer work arrays associated with matrix
    !               operations in chord iteration (MITER .ne. 0).
    !  F          = Dummy name for user supplied routine for f.
    !  JAC        = Dummy name for user supplied Jacobian routine.
    !  PDUM       = Unused dummy subroutine name.  Included for uniformity
    !               over collection of integrators.
    !  NFLAG      = Input/output flag, with values and meanings as follows:
    !               INPUT
    !                   0 first call for this time step.
    !                  -1 convergence failure in previous call to DVNLSD.
    !                  -2 error test failure in DVSTEP.
    !               OUTPUT
    !                   0 successful completion of nonlinear solver.
    !                  -1 convergence failure or singular matrix.
    !                  -2 unrecoverable error in matrix preprocessing
    !                     (cannot occur here).
    !                  -3 unrecoverable error in solution (cannot occur
    !                     here).
    !  RPAR, IPAR = Dummy names for user's real and integer work arrays.
    ! 
    !  IPUP       = Own variable flag with values and meanings as follows:
    !               0,            do not update the Newton matrix.
    !               MITER .ne. 0, update Newton matrix, because it is the
    !                             initial step, order was changed, the error
    !                             test failed, or an update is indicated by
    !                             the scalar RC or step counter NST.
    ! 
    !  For more details, see comments in driver subroutine.
    ! -----------------------------------------------------------------------
    !

    implicit none
  
    type(dvode_t) :: vstate
    type(rwork_t) :: rwork
    real(dp_t)    :: Y(vstate % N)
    real(dp_t) :: RPAR(:)
    integer    :: IWM(:), NFLAG, IPAR(:)
    
    real(dp_t) :: CSCALE, DCON, DEL, DELP
    integer    :: I, IERPJ, IERSL, M

    ! Parameter declarations
    real(dp_t), parameter :: CCMAX = 0.3D0
    real(dp_t), parameter :: CRDOWN = 0.3D0
    real(dp_t), parameter :: RDIV  = 2.0D0
    integer, parameter :: MAXCOR = 3
    integer, parameter :: MSBP = 20

    ! -----------------------------------------------------------------------
    !  On the first step, on a change of method order, or after a
    !  nonlinear convergence failure with NFLAG = -2, set IPUP = MITER
    !  to force a Jacobian update when MITER .ne. 0.
    ! -----------------------------------------------------------------------
    IF (vstate % JSTART .EQ. 0) vstate % NSLP = 0
    IF (NFLAG .EQ. 0) vstate % ICF = 0
    IF (NFLAG .EQ. -2) vstate % IPUP = vstate % MITER
    IF ( (vstate % JSTART .EQ. 0) .OR. (vstate % JSTART .EQ. -1) ) vstate % IPUP = vstate % MITER
    ! If this is functional iteration, set CRATE .eq. 1 and drop to 220
    IF (vstate % MITER .EQ. 0) THEN
       vstate % CRATE = ONE
       GO TO 220
    ENDIF
    ! -----------------------------------------------------------------------
    !  RC is the ratio of new to old values of the coefficient H/EL(2)=h/l1.
    !  When RC differs from 1 by more than CCMAX, IPUP is set to MITER
    !  to force DVJAC to be called, if a Jacobian is involved.
    !  In any case, DVJAC is called at least every MSBP steps.
    ! -----------------------------------------------------------------------
    vstate % DRC = ABS(vstate % RC-ONE)
    IF (vstate % DRC .GT. CCMAX .OR. vstate % NST .GE. vstate % NSLP+MSBP) vstate % IPUP = vstate % MITER
    ! -----------------------------------------------------------------------
    !  Up to MAXCOR corrector iterations are taken.  A convergence test is
    !  made on the r.m.s. norm of each correction, weighted by the error
    !  weight vector EWT.  The sum of the corrections is accumulated in the
    !  vector ACOR(i).  The YH array is not altered in the corrector loop.
    ! -----------------------------------------------------------------------
220 continue
    M = 0
    DELP = ZERO

    CALL DCOPY(vstate % N, rwork % yh(:,1), 1, Y, 1)
    CALL f_rhs (vstate % N, vstate % TN, Y, rwork % savf, RPAR, IPAR)
    vstate % NFE = vstate % NFE + 1
    IF (vstate % IPUP .LE. 0) GO TO 250
    ! -----------------------------------------------------------------------
    !  If indicated, the matrix P = I - h*rl1*J is reevaluated and
    !  preprocessed before starting the corrector iteration.  IPUP is set
    !  to 0 as an indicator that this has been done.
    ! -----------------------------------------------------------------------
    CALL DVJAC (Y, IWM, IERPJ, RPAR, IPAR, rwork, vstate)
    vstate % IPUP = 0
    vstate % RC = ONE
    vstate % DRC = ZERO
    vstate % CRATE = ONE
    vstate % NSLP = vstate % NST
    ! If matrix is singular, take error return to force cut in step size. --
    IF (IERPJ .NE. 0) GO TO 430
250 continue
    do I = 1,vstate % N
       rwork % acor(I) = ZERO
    end do
    ! This is a looping point for the corrector iteration. -----------------
270 IF (vstate % MITER .NE. 0) GO TO 350
    ! -----------------------------------------------------------------------
    !  In the case of functional iteration, update Y directly from
    !  the result of the last function evaluation.
    ! -----------------------------------------------------------------------
    do I = 1,vstate % N
       rwork % SAVF(I) = vstate % RL1*(vstate % H * rwork % SAVF(I) - rwork % YH(I,2))
    end do
    do I = 1,vstate % N
       Y(I) = rwork % SAVF(I) - rwork % ACOR(I)
    end do
    DEL = DVNORM (vstate % N, Y, rwork % EWT)
    do I = 1,vstate % N
       Y(I) = rwork % YH(I,1) + rwork % SAVF(I)
    end do
    CALL DCOPY(vstate % N, rwork % SAVF, 1, rwork % ACOR, 1)

    GO TO 400
    ! -----------------------------------------------------------------------
    !  In the case of the chord method, compute the corrector error,
    !  and solve the linear system with that as right-hand side and
    !  P as coefficient matrix.  The correction is scaled by the factor
    !  2/(1+RC) to account for changes in h*rl1 since the last DVJAC call.
    ! -----------------------------------------------------------------------
350 continue
    do I = 1,vstate % N
       Y(I) = (vstate % RL1*vstate % H) * rwork % SAVF(I) - &
            (vstate % RL1 * rwork % YH(I,2) + rwork % ACOR(I))
    end do
    CALL DVSOL (rwork % wm, IWM, Y, IERSL, vstate)
    vstate % NNI = vstate % NNI + 1
    IF (IERSL .GT. 0) GO TO 410
    IF (vstate % METH .EQ. 2 .AND. vstate % RC .NE. ONE) THEN
       CSCALE = TWO/(ONE + vstate % RC)
       CALL DSCAL (vstate % N, CSCALE, Y, 1)
    ENDIF
    DEL = DVNORM (vstate % N, Y, rwork % EWT)
    call daxpy(vstate % N, ONE, Y, 1, rwork % acor, 1)

    do I = 1,vstate % N
       Y(I) = rwork % YH(I,1) + rwork % ACOR(I)
    end do
    ! -----------------------------------------------------------------------
    !  Test for convergence.  If M .gt. 0, an estimate of the convergence
    !  rate constant is stored in CRATE, and this is used in the test.
    ! -----------------------------------------------------------------------
400 continue
    IF (M .NE. 0) vstate % CRATE = MAX(CRDOWN*vstate % CRATE,DEL/DELP)
    DCON = DEL*MIN(ONE,vstate % CRATE)/vstate % TQ(4)
    IF (DCON .LE. ONE) GO TO 450
    M = M + 1
    IF (M .EQ. MAXCOR) GO TO 410
    IF (M .GE. 2 .AND. DEL .GT. RDIV*DELP) GO TO 410
    DELP = DEL
    CALL f_rhs (vstate % N, vstate % TN, Y, rwork % SAVF, RPAR, IPAR)
    vstate % NFE = vstate % NFE + 1
    GO TO 270
    
410 IF (vstate % MITER .EQ. 0 .OR. vstate % JCUR .EQ. 1) GO TO 430
    vstate % ICF = 1
    vstate % IPUP = vstate % MITER
    GO TO 220
    
430 CONTINUE
    NFLAG = -1
    vstate % ICF = 2
    vstate % IPUP = vstate % MITER
    RETURN

    ! Return for successful step. ------------------------------------------
450 continue
    NFLAG = 0
    vstate % JCUR = 0
    vstate % ICF = 0
    IF (M .EQ. 0) vstate % ACNRM = DEL
    IF (M .GT. 0) vstate % ACNRM = DVNORM (vstate % N, rwork % ACOR, rwork % EWT)
    RETURN
  end subroutine dvnlsd

#ifdef CUDA  
  attributes(device) &
#endif  
  subroutine dvjust(IORD, rwork, vstate)

    !$acc routine seq
    
    ! -----------------------------------------------------------------------
    !  Call sequence input -- YH, LDYH, IORD
    !  Call sequence output -- YH
    !  COMMON block input -- NQ, METH, LMAX, HSCAL, TAU(13), N
    !  COMMON block variables accessed:
    !      /DVOD01/ -- HSCAL, TAU(13), LMAX, METH, N, NQ,
    ! 
    !  Subroutines called by DVJUST: DAXPY
    !  Function routines called by DVJUST: None
    ! -----------------------------------------------------------------------
    !  This subroutine adjusts the YH array on reduction of order,
    !  and also when the order is increased for the stiff option (METH = 2).
    !  Communication with DVJUST uses the following:
    !  IORD  = An integer flag used when METH = 2 to indicate an order
    !          increase (IORD = +1) or an order decrease (IORD = -1).
    !  HSCAL = Step size H used in scaling of Nordsieck array YH.
    !          (If IORD = +1, DVJUST assumes that HSCAL = TAU(1).)
    !  See References 1 and 2 for details.
    ! -----------------------------------------------------------------------
    !

    implicit none
  
    type(dvode_t) :: vstate
    type(rwork_t) :: rwork
    
    integer    :: IORD

    real(dp_t) :: ALPH0, ALPH1, HSUM, PROD, T1, XI,XIOLD
    integer    :: I, IBACK, J, JP1, LP1, NQM1, NQM2, NQP1

    IF ((vstate % NQ .EQ. 2) .AND. (IORD .NE. 1)) RETURN
    NQM1 = vstate % NQ - 1
    NQM2 = vstate % NQ - 2
    GO TO (100, 200), vstate % METH
    ! -----------------------------------------------------------------------
    !  Nonstiff option...
    !  Check to see if the order is being increased or decreased.
    ! -----------------------------------------------------------------------
100 CONTINUE
    IF (IORD .EQ. 1) GO TO 180
    ! Order decrease. ------------------------------------------------------
    do J = 1, vstate % LMAX
       vstate % EL(J) = ZERO
    end do
    vstate % EL(2) = ONE
    HSUM = ZERO
    do J = 1, NQM2
       ! Construct coefficients of x*(x+xi(1))*...*(x+xi(j)). -----------------
       HSUM = HSUM + vstate % TAU(J)
       XI = HSUM/vstate % HSCAL
       JP1 = J + 1
       do IBACK = 1, JP1
          I = (J + 3) - IBACK
          vstate % EL(I) = vstate % EL(I)*XI + vstate % EL(I-1)
       end do
    end do
    ! Construct coefficients of integrated polynomial. ---------------------
    do J = 2, NQM1
       vstate % EL(J+1) = REAL(vstate % NQ) * vstate % EL(J)/REAL(J)
    end do
    ! Subtract correction terms from YH array. -----------------------------
    do J = 3, vstate % NQ
       do I = 1, vstate % N
          rwork % YH(I,J) = rwork % YH(I,J) - &
               rwork % YH(I,vstate % L) * vstate % EL(J)
       end do
    end do
    RETURN
    ! Order increase. ------------------------------------------------------
    ! Zero out next column in YH array. ------------------------------------
180 CONTINUE
    LP1 = vstate % L + 1
    do I = 1, vstate % N
       rwork % YH(I,LP1) = ZERO
    end do
    RETURN
    ! -----------------------------------------------------------------------
    !  Stiff option...
    !  Check to see if the order is being increased or decreased.
    ! -----------------------------------------------------------------------
200 CONTINUE
    IF (IORD .EQ. 1) GO TO 300
    ! Order decrease. ------------------------------------------------------
    do J = 1, vstate % LMAX
       vstate % EL(J) = ZERO
    end do
    vstate % EL(3) = ONE
    HSUM = ZERO
    do J = 1,NQM2
       ! Construct coefficients of x*x*(x+xi(1))*...*(x+xi(j)). ---------------
       HSUM = HSUM + vstate % TAU(J)
       XI = HSUM/vstate % HSCAL
       JP1 = J + 1
       do IBACK = 1, JP1
          I = (J + 4) - IBACK
          vstate % EL(I) = vstate % EL(I) * XI + vstate % EL(I-1)
       end do
    end do
    ! Subtract correction terms from YH array. -----------------------------
    do J = 3,vstate % NQ
       do I = 1, vstate % N
          rwork % YH(I,J) = rwork % YH(I,J) - &
               rwork % YH(I,vstate % L) * vstate % EL(J)
       end do
    end do
    RETURN
    ! Order increase. ------------------------------------------------------
300 continue
    do J = 1, vstate % LMAX
       vstate % EL(J) = ZERO
    end do
    vstate % EL(3) = ONE
    ALPH0 = -ONE
    ALPH1 = ONE
    PROD = ONE
    XIOLD = ONE
    HSUM = vstate % HSCAL
    IF (vstate % NQ .EQ. 1) GO TO 340
    do J = 1, NQM1
       ! Construct coefficients of x*x*(x+xi(1))*...*(x+xi(j)). ---------------
       JP1 = J + 1
       HSUM = HSUM + vstate % TAU(JP1)
       XI = HSUM/vstate % HSCAL
       PROD = PROD*XI
       ALPH0 = ALPH0 - ONE/REAL(JP1)
       ALPH1 = ALPH1 + ONE/XI
       do IBACK = 1, JP1
          I = (J + 4) - IBACK
          vstate % EL(I) = vstate % EL(I) * XIOLD + vstate % EL(I-1)
       end do
       XIOLD = XI
    end do
340 CONTINUE
    T1 = (-ALPH0 - ALPH1)/PROD
    ! Load column L + 1 in YH array. ---------------------------------------
    LP1 = vstate % L + 1
    do I = 1, vstate % N
       rwork % YH(I,LP1) = T1 * rwork % YH(I,vstate % LMAX)
    end do
    ! Add correction terms to YH array. ------------------------------------
    NQP1 = vstate % NQ + 1
    do J = 3, NQP1
       CALL DAXPY(vstate % N, vstate % EL(J), &
            rwork % YH(1:vstate % N, LP1), 1, rwork % YH(1:vstate % N, J), 1)
    end do
    RETURN
  end subroutine dvjust

#ifdef CUDA  
  attributes(device) &
#endif  
  subroutine dvset(vstate)

    !$acc routine seq
    
    ! -----------------------------------------------------------------------
    !  Call sequence communication: None
    !  COMMON block variables accessed:
    !      /DVOD01/ -- EL(13), H, TAU(13), TQ(5), L(= NQ + 1),
    !                  METH, NQ, NQWAIT
    ! 
    !  Subroutines called by DVSET: None
    !  Function routines called by DVSET: None
    ! -----------------------------------------------------------------------
    !  DVSET is called by DVSTEP and sets coefficients for use there.
    ! 
    !  For each order NQ, the coefficients in EL are calculated by use of
    !   the generating polynomial lambda(x), with coefficients EL(i).
    !       lambda(x) = EL(1) + EL(2)*x + ... + EL(NQ+1)*(x**NQ).
    !  For the backward differentiation formulas,
    !                                      NQ-1
    !       lambda(x) = (1 + x/xi*(NQ)) * product (1 + x/xi(i) ) .
    !                                      i = 1
    !  For the Adams formulas,
    !                               NQ-1
    !       (d/dx) lambda(x) = c * product (1 + x/xi(i) ) ,
    !                               i = 1
    !       lambda(-1) = 0,    lambda(0) = 1,
    !  where c is a normalization constant.
    !  In both cases, xi(i) is defined by
    !       H*xi(i) = t sub n  -  t sub (n-i)
    !               = H + TAU(1) + TAU(2) + ... TAU(i-1).
    ! 
    ! 
    !  In addition to variables described previously, communication
    !  with DVSET uses the following:
    !    TAU    = A vector of length 13 containing the past NQ values
    !             of H.
    !    EL     = A vector of length 13 in which vset stores the
    !             coefficients for the corrector formula.
    !    TQ     = A vector of length 5 in which vset stores constants
    !             used for the convergence test, the error test, and the
    !             selection of H at a new order.
    !    METH   = The basic method indicator.
    !    NQ     = The current order.
    !    L      = NQ + 1, the length of the vector stored in EL, and
    !             the number of columns of the YH array being used.
    !    NQWAIT = A counter controlling the frequency of order changes.
    !             An order change is about to be considered if NQWAIT = 1.
    ! -----------------------------------------------------------------------
    !

    implicit none
  
    type(dvode_t), intent(inout) :: vstate
    real(dp_t) :: EM(13)
    real(dp_t) :: AHATN0, ALPH0, CNQM1, CSUM, ELP
    real(dp_t) :: EM0, FLOTI, FLOTL, FLOTNQ, HSUM, RXI, RXIS, S
    real(dp_t) :: T1, T2, T3, T4, T5, T6, XI
    integer    :: I, IBACK, J, JP1, NQM1, NQM2

    ! Parameter declaration
    real(dp_t), parameter :: CORTES = 0.1D0

    FLOTL = REAL(vstate % L)
    NQM1 = vstate % NQ - 1
    NQM2 = vstate % NQ - 2
    GO TO (100, 200), vstate % METH

    !  Set coefficients for Adams methods. ----------------------------------
100 IF (vstate % NQ .NE. 1) GO TO 110
    vstate % EL(1) = ONE
    vstate % EL(2) = ONE
    vstate % TQ(1) = ONE
    vstate % TQ(2) = TWO
    vstate % TQ(3) = SIX*vstate % TQ(2)
    vstate % TQ(5) = ONE
    GO TO 300
110 continue
    HSUM = vstate % H
    EM(1) = ONE
    FLOTNQ = FLOTL - ONE
    do I = 2, vstate % L
       EM(I) = ZERO
    end do
    do J = 1, NQM1
       IF ((J .NE. NQM1) .OR. (vstate % NQWAIT .NE. 1)) GO TO 130
       S = ONE
       CSUM = ZERO
       do I = 1, NQM1
          CSUM = CSUM + S*EM(I)/REAL(I+1)
          S = -S
       end do
       vstate % TQ(1) = EM(NQM1)/(FLOTNQ*CSUM)
130    continue
       RXI = vstate % H/HSUM
       do IBACK = 1, J
          I = (J + 2) - IBACK
          EM(I) = EM(I) + EM(I-1)*RXI
       end do
       HSUM = HSUM + vstate % TAU(J)
    end do
    ! Compute integral from -1 to 0 of polynomial and of x times it. -------
    S = ONE
    EM0 = ZERO
    CSUM = ZERO
    do I = 1, vstate % NQ
       FLOTI = REAL(I)
       EM0 = EM0 + S*EM(I)/FLOTI
       CSUM = CSUM + S*EM(I)/(FLOTI+ONE)
       S = -S
    end do
    ! In EL, form coefficients of normalized integrated polynomial. --------
    S = ONE/EM0
    vstate % EL(1) = ONE
    do I = 1, vstate % NQ
       vstate % EL(I+1) = S*EM(I)/REAL(I)
    end do
    XI = HSUM/vstate % H
    vstate % TQ(2) = XI*EM0/CSUM
    vstate % TQ(5) = XI/vstate % EL(vstate % L)
    IF (vstate % NQWAIT .NE. 1) GO TO 300
    ! For higher order control constant, multiply polynomial by 1+x/xi(q). -
    RXI = ONE/XI
    do IBACK = 1, vstate % NQ
       I = (vstate % L + 1) - IBACK
       EM(I) = EM(I) + EM(I-1)*RXI
    end do
    ! Compute integral of polynomial. --------------------------------------
    S = ONE
    CSUM = ZERO
    do I = 1, vstate % L
       CSUM = CSUM + S*EM(I)/REAL(I+1)
       S = -S
    end do
    vstate % TQ(3) = FLOTL*EM0/CSUM
    GO TO 300

    ! Set coefficients for BDF methods. ------------------------------------
200 continue
    do I = 3, vstate % L
       vstate % EL(I) = ZERO
    end do
    vstate % EL(1) = ONE
    vstate % EL(2) = ONE
    ALPH0 = -ONE
    AHATN0 = -ONE
    HSUM = vstate % H
    RXI = ONE
    RXIS = ONE
    IF (vstate % NQ .EQ. 1) GO TO 240
    do J = 1, NQM2
       ! In EL, construct coefficients of (1+x/xi(1))*...*(1+x/xi(j+1)). ------
       HSUM = HSUM + vstate % TAU(J)
       RXI = vstate % H/HSUM
       JP1 = J + 1
       ALPH0 = ALPH0 - ONE/REAL(JP1)
       do IBACK = 1, JP1
          I = (J + 3) - IBACK
          vstate % EL(I) = vstate % EL(I) + vstate % EL(I-1)*RXI
       end do
    end do
    ALPH0 = ALPH0 - ONE/REAL(vstate % NQ)
    RXIS = -vstate % EL(2) - ALPH0
    HSUM = HSUM + vstate % TAU(NQM1)
    RXI = vstate % H/HSUM
    AHATN0 = -vstate % EL(2) - RXI
    do IBACK = 1, vstate % NQ
       I = (vstate % NQ + 2) - IBACK
       vstate % EL(I) = vstate % EL(I) + vstate % EL(I-1)*RXIS
    end do
240 continue
    T1 = ONE - AHATN0 + ALPH0
    T2 = ONE + REAL(vstate % NQ)*T1
    vstate % TQ(2) = ABS(ALPH0*T2/T1)
    vstate % TQ(5) = ABS(T2/(vstate % EL(vstate % L)*RXI/RXIS))
    IF (vstate % NQWAIT .NE. 1) GO TO 300
    CNQM1 = RXIS/vstate % EL(vstate % L)
    T3 = ALPH0 + ONE/REAL(vstate % NQ)
    T4 = AHATN0 + RXI
    ELP = T3/(ONE - T4 + T3)
    vstate % TQ(1) = ABS(ELP/CNQM1)
    HSUM = HSUM + vstate % TAU(vstate % NQ)
    RXI = vstate % H/HSUM
    T5 = ALPH0 - ONE/REAL(vstate % NQ+1)
    T6 = AHATN0 - RXI
    ELP = T2/(ONE - T6 + T5)
    vstate % TQ(3) = ABS(ELP*RXI*(FLOTL + ONE)*T5)
300 continue
    vstate % TQ(4) = CORTES*vstate % TQ(2)
    RETURN
  end subroutine dvset

#ifdef CUDA  
  attributes(device) &
#endif  
  subroutine dvstep(Y, IWM, RPAR, IPAR, rwork, vstate)

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
    !  Y      = An array of length N used for the dependent variable vector.
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

    implicit none

    type(dvode_t) :: vstate
    type(rwork_t) :: rwork
    real(dp_t) :: Y(vstate % N)
    real(dp_t) :: RPAR(:)
    real(dp_t) :: yhscratch(vstate % N * vstate % LMAX)
    integer    :: IWM(:), IPAR(:)
    
    real(dp_t) :: CNQUOT, DDN, DSM, DUP, TOLD
    real(dp_t) :: ETAQ, ETAQM1, ETAQP1, FLOTL, R
    integer    :: I, I1, I2, IBACK, J, JB, NCF, NFLAG

    ! Parameter declarations
    integer, parameter :: KFC = -3
    integer, parameter :: KFH = -7
    integer, parameter :: MXNCF = 10
    real(dp_t), parameter :: ADDON = 1.0D-6
    real(dp_t), parameter :: BIAS1 = 6.0D0
    real(dp_t), parameter :: BIAS2 = 6.0D0
    real(dp_t), parameter :: BIAS3 = 10.0D0
    real(dp_t), parameter :: ETACF = 0.25D0
    real(dp_t), parameter :: ETAMIN = 0.1D0
    real(dp_t), parameter :: ETAMXF = 0.2D0
    real(dp_t), parameter :: ETAMX1 = 1.0D4
    real(dp_t), parameter :: ETAMX2 = 10.0D0
    real(dp_t), parameter :: ETAMX3 = 10.0D0
    real(dp_t), parameter :: ONEPSM = 1.00001D0
    real(dp_t), parameter :: THRESH = 1.5D0

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
    vstate % NQNYH = vstate % NQ * vstate % NYH
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
    IF (vstate % N .EQ. vstate % NYH) GO TO 120
    I1 = 1 + (vstate % NEWQ + 1)*vstate % NYH
    I2 = (vstate % MAXORD + 1)*vstate % NYH
    IF (I1 .GT. I2) GO TO 120
    rwork % YH(:, vstate % NEWQ + 1:) = ZERO
120 IF (vstate % NEWQ .LE. vstate % MAXORD) GO TO 140
    FLOTL = REAL(vstate % LMAX)
    IF (vstate % MAXORD .LT. vstate % NQ-1) THEN
       DDN = DVNORM (vstate % N, rwork % SAVF, rwork % EWT)/vstate % TQ(1)
       vstate % ETA = ONE/((BIAS1*DDN)**(ONE/FLOTL) + ADDON)
    ENDIF
    IF (vstate % MAXORD .EQ. vstate % NQ .AND. vstate % NEWQ .EQ. vstate % NQ+1) vstate % ETA = ETAQ
    IF (vstate % MAXORD .EQ. vstate % NQ-1 .AND. vstate % NEWQ .EQ. vstate % NQ+1) THEN
       vstate % ETA = ETAQM1
       CALL DVJUST (-1, rwork, vstate)
    ENDIF
    IF (vstate % MAXORD .EQ. vstate % NQ-1 .AND. vstate % NEWQ .EQ. vstate % NQ) THEN
       DDN = DVNORM (vstate % N, rwork % SAVF, rwork % EWT)/vstate % TQ(1)
       vstate % ETA = ONE/((BIAS1*DDN)**(ONE/FLOTL) + ADDON)
       CALL DVJUST (-1, rwork, vstate)
    ENDIF
    vstate % ETA = MIN(vstate % ETA,ONE)
    vstate % NQ = vstate % MAXORD
    vstate % L = vstate % LMAX
140 continue
    IF (vstate % KUTH .EQ. 1) vstate % ETA = MIN(vstate % ETA,ABS(vstate % H/vstate % HSCAL))
    IF (vstate % KUTH .EQ. 0) vstate % ETA = MAX(vstate % ETA,vstate % HMIN/ABS(vstate % HSCAL))
    vstate % ETA = vstate % ETA/MAX(ONE,ABS(vstate % HSCAL)*vstate % HMXI*vstate % ETA)
    vstate % NEWH = 1
    vstate % NQWAIT = vstate % L
    IF (vstate % NEWQ .LE. vstate % MAXORD) GO TO 50
    ! Rescale the history array for a change in H by a factor of ETA. ------
150 continue
    R = ONE

    do J = 2, vstate % L
       R = R * vstate % ETA
       CALL DSCAL(vstate % N, R, rwork % YH(1:vstate % N,J), 1)
    end do
    vstate % H = vstate % HSCAL * vstate % ETA
    vstate % HSCAL = vstate % H
    vstate % RC = vstate % RC * vstate % ETA
    vstate % NQNYH = vstate % NQ*vstate % NYH
    ! -----------------------------------------------------------------------
    !  This section computes the predicted values by effectively
    !  multiplying the YH array by the Pascal triangle matrix.
    !  DVSET is called to calculate all integration coefficients.
    !  RC is the ratio of new to old values of the coefficient H/EL(2)=h/l1.
    ! -----------------------------------------------------------------------
200 continue
    vstate % TN = vstate % TN + vstate % H
    !DON - begin
    ! optimize this multiplication and check the math as it makes no sense
    ! Make the 2-D YH array into 1-D yhscratch
    do I = 1, vstate % NYH
       do JB = 1, vstate % LMAX
          I1 = I + (JB-1) * vstate % NYH
          yhscratch(I1) = rwork % YH(I, JB)
       end do
    end do
    I1 = vstate % NQNYH + 1
    do JB = 1, vstate % NQ
       I1 = I1 - vstate % NYH
       do I = I1, vstate % NQNYH
          yhscratch(I) = yhscratch(I) + yhscratch(I+vstate % NYH)
       end do
    end do
    ! Make the 1-D yhscratch array into 2-D YH
    do I = 1, vstate % NYH
       do JB = 1, vstate % LMAX
          I1 = I + (JB-1) * vstate % NYH
          rwork % YH(I, JB) = yhscratch(I1)
       end do
    end do
    !DON - end
    CALL DVSET(vstate)
    vstate % RL1 = ONE/vstate % EL(2)
    vstate % RC = vstate % RC * (vstate % RL1/vstate % PRL1)
    vstate % PRL1 = vstate % RL1
    ! 
    !  Call the nonlinear system solver. ------------------------------------
    !
    CALL dvnlsd (Y, IWM, NFLAG, RPAR, IPAR, rwork, vstate)

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
    !DON - begin
    ! optimize this multiplication and check the math as it makes no sense
    ! Make the 2-D YH array into 1-D yhscratch
    do I = 1, vstate % NYH
       do JB = 1, vstate % LMAX
          I1 = I + (JB-1) * vstate % NYH
          yhscratch(I1) = rwork % YH(I, JB)
       end do
    end do
    I1 = vstate % NQNYH + 1
    do JB = 1, vstate % NQ
       I1 = I1 - vstate % NYH
       do I = I1, vstate % NQNYH
          yhscratch(I) = yhscratch(I) - yhscratch(I+vstate % NYH)
       end do
    end do
    ! Make the 1-D yhscratch array into 2-D YH
    do I = 1, vstate % NYH
       do JB = 1, vstate % LMAX
          I1 = I + (JB-1) * vstate % NYH
          rwork % YH(I, JB) = yhscratch(I1)
       end do
    end do
    !DON - end
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
       CALL DAXPY(vstate % N, vstate % EL(J), rwork % acor, 1, rwork % yh(:,J), 1)
    end do
    vstate % NQWAIT = vstate % NQWAIT - 1
    IF ((vstate % L .EQ. vstate % LMAX) .OR. (vstate % NQWAIT .NE. 1)) GO TO 490
    CALL DCOPY(vstate % N, rwork % acor, 1, rwork % yh(:,vstate % LMAX), 1)
    
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
    !DON - begin
    ! optimize this multiplication and check the math as it makes no sense
    ! Make the 2-D YH array into 1-D yhscratch
    do I = 1, vstate % NYH
       do JB = 1, vstate % LMAX
          I1 = I + (JB-1) * vstate % NYH
          yhscratch(I1) = rwork % YH(I, JB)
       end do
    end do
    I1 = vstate % NQNYH + 1
    do JB = 1, vstate % NQ
       I1 = I1 - vstate % NYH
       do I = I1, vstate % NQNYH
          yhscratch(I) = yhscratch(I) - yhscratch(I+vstate % NYH)
       end do
    end do
    ! Make the 1-D yhscratch array into 2-D YH
    do I = 1, vstate % NYH
       do JB = 1, vstate % LMAX
          I1 = I + (JB-1) * vstate % NYH
          rwork % YH(I, JB) = yhscratch(I1)
       end do
    end do
    !DON - end
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
    CALL f_rhs (vstate % N, vstate % TN, Y, rwork % savf, RPAR, IPAR)
    vstate % NFE = vstate % NFE + 1
    do I = 1, vstate % N
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
    DDN = DVNORM (vstate % N, rwork % yh(:,vstate % L), rwork % ewt)/vstate % TQ(1)
    ETAQM1 = ONE/((BIAS1*DDN)**(ONE/(FLOTL - ONE)) + ADDON)
570 continue
    ETAQP1 = ZERO
    IF (vstate % L .EQ. vstate % LMAX) GO TO 580
    ! Compute ratio of new H to current H at current order plus one. -------
    CNQUOT = (vstate % TQ(5)/vstate % CONP)*(vstate % H/vstate % TAU(2))**vstate % L
    do I = 1, vstate % N
       rwork % savf(I) = rwork % acor(I) - CNQUOT * rwork % yh(I,vstate % LMAX)
    end do
    DUP = DVNORM (vstate % N, rwork % savf, rwork % ewt)/vstate % TQ(3)
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
    CALL DCOPY(vstate % N, rwork % acor, 1, rwork % yh(:,vstate % LMAX), 1)
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
    CALL DSCAL(vstate % N, R, rwork % acor, 1)

720 continue
    vstate % JSTART = 1
    RETURN
  end subroutine dvstep
      
end module dvode_module
