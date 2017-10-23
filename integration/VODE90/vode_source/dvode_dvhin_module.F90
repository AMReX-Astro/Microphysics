module dvode_dvhin_module

  use dvode_constants_module
  use dvode_dvnorm_module

  implicit none

contains

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

    use rpar_indices
    use vode_rhs_module, only: f_rhs
    use bl_types, only: dp_t
  
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

end module dvode_dvhin_module
