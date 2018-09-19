module dvode_dvhin_module

  use vode_rhs_module, only: f_rhs, jac
  use vode_type_module, only: rwork_t
  use vode_parameters_module, only: VODE_LMAX, VODE_NEQS, VODE_LIW,   &
                                    VODE_LENWM, VODE_MAXORD, VODE_ITOL
  use dvode_type_module, only: dvode_t
  use amrex_fort_module, only: rt => amrex_real

  use dvode_constants_module

  implicit none

contains

  subroutine dvhin(vstate, rwork, H0, NITER, IER)
  
    !$acc routine seq
    
    ! -----------------------------------------------------------------------
    !  Call sequence input -- N, T0, Y0, YDOT, F, RPAR, TOUT, UROUND,
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
    !  Y0     = Array of initial conditions, input.
    !  YDOT   = Array of initial first derivatives, input.
    !  [NOTE: Y0 = YH(:,1) and YDOT = YH(:,2) in this subroutine now]
    !
    !  F      = Name of subroutine for right-hand side f(t,y), input.
    !  RPAR = Dummy names for user's real and integer work arrays.
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
    !
    ! Note: the following variable replacements have been made ---
    !    t0 = vstate % t
    !    yh = rwork % yh
    !    rpar = vstate % rpar
    !    tout = vstate % tout
    !    uround = vstate % uround
    !    ewt = rwork % ewt
    !    itol = VODE_ITOL
    !    atol = vstate % atol
    !    y = vstate % y
    !    temp = rwork % acor
  
    use vode_rhs_module, only: f_rhs
    use dvode_dvnorm_module, only: dvnorm ! function

    implicit none

    ! Declare arguments
    type(dvode_t), intent(inout) :: vstate
    type(rwork_t), intent(inout) :: rwork
    real(rt), intent(inout) :: H0
    integer,    intent(  out) :: NITER, IER

    ! Declare local variables
    real(rt) :: AFI, ATOLI, DELYI, H, HG, HLB, HNEW, HRAT
    real(rt) :: HUB, T1, TDIST, TROUND, YDDNRM
    integer    :: I, ITER

    real(rt), parameter :: PT1 = 0.1D0

    !$gpu

    NITER = 0
    TDIST = ABS(vstate % TOUT - vstate % T)
    TROUND = vstate % UROUND*MAX(ABS(vstate % T),ABS(vstate % TOUT))
    IF (TDIST .LT. TWO*TROUND) GO TO 100

    ! Set a lower bound on h based on the roundoff level in vstate % T and vstate % TOUT. ---
    HLB = HUN*TROUND
    ! Set an upper bound on h based on vstate % TOUT-vstate % T and the initial Y and YDOT. -
    HUB = PT1*TDIST
    ATOLI = vstate % ATOL(1)
    
    do I = 1, VODE_NEQS
       IF (VODE_ITOL .EQ. 2 .OR. VODE_ITOL .EQ. 4) ATOLI = vstate % ATOL(I)
       DELYI = PT1*ABS(rwork % YH(I,1)) + ATOLI
       AFI = ABS(rwork % YH(I,2))
       IF (AFI*HUB .GT. DELYI) HUB = DELYI/AFI
    end do

    ! Set initial guess for h as geometric mean of upper and lower bounds. -
    ITER = 0
    HG = SQRT(HLB*HUB)
    ! If the bounds have crossed, exit with the mean value. ----------------
    IF (HUB .LT. HLB) THEN
       H0 = HG
       GO TO 90
    ENDIF

    ! Looping point for iteration. -----------------------------------------
50  CONTINUE
    ! Estimate the second derivative as a difference quotient in f. --------
    H = SIGN (HG, vstate % TOUT - vstate % T)
    T1 = vstate % T + H
    do I = 1, VODE_NEQS
       vstate % Y(I) = rwork % YH(I,1) + H*rwork % YH(I,2)
    end do
    CALL f_rhs(T1, vstate % Y, rwork % ACOR, vstate % RPAR)
    do I = 1, VODE_NEQS
       rwork % ACOR(I) = (rwork % ACOR(I) - rwork % YH(I,2))/H
    end do
    YDDNRM = DVNORM(rwork % ACOR, rwork % EWT)
    ! Get the corresponding new value of h. --------------------------------
    IF (YDDNRM*HUB*HUB .GT. TWO) THEN
       HNEW = SQRT(TWO/YDDNRM)
    ELSE
       HNEW = SQRT(HG*HUB)
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
    H0 = SIGN(H0, vstate % TOUT - vstate % T)
    NITER = ITER
    IER = 0
    RETURN
    ! Error return for vstate % TOUT - vstate % T too small. --------------------------------
100 continue
    IER = -1
    RETURN
  end subroutine dvhin

end module dvode_dvhin_module
