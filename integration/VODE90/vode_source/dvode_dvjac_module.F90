module dvode_dvjac_module

  use vode_rhs_module, only: f_rhs, jac
  use dvode_dacopy_module
  use dvode_constants_module
  use blas_module
  use linpack_module
  use vode_type_module, only: rwork_t
  use dvode_type_module, only: dvode_t
  use bl_types, only: dp_t
  use dvode_dvnorm_module, only: dvnorm
  
  implicit none

contains

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

end module dvode_dvjac_module
