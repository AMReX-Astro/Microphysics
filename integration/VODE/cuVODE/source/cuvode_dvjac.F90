module cuvode_dvjac_module

  use cuvode_parameters_module, only: VODE_LMAX, VODE_NEQS, VODE_LIW,   &
                                      VODE_LENWM, VODE_MAXORD, VODE_ITOL
  use cuvode_types_module, only: dvode_t, rwork_t
  use amrex_fort_module, only: rt => amrex_real
  use linpack_module
  use cuvode_dacopy_module

  use cuvode_constants_module

  implicit none

contains

#if defined(AMREX_USE_CUDA) && !defined(AMREX_USE_GPU_PRAGMA)
  attributes(device) &
#endif
  subroutine dvjac(IWM, IERPJ, rwork, vstate)

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
    !  Y          = Array containing predicted values on entry.
    !  YH         = The Nordsieck array, an LDYH by LMAX array, input.
    !  LDYH       = A constant .ge. N, the first dimension of YH, input.
    !  EWT        = An error weight array of length N.
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
    integer,       intent(  out) :: IERPJ

    ! Declare local variables
    real(rt) :: CON, DI, FAC, HRL1, R, R0, SRUR, YI, YJ, YJJ
    integer    :: I, I1, I2, IER, II, J, J1, JJ, JOK, LENP, MBA, MBAND
    integer    :: MEB1, MEBAND, ML, ML3, MU, NP1

    ! Parameter declarations
    real(rt), parameter :: PT1 = 0.1e0_rt

    !$gpu

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
       LENP = VODE_NEQS * VODE_NEQS
       do I = 1,LENP
          rwork % WM(I+2) = ZERO
       end do
       CALL JAC (vstate % TN, vstate % Y, 0, 0, &
            rwork % WM(3:3 + VODE_NEQS**2 - 1), VODE_NEQS, vstate % RPAR)
       if (vstate % JSV .EQ. 1) then
          do I = 0, LENP-1
             rwork % WM(vstate % LOCJS + I) = rwork % WM(3 + I)
          end do
       endif
    ENDIF

    IF (JOK .EQ. -1 .AND. vstate % MITER .EQ. 2) THEN
       ! If MITER = 2, make N calls to F to approximate the Jacobian. ---------
       vstate % NJE = vstate % NJE + 1
       vstate % NSLJ = vstate % NST
       vstate % JCUR = 1
       FAC = DVNORM (rwork % SAVF, rwork % EWT)
       R0 = THOU*ABS(vstate % H) * vstate % UROUND * REAL(VODE_NEQS)*FAC
       IF (R0 .EQ. ZERO) R0 = ONE
       SRUR = rwork % WM(1)
       J1 = 2
       do J = 1,VODE_NEQS
          YJ = vstate % Y(J)
          R = MAX(SRUR*ABS(YJ),R0/rwork % EWT(J))
          vstate % Y(J) = vstate % Y(J) + R
          FAC = ONE/R
          CALL f_rhs (vstate % TN, vstate % Y, rwork % acor, vstate % RPAR)
          do I = 1,VODE_NEQS
             rwork % WM(I+J1) = (rwork % acor(I) - rwork % SAVF(I))*FAC
          end do
          vstate % Y(J) = YJ
          J1 = J1 + VODE_NEQS
       end do
       vstate % NFE = vstate % NFE + VODE_NEQS
       LENP = VODE_NEQS * VODE_NEQS
       if (vstate % JSV .EQ. 1) then
          do I = 0, LENP-1
             rwork % WM(vstate % LOCJS + I) = rwork % WM(3 + I)
          end do
       end if
    ENDIF

    IF (JOK .EQ. 1 .AND. (vstate % MITER .EQ. 1 .OR. vstate % MITER .EQ. 2)) THEN
       vstate % JCUR = 0
       LENP = VODE_NEQS * VODE_NEQS
       do I = 0, LENP - 1
           rwork % WM(3 + I) = rwork % WM(vstate % LOCJS + I)
       end do
    ENDIF

    IF (vstate % MITER .EQ. 1 .OR. vstate % MITER .EQ. 2) THEN
       ! Multiply Jacobian by scalar, add identity, and do LU decomposition. --
       CON = -HRL1
       rwork % WM(3:3+LENP-1) = rwork % WM(3:3+LENP-1) * CON
       J = 3
       NP1 = VODE_NEQS + 1
       do I = 1,VODE_NEQS
          rwork % WM(J) = rwork % WM(J) + ONE
          J = J + NP1
       end do
       vstate % NLU = vstate % NLU + 1
       CALL DGEFA (rwork % WM(3:3 + VODE_NEQS**2 - 1), IWM(31:31 + VODE_NEQS - 1), IER)
       IF (IER .NE. 0) IERPJ = 1
       RETURN
    ENDIF
    ! End of code block for MITER = 1 or 2. --------------------------------

  end subroutine dvjac

end module cuvode_dvjac_module
