module cuvode_dvjac_module

  use cuvode_parameters_module, only: VODE_LMAX, VODE_NEQS, VODE_LIW, VODE_MAXORD
  use cuvode_types_module, only: dvode_t
  use amrex_fort_module, only: rt => amrex_real
  use linpack_module
  use cuvode_dacopy_module

  use cuvode_constants_module

  implicit none

contains

#if defined(AMREX_USE_CUDA) && !defined(AMREX_USE_GPU_PRAGMA)
  attributes(device) &
#endif
  subroutine dvjac(pivot, IERPJ, vstate)

    !$acc routine seq
    
    ! -----------------------------------------------------------------------
    !  DVJAC is called by DVNLSD to compute and process the matrix
    !  P = I - h*rl1*J , where J is an approximation to the Jacobian.
    !  Here J is computed by the user-supplied routine JAC if
    !  MITER = 1 or 4, or by finite differencing if MITER = 2, 3, or 5.
    !  If MITER = 3, a diagonal approximation to J is used.
    !  If JSV = -1, J is computed from scratch in all cases.
    !  If JSV = 1 and MITER = 1, 2, 4, or 5, and if the saved value of J is
    !  considered acceptable, then P is constructed from the saved J.
    !  J is replaced by P.  P is then subjected to LU decomposition
    ! in preparation for later solution of linear systems with P as
    ! coefficient matrix. This is done by DGEFA.
    ! 
    !  Communication with DVJAC is done with the following variables.  (For
    !  more details, please see the comments in the driver subroutine.)
    !  Y          = Array containing predicted values on entry.
    !  YH         = The Nordsieck array, an LDYH by LMAX array, input.
    !  LDYH       = A constant .ge. N, the first dimension of YH, input.
    !  EWT        = An error weight array of length N.
    !  SAVF       = Array containing f evaluated at predicted y, input.
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
    integer,       intent(inout) :: pivot(VODE_NEQS)
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
          vstate % JAC(I) = ZERO
       end do
       CALL JAC (vstate % TN, vstate, 0, 0, vstate % jac, VODE_NEQS)
       if (vstate % JSV .EQ. 1) then
          do i = 1, LENP
             vstate % jac_save(i) = vstate % jac(i)
          end do
       endif
    ENDIF

    IF (JOK .EQ. -1 .AND. vstate % MITER .EQ. 2) THEN
       ! If MITER = 2, make N calls to F to approximate the Jacobian. ---------
       vstate % NJE = vstate % NJE + 1
       vstate % NSLJ = vstate % NST
       vstate % JCUR = 1
       FAC = DVNORM (vstate % SAVF, vstate % EWT)
       R0 = THOU*ABS(vstate % H) * vstate % UROUND * REAL(VODE_NEQS)*FAC
       IF (R0 .EQ. ZERO) R0 = ONE
       J1 = 0
       do J = 1,VODE_NEQS
          YJ = vstate % Y(J)
          R = MAX(vstate % SRUR * ABS(YJ), R0 / vstate % EWT(J))
          vstate % Y(J) = vstate % Y(J) + R
          FAC = ONE/R
          CALL f_rhs (vstate % TN, vstate, vstate % acor)
          do I = 1,VODE_NEQS
             vstate % jac(I+J1) = (vstate % acor(I) - vstate % SAVF(I))*FAC
          end do
          vstate % Y(J) = YJ
          J1 = J1 + VODE_NEQS
       end do
       vstate % NFE = vstate % NFE + VODE_NEQS
       LENP = VODE_NEQS * VODE_NEQS
       if (vstate % JSV .EQ. 1) then
          do i = 1, LENP
             vstate % jac_save(i) = vstate % jac(i)
          end do
       end if
    ENDIF

    IF (JOK .EQ. 1) THEN
       vstate % JCUR = 0
       LENP = VODE_NEQS * VODE_NEQS
       do i = 1, LENP
           vstate % jac(i) = vstate % jac_save(i)
       end do
    ENDIF

    ! Multiply Jacobian by scalar, add identity, and do LU decomposition. --
    CON = -HRL1
    vstate % jac(:) = vstate % jac(:) * CON
    J = 1
    NP1 = VODE_NEQS + 1
    do I = 1,VODE_NEQS
       vstate % jac(J) = vstate % jac(J) + ONE
       J = J + NP1
    end do
    vstate % NLU = vstate % NLU + 1
    CALL DGEFA (vstate % jac, pivot, IER)
    IF (IER .NE. 0) IERPJ = 1

  end subroutine dvjac

end module cuvode_dvjac_module
