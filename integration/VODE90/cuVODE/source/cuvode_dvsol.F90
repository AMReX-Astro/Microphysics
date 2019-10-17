module cuvode_dvsol_module

  use cuvode_parameters_module, only: VODE_LMAX, VODE_NEQS, VODE_LIW,   &
                                      VODE_LENWM, VODE_MAXORD, VODE_ITOL
  use cuvode_types_module, only: dvode_t
  use amrex_fort_module, only: rt => amrex_real
  use linpack_module

  use cuvode_constants_module

  implicit none

contains

#if defined(AMREX_USE_CUDA) && !defined(AMREX_USE_GPU_PRAGMA)
  attributes(device) &
#endif
  subroutine dvsol(WM, IWM, IERSL, vstate)

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
    !  X     = The right-hand side array on input, and the solution array
    !          on output, of length N.
    !  IERSL = Output flag.  IERSL = 0 if no trouble occurred.
    !          IERSL = 1 if a singular matrix arose with MITER = 3.
    ! -----------------------------------------------------------------------
    ! 

    implicit none

    ! Declare arguments
    type(dvode_t), intent(inout) :: vstate
    real(rt),    intent(inout) :: WM(VODE_LENWM)
    integer,       intent(in   ) :: IWM(VODE_LIW)
    integer,       intent(  out) :: IERSL

    ! Declare local variables
    integer    :: I, MEBAND, ML, MU
    real(rt) :: DI, HRL1, PHRL1, R

    !$gpu

    IERSL = 0
    GO TO (100, 100, 300, 400, 400), vstate % MITER
100 continue
    CALL DGESL (WM(3:3 + VODE_NEQS**2 - 1), IWM(31:31 + VODE_NEQS - 1), vstate % Y(:))
    RETURN

300 continue
    PHRL1 = WM(2)
    HRL1 = vstate % H*vstate % RL1
    WM(2) = HRL1
    IF (HRL1 .EQ. PHRL1) GO TO 330
    R = HRL1/PHRL1
    do I = 1,VODE_NEQS
       DI = ONE - R*(ONE - ONE/WM(I+2))
       IF (ABS(DI) .EQ. ZERO) GO TO 390
       WM(I+2) = ONE/DI
    end do

330 continue
    do I = 1,VODE_NEQS
       vstate % Y(I) = WM(I+2)*vstate % Y(I)
    end do
    RETURN
390 continue
    IERSL = 1
    RETURN

400 continue
    ML = IWM(1)
    MU = IWM(2)
    MEBAND = 2*ML + MU + 1
    CALL DGBSL (WM(3:3 + MEBAND * VODE_NEQS - 1), MEBAND, &
         ML, MU, IWM(31:31 + VODE_NEQS - 1), vstate % Y(:), 0)
    RETURN
  end subroutine dvsol

end module cuvode_dvsol_module
