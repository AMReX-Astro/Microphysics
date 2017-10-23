module dvode_dvsol_module

  use dvode_constants_module
  use linpack_module

  implicit none

contains

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

    use dvode_type_module, only: dvode_t
    use bl_types, only: dp_t

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

end module dvode_dvsol_module
