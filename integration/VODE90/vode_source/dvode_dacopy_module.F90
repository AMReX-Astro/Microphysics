module dvode_dacopy_module

  use dvode_constants_module

  implicit none

contains

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
    !  Subroutines called by DACOPY: DCOPYN
    !  Function routines called by DACOPY: None
    ! -----------------------------------------------------------------------
    !  This routine copies one rectangular array, A, to another, B,
    !  where A and B may have different row dimensions, NROWA and NROWB.
    !  The data copied consists of NROW rows and NCOL columns.
    ! -----------------------------------------------------------------------

    use bl_types, only: dp_t
    use blas_module
    use linpack_module

    implicit none

    integer    :: NROW, NCOL, NROWA, NROWB
    real(dp_t) :: A(NROWA,NCOL), B(NROWB,NCOL)
    integer    :: IC

    do IC = 1,NCOL
       CALL DCOPYN (NROW, A(:,IC), 1, B(:,IC), 1)
    end do
    RETURN
  end subroutine dacopy

end module dvode_dacopy_module
