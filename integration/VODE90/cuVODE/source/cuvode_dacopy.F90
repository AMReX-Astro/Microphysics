module cuvode_dacopy_module

  use amrex_fort_module, only: rt => amrex_real
  use cuvode_constants_module

  implicit none

contains

#if defined(AMREX_USE_CUDA) && !defined(AMREX_USE_GPU_PRAGMA)
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

    ! Declare arguments
    integer,     intent(in   ) :: NROW, NCOL, NROWA, NROWB
    real(rt),  intent(in   ) :: A(NROWA,NCOL)
    real(rt),  intent(inout) :: B(NROWB,NCOL)

    ! Declare local variables
    integer    :: IC

    !$gpu

    do IC = 1,NCOL
       B(1:NROW,IC) = A(1:NROW,IC)
    end do
    RETURN
  end subroutine dacopy

end module cuvode_dacopy_module
