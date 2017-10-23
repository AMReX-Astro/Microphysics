module dvode_dumach_module

  use dvode_constants_module

  implicit none

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

    use bl_types, only: dp_t

    implicit none
  
    real(dp_t) :: U, dum
    dum = EPSILON(U)
  end function dumach

end module dvode_dumach_module
