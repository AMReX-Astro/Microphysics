module cuvode_dewset_module
  
  use cuvode_parameters_module, only: VODE_LMAX, VODE_NEQS, VODE_LIW,   &
                                      VODE_LENWM, VODE_MAXORD, VODE_ITOL
  use cuvode_types_module, only: dvode_t, rwork_t

  use cuvode_constants_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

contains

#if defined(AMREX_USE_CUDA) && !defined(AMREX_USE_GPU_PRAGMA)
  attributes(device) &
#endif
  subroutine dewset(vstate, rwork)

    !$acc routine seq
    
    ! ***BEGIN PROLOGUE  DEWSET
    ! ***SUBSIDIARY
    ! ***PURPOSE  Set error weight array.
    ! ***TYPE      DOUBLE PRECISION (SEWSET-S, DEWSET-D)
    ! ***AUTHOR  Hindmarsh, Alan C., (LLNL)
    ! ***DESCRIPTION
    ! 
    !   This subroutine sets the error weight array EWT according to
    !       EWT(i) = RTOL(i)*ABS(YCUR(i)) + ATOL(i),  i = 1,...,N,
    !   with the subscript on RTOL and/or ATOL possibly replaced by 1 above,
    !   depending on the value of ITOL.
    ! 
    ! ***SEE ALSO  DLSODE
    ! ***ROUTINES CALLED  (NONE)
    ! ***REVISION HISTORY  (YYMMDD)
    !    791129  DATE WRITTEN
    !    890501  Modified prologue to SLATEC/LDOC format.  (FNF)
    !    890503  Minor cosmetic changes.  (FNF)
    !    930809  Renamed to allow single/double precision versions. (ACH)
    ! ***END PROLOGUE  DEWSET
    ! **End

    implicit none

    ! Declare arguments
    type(dvode_t), intent(in   ) :: vstate
    type(rwork_t), intent(inout) :: rwork

    ! Declare local variables
    integer    :: I, N

    !$gpu

    if (VODE_ITOL == 1) then
       do I = 1,VODE_NEQS
          rwork % EWT(I) = vstate % RTOL(1)*ABS(rwork % YH(I,1)) + vstate % ATOL(1)
       end do
    else if (VODE_ITOL == 2) then
       do I = 1,VODE_NEQS
          rwork % EWT(I) = vstate % RTOL(1)*ABS(rwork % YH(I,1)) + vstate % ATOL(I)
       end do
    else if (VODE_ITOL == 3) then
       do I = 1,VODE_NEQS
          rwork % EWT(I) = vstate % RTOL(I)*ABS(rwork % YH(I,1)) + vstate % ATOL(1)
       end do
    else if (VODE_ITOL == 4) then
       do I = 1,VODE_NEQS
          rwork % EWT(I) = vstate % RTOL(I)*ABS(rwork % YH(I,1)) + vstate % ATOL(I)
       end do
    end if

    RETURN
  end subroutine dewset

end module cuvode_dewset_module
