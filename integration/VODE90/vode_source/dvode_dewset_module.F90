module dvode_dewset_module

  use dvode_constants_module

  implicit none

contains

#ifdef CUDA
  attributes(device) &
#endif  
  subroutine dewset(N, ITOL, RTOL, ATOL, YCUR, EWT)

    !$acc routine seq
    
    ! ***BEGIN PROLOGUE  DEWSET
    ! ***SUBSIDIARY
    ! ***PURPOSE  Set error weight vector.
    ! ***TYPE      DOUBLE PRECISION (SEWSET-S, DEWSET-D)
    ! ***AUTHOR  Hindmarsh, Alan C., (LLNL)
    ! ***DESCRIPTION
    ! 
    !   This subroutine sets the error weight vector EWT according to
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

    use bl_types, only: dp_t

    implicit none
  
    integer    :: I, N, ITOL
    real(dp_t) :: RTOL(N), ATOL(N)
    real(dp_t) :: YCUR(N), EWT(N)

    GO TO (10, 20, 30, 40), ITOL
10  CONTINUE
    do I = 1,N
       EWT(I) = RTOL(1)*ABS(YCUR(I)) + ATOL(1)
    end do
    RETURN
20  CONTINUE
    do I = 1,N
       EWT(I) = RTOL(1)*ABS(YCUR(I)) + ATOL(I)
    end do
    RETURN
30  CONTINUE
    do I = 1,N
       EWT(I) = RTOL(I)*ABS(YCUR(I)) + ATOL(1)
    end do
    RETURN
40  CONTINUE

    do I = 1,N
       EWT(I) = RTOL(I)*ABS(YCUR(I)) + ATOL(I)
    end do
    RETURN
  end subroutine dewset

end module dvode_dewset_module
