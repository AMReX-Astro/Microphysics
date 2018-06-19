module dvode_dewset_module
  
  use vode_type_module, only: rwork_t
  use vode_parameters_module, only: VODE_LMAX, VODE_NEQS, VODE_LIW,   &
                                    VODE_LENWM, VODE_MAXORD, VODE_ITOL
  use dvode_type_module, only: dvode_t

  use dvode_constants_module

  implicit none

contains

  AMREX_DEVICE subroutine dewset(vstate, rwork)

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

    implicit none

    ! Declare arguments
    type(dvode_t), intent(in   ) :: vstate
    type(rwork_t), intent(inout) :: rwork

    ! Declare local variables
    integer    :: I, N

    GO TO (10, 20, 30, 40), VODE_ITOL
10  CONTINUE
    do I = 1,VODE_NEQS
       rwork % EWT(I) = vstate % RTOL(1)*ABS(rwork % YH(I,1)) + vstate % ATOL(1)
    end do
    RETURN
20  CONTINUE
    do I = 1,VODE_NEQS
       rwork % EWT(I) = vstate % RTOL(1)*ABS(rwork % YH(I,1)) + vstate % ATOL(I)
    end do
    RETURN
30  CONTINUE
    do I = 1,VODE_NEQS
       rwork % EWT(I) = vstate % RTOL(I)*ABS(rwork % YH(I,1)) + vstate % ATOL(1)
    end do
    RETURN
40  CONTINUE

    do I = 1,VODE_NEQS
       rwork % EWT(I) = vstate % RTOL(I)*ABS(rwork % YH(I,1)) + vstate % ATOL(I)
    end do
    RETURN
  end subroutine dewset

end module dvode_dewset_module
