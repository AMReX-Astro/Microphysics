module cuvode_dvjust_module

  use cuvode_parameters_module, only: VODE_LMAX, VODE_NEQS, VODE_LIW,   &
                                      VODE_LENWM, VODE_MAXORD, VODE_ITOL
  use cuvode_types_module, only: dvode_t, rwork_t
  use amrex_fort_module, only: rt => amrex_real
  use blas_module

  use cuvode_constants_module

  implicit none

contains

#if defined(AMREX_USE_CUDA) && !defined(AMREX_USE_GPU_PRAGMA)
  attributes(device) &
#endif
  subroutine dvjust(IORD, rwork, vstate)

    !$acc routine seq
    
    ! -----------------------------------------------------------------------
    !  Call sequence input -- YH, LDYH, IORD
    !  Call sequence output -- YH
    !  COMMON block input -- NQ, METH, LMAX, HSCAL, TAU(13), N
    !  COMMON block variables accessed:
    !      /DVOD01/ -- HSCAL, TAU(13), LMAX, METH, N, NQ,
    ! 
    !  Subroutines called by DVJUST: DAXPY
    !  Function routines called by DVJUST: None
    ! -----------------------------------------------------------------------
    !  This subroutine adjusts the YH array on reduction of order,
    !  and also when the order is increased for the stiff option (METH = 2).
    !  Communication with DVJUST uses the following:
    !  IORD  = An integer flag used when METH = 2 to indicate an order
    !          increase (IORD = +1) or an order decrease (IORD = -1).
    !  HSCAL = Step size H used in scaling of Nordsieck array YH.
    !          (If IORD = +1, DVJUST assumes that HSCAL = TAU(1).)
    !  See References 1 and 2 for details.
    ! -----------------------------------------------------------------------
    !

    implicit none

    ! Declare arguments
    type(dvode_t), intent(inout) :: vstate
    type(rwork_t), intent(inout) :: rwork
    integer,       intent(in   ) :: IORD

    ! Declare local variables
    real(rt) :: ALPH0, ALPH1, HSUM, PROD, T1, XI,XIOLD
    integer    :: I, IBACK, J, JP1, LP1, NQM1, NQM2, NQP1

    !$gpu

    IF ((vstate % NQ .EQ. 2) .AND. (IORD .NE. 1)) RETURN
    NQM1 = vstate % NQ - 1
    NQM2 = vstate % NQ - 2

    select case (vstate % METH)
    case (1)
       go to 100
    case (2)
       go to 200
    end select

    ! -----------------------------------------------------------------------
    !  Nonstiff option...
    !  Check to see if the order is being increased or decreased.
    ! -----------------------------------------------------------------------
100 CONTINUE
    IF (IORD .EQ. 1) GO TO 180
    ! Order decrease. ------------------------------------------------------
    do J = 1, VODE_LMAX
       vstate % EL(J) = ZERO
    end do
    vstate % EL(2) = ONE
    HSUM = ZERO
    do J = 1, NQM2
       ! Construct coefficients of x*(x+xi(1))*...*(x+xi(j)). -----------------
       HSUM = HSUM + vstate % TAU(J)
       XI = HSUM/vstate % HSCAL
       JP1 = J + 1
       do IBACK = 1, JP1
          I = (J + 3) - IBACK
          vstate % EL(I) = vstate % EL(I)*XI + vstate % EL(I-1)
       end do
    end do
    ! Construct coefficients of integrated polynomial. ---------------------
    do J = 2, NQM1
       vstate % EL(J+1) = REAL(vstate % NQ) * vstate % EL(J)/REAL(J)
    end do
    ! Subtract correction terms from YH array. -----------------------------
    do J = 3, vstate % NQ
       do I = 1, VODE_NEQS
          rwork % YH(I,J) = rwork % YH(I,J) - &
               rwork % YH(I,vstate % L) * vstate % EL(J)
       end do
    end do
    RETURN
    ! Order increase. ------------------------------------------------------
    ! Zero out next column in YH array. ------------------------------------
180 CONTINUE
    LP1 = vstate % L + 1
    do I = 1, VODE_NEQS
       rwork % YH(I,LP1) = ZERO
    end do
    RETURN
    ! -----------------------------------------------------------------------
    !  Stiff option...
    !  Check to see if the order is being increased or decreased.
    ! -----------------------------------------------------------------------
200 CONTINUE
    IF (IORD .EQ. 1) GO TO 300
    ! Order decrease. ------------------------------------------------------
    do J = 1, VODE_LMAX
       vstate % EL(J) = ZERO
    end do
    vstate % EL(3) = ONE
    HSUM = ZERO
    do J = 1,NQM2
       ! Construct coefficients of x*x*(x+xi(1))*...*(x+xi(j)). ---------------
       HSUM = HSUM + vstate % TAU(J)
       XI = HSUM/vstate % HSCAL
       JP1 = J + 1
       do IBACK = 1, JP1
          I = (J + 4) - IBACK
          vstate % EL(I) = vstate % EL(I) * XI + vstate % EL(I-1)
       end do
    end do
    ! Subtract correction terms from YH array. -----------------------------
    do J = 3,vstate % NQ
       do I = 1, VODE_NEQS
          rwork % YH(I,J) = rwork % YH(I,J) - &
               rwork % YH(I,vstate % L) * vstate % EL(J)
       end do
    end do
    RETURN
    ! Order increase. ------------------------------------------------------
300 continue
    do J = 1, VODE_LMAX
       vstate % EL(J) = ZERO
    end do
    vstate % EL(3) = ONE
    ALPH0 = -ONE
    ALPH1 = ONE
    PROD = ONE
    XIOLD = ONE
    HSUM = vstate % HSCAL
    IF (vstate % NQ .EQ. 1) GO TO 340
    do J = 1, NQM1
       ! Construct coefficients of x*x*(x+xi(1))*...*(x+xi(j)). ---------------
       JP1 = J + 1
       HSUM = HSUM + vstate % TAU(JP1)
       XI = HSUM/vstate % HSCAL
       PROD = PROD*XI
       ALPH0 = ALPH0 - ONE/REAL(JP1)
       ALPH1 = ALPH1 + ONE/XI
       do IBACK = 1, JP1
          I = (J + 4) - IBACK
          vstate % EL(I) = vstate % EL(I) * XIOLD + vstate % EL(I-1)
       end do
       XIOLD = XI
    end do
340 CONTINUE
    T1 = (-ALPH0 - ALPH1)/PROD
    ! Load column L + 1 in YH array. ---------------------------------------
    LP1 = vstate % L + 1
    do I = 1, VODE_NEQS
       rwork % YH(I,LP1) = T1 * rwork % YH(I,VODE_LMAX)
    end do
    ! Add correction terms to YH array. ------------------------------------
    NQP1 = vstate % NQ + 1
    do J = 3, NQP1
       CALL DAXPYN(VODE_NEQS, vstate % EL(J), &
            rwork % YH(1:VODE_NEQS, LP1), 1, rwork % YH(1:VODE_NEQS, J), 1)
    end do
    RETURN
  end subroutine dvjust

end module cuvode_dvjust_module
