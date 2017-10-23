module dvode_dvset_module

  use dvode_constants_module

  implicit none

contains

#ifdef CUDA  
  attributes(device) &
#endif  
  subroutine dvset(vstate)

    !$acc routine seq
    
    ! -----------------------------------------------------------------------
    !  Call sequence communication: None
    !  COMMON block variables accessed:
    !      /DVOD01/ -- EL(13), H, TAU(13), TQ(5), L(= NQ + 1),
    !                  METH, NQ, NQWAIT
    ! 
    !  Subroutines called by DVSET: None
    !  Function routines called by DVSET: None
    ! -----------------------------------------------------------------------
    !  DVSET is called by DVSTEP and sets coefficients for use there.
    ! 
    !  For each order NQ, the coefficients in EL are calculated by use of
    !   the generating polynomial lambda(x), with coefficients EL(i).
    !       lambda(x) = EL(1) + EL(2)*x + ... + EL(NQ+1)*(x**NQ).
    !  For the backward differentiation formulas,
    !                                      NQ-1
    !       lambda(x) = (1 + x/xi*(NQ)) * product (1 + x/xi(i) ) .
    !                                      i = 1
    !  For the Adams formulas,
    !                               NQ-1
    !       (d/dx) lambda(x) = c * product (1 + x/xi(i) ) ,
    !                               i = 1
    !       lambda(-1) = 0,    lambda(0) = 1,
    !  where c is a normalization constant.
    !  In both cases, xi(i) is defined by
    !       H*xi(i) = t sub n  -  t sub (n-i)
    !               = H + TAU(1) + TAU(2) + ... TAU(i-1).
    ! 
    ! 
    !  In addition to variables described previously, communication
    !  with DVSET uses the following:
    !    TAU    = A vector of length 13 containing the past NQ values
    !             of H.
    !    EL     = A vector of length 13 in which vset stores the
    !             coefficients for the corrector formula.
    !    TQ     = A vector of length 5 in which vset stores constants
    !             used for the convergence test, the error test, and the
    !             selection of H at a new order.
    !    METH   = The basic method indicator.
    !    NQ     = The current order.
    !    L      = NQ + 1, the length of the vector stored in EL, and
    !             the number of columns of the YH array being used.
    !    NQWAIT = A counter controlling the frequency of order changes.
    !             An order change is about to be considered if NQWAIT = 1.
    ! -----------------------------------------------------------------------
    !

    use vode_type_module, only: rwork_t
    use dvode_type_module, only: dvode_t
    use bl_types, only: dp_t

    implicit none
  
    type(dvode_t), intent(inout) :: vstate
    real(dp_t) :: EM(13)
    real(dp_t) :: AHATN0, ALPH0, CNQM1, CSUM, ELP
    real(dp_t) :: EM0, FLOTI, FLOTL, FLOTNQ, HSUM, RXI, RXIS, S
    real(dp_t) :: T1, T2, T3, T4, T5, T6, XI
    integer    :: I, IBACK, J, JP1, NQM1, NQM2

    ! Parameter declaration
    real(dp_t), parameter :: CORTES = 0.1D0

    FLOTL = REAL(vstate % L)
    NQM1 = vstate % NQ - 1
    NQM2 = vstate % NQ - 2
    GO TO (100, 200), vstate % METH

    !  Set coefficients for Adams methods. ----------------------------------
100 IF (vstate % NQ .NE. 1) GO TO 110
    vstate % EL(1) = ONE
    vstate % EL(2) = ONE
    vstate % TQ(1) = ONE
    vstate % TQ(2) = TWO
    vstate % TQ(3) = SIX*vstate % TQ(2)
    vstate % TQ(5) = ONE
    GO TO 300
110 continue
    HSUM = vstate % H
    EM(1) = ONE
    FLOTNQ = FLOTL - ONE
    do I = 2, vstate % L
       EM(I) = ZERO
    end do
    do J = 1, NQM1
       IF ((J .NE. NQM1) .OR. (vstate % NQWAIT .NE. 1)) GO TO 130
       S = ONE
       CSUM = ZERO
       do I = 1, NQM1
          CSUM = CSUM + S*EM(I)/REAL(I+1)
          S = -S
       end do
       vstate % TQ(1) = EM(NQM1)/(FLOTNQ*CSUM)
130    continue
       RXI = vstate % H/HSUM
       do IBACK = 1, J
          I = (J + 2) - IBACK
          EM(I) = EM(I) + EM(I-1)*RXI
       end do
       HSUM = HSUM + vstate % TAU(J)
    end do
    ! Compute integral from -1 to 0 of polynomial and of x times it. -------
    S = ONE
    EM0 = ZERO
    CSUM = ZERO
    do I = 1, vstate % NQ
       FLOTI = REAL(I)
       EM0 = EM0 + S*EM(I)/FLOTI
       CSUM = CSUM + S*EM(I)/(FLOTI+ONE)
       S = -S
    end do
    ! In EL, form coefficients of normalized integrated polynomial. --------
    S = ONE/EM0
    vstate % EL(1) = ONE
    do I = 1, vstate % NQ
       vstate % EL(I+1) = S*EM(I)/REAL(I)
    end do
    XI = HSUM/vstate % H
    vstate % TQ(2) = XI*EM0/CSUM
    vstate % TQ(5) = XI/vstate % EL(vstate % L)
    IF (vstate % NQWAIT .NE. 1) GO TO 300
    ! For higher order control constant, multiply polynomial by 1+x/xi(q). -
    RXI = ONE/XI
    do IBACK = 1, vstate % NQ
       I = (vstate % L + 1) - IBACK
       EM(I) = EM(I) + EM(I-1)*RXI
    end do
    ! Compute integral of polynomial. --------------------------------------
    S = ONE
    CSUM = ZERO
    do I = 1, vstate % L
       CSUM = CSUM + S*EM(I)/REAL(I+1)
       S = -S
    end do
    vstate % TQ(3) = FLOTL*EM0/CSUM
    GO TO 300

    ! Set coefficients for BDF methods. ------------------------------------
200 continue
    do I = 3, vstate % L
       vstate % EL(I) = ZERO
    end do
    vstate % EL(1) = ONE
    vstate % EL(2) = ONE
    ALPH0 = -ONE
    AHATN0 = -ONE
    HSUM = vstate % H
    RXI = ONE
    RXIS = ONE
    IF (vstate % NQ .EQ. 1) GO TO 240
    do J = 1, NQM2
       ! In EL, construct coefficients of (1+x/xi(1))*...*(1+x/xi(j+1)). ------
       HSUM = HSUM + vstate % TAU(J)
       RXI = vstate % H/HSUM
       JP1 = J + 1
       ALPH0 = ALPH0 - ONE/REAL(JP1)
       do IBACK = 1, JP1
          I = (J + 3) - IBACK
          vstate % EL(I) = vstate % EL(I) + vstate % EL(I-1)*RXI
       end do
    end do
    ALPH0 = ALPH0 - ONE/REAL(vstate % NQ)
    RXIS = -vstate % EL(2) - ALPH0
    HSUM = HSUM + vstate % TAU(NQM1)
    RXI = vstate % H/HSUM
    AHATN0 = -vstate % EL(2) - RXI
    do IBACK = 1, vstate % NQ
       I = (vstate % NQ + 2) - IBACK
       vstate % EL(I) = vstate % EL(I) + vstate % EL(I-1)*RXIS
    end do
240 continue
    T1 = ONE - AHATN0 + ALPH0
    T2 = ONE + REAL(vstate % NQ)*T1
    vstate % TQ(2) = ABS(ALPH0*T2/T1)
    vstate % TQ(5) = ABS(T2/(vstate % EL(vstate % L)*RXI/RXIS))
    IF (vstate % NQWAIT .NE. 1) GO TO 300
    CNQM1 = RXIS/vstate % EL(vstate % L)
    T3 = ALPH0 + ONE/REAL(vstate % NQ)
    T4 = AHATN0 + RXI
    ELP = T3/(ONE - T4 + T3)
    vstate % TQ(1) = ABS(ELP/CNQM1)
    HSUM = HSUM + vstate % TAU(vstate % NQ)
    RXI = vstate % H/HSUM
    T5 = ALPH0 - ONE/REAL(vstate % NQ+1)
    T6 = AHATN0 - RXI
    ELP = T2/(ONE - T6 + T5)
    vstate % TQ(3) = ABS(ELP*RXI*(FLOTL + ONE)*T5)
300 continue
    vstate % TQ(4) = CORTES*vstate % TQ(2)
    RETURN
  end subroutine dvset

end module dvode_dvset_module
