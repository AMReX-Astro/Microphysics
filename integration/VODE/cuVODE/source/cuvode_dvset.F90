module cuvode_dvset_module

  use cuvode_types_module, only: dvode_t
  use amrex_fort_module, only: rt => amrex_real

  implicit none

contains

#if defined(AMREX_USE_CUDA) && !defined(AMREX_USE_GPU_PRAGMA)
  attributes(device) &
#endif
  subroutine dvset(vstate)

    !$acc routine seq
    
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
    !    TAU    = A array of length 13 containing the past NQ values
    !             of H.
    !    EL     = A array of length 13 in which vset stores the
    !             coefficients for the corrector formula.
    !    TQ     = A array of length 5 in which vset stores constants
    !             used for the convergence test, the error test, and the
    !             selection of H at a new order.
    !    NQ     = The current order.
    !    L      = NQ + 1, the length of the array stored in EL, and
    !             the number of columns of the YH array being used.
    !    NQWAIT = A counter controlling the frequency of order changes.
    !             An order change is about to be considered if NQWAIT = 1.
    ! -----------------------------------------------------------------------

    implicit none

    ! Declare arguments
    type(dvode_t), intent(inout) :: vstate

    ! Declare local variables
    real(rt) :: EM(13)
    real(rt) :: AHATN0, ALPH0, CNQM1, CSUM, ELP
    real(rt) :: EM0, FLOTI, FLOTL, FLOTNQ, HSUM, RXI, RXIS, S
    real(rt) :: T1, T2, T3, T4, T5, T6, XI
    integer    :: I, IBACK, J, JP1, NQM1, NQM2

    ! Parameter declaration
    real(rt), parameter :: CORTES = 0.1e0_rt

    !$gpu

    FLOTL = REAL(vstate % L)
    NQM1 = vstate % NQ - 1
    NQM2 = vstate % NQ - 2

    ! Set coefficients for BDF methods. ------------------------------------
    do I = 3, vstate % L
       vstate % EL(I) = 0.0_rt
    end do
    vstate % EL(1) = 1.0_rt
    vstate % EL(2) = 1.0_rt
    ALPH0 = -1.0_rt
    AHATN0 = -1.0_rt
    HSUM = vstate % H
    RXI = 1.0_rt
    RXIS = 1.0_rt
    IF (vstate % NQ /= 1) then
       do J = 1, NQM2
          ! In EL, construct coefficients of (1+x/xi(1))*...*(1+x/xi(j+1)). ------
          HSUM = HSUM + vstate % TAU(J)
          RXI = vstate % H/HSUM
          JP1 = J + 1
          ALPH0 = ALPH0 - 1.0_rt/REAL(JP1)
          do IBACK = 1, JP1
             I = (J + 3) - IBACK
             vstate % EL(I) = vstate % EL(I) + vstate % EL(I-1)*RXI
          end do
       end do
       ALPH0 = ALPH0 - 1.0_rt/REAL(vstate % NQ)
       RXIS = -vstate % EL(2) - ALPH0
       HSUM = HSUM + vstate % TAU(NQM1)
       RXI = vstate % H/HSUM
       AHATN0 = -vstate % EL(2) - RXI
       do IBACK = 1, vstate % NQ
          I = (vstate % NQ + 2) - IBACK
          vstate % EL(I) = vstate % EL(I) + vstate % EL(I-1)*RXIS
       end do
    end IF

    T1 = 1.0_rt - AHATN0 + ALPH0
    T2 = 1.0_rt + REAL(vstate % NQ)*T1
    vstate % TQ(2) = ABS(ALPH0*T2/T1)
    vstate % TQ(5) = ABS(T2/(vstate % EL(vstate % L)*RXI/RXIS))

    IF (vstate % NQWAIT == 1) then
       CNQM1 = RXIS/vstate % EL(vstate % L)
       T3 = ALPH0 + 1.0_rt/REAL(vstate % NQ)
       T4 = AHATN0 + RXI
       ELP = T3/(1.0_rt - T4 + T3)
       vstate % TQ(1) = ABS(ELP/CNQM1)
       HSUM = HSUM + vstate % TAU(vstate % NQ)
       RXI = vstate % H/HSUM
       T5 = ALPH0 - 1.0_rt/REAL(vstate % NQ+1)
       T6 = AHATN0 - RXI
       ELP = T2/(1.0_rt - T6 + T5)
       vstate % TQ(3) = ABS(ELP*RXI*(FLOTL + 1.0_rt)*T5)
    end IF
    vstate % TQ(4) = CORTES*vstate % TQ(2)
    RETURN
  end subroutine dvset

end module cuvode_dvset_module
