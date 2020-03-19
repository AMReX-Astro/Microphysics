module cuvode_dvindy_module

  use cuvode_parameters_module, only: VODE_LMAX, VODE_NEQS, VODE_MAXORD
  use cuvode_types_module, only: dvode_t
  use amrex_fort_module, only: rt => amrex_real

  implicit none

contains

#if defined(AMREX_USE_CUDA) && !defined(AMREX_USE_GPU_PRAGMA)
  attributes(device) &
#endif
  subroutine dvindy(vstate, IFLAG)

    !$acc routine seq

    ! -----------------------------------------------------------------------
    !  DVINDY computes interpolated values of the K-th derivative of the
    !  dependent variable array y, and stores it in DKY.  This routine
    !  is called within the package with K = 0 and T = TOUT, but may
    !  also be called by the user for any K up to the current order.
    !  (See detailed instructions in the usage documentation.)
    ! -----------------------------------------------------------------------
    !  The computed values in DKY are gotten by interpolation using the
    !  Nordsieck history array YH.  This array corresponds uniquely to a
    !  array-valued polynomial of degree NQCUR or less, and DKY is set
    !  to the K-th derivative of this polynomial at T.
    !  The formula for DKY is:
    !               q
    !   DKY(i)  =  sum  c(j,K) * (T - TN)**(j-K) * H**(-j) * YH(i,j+1)
    !              j=K
    !  where  c(j,K) = j*(j-1)*...*(j-K+1), q = NQCUR, TN = TCUR, H = HCUR.
    !  The quantities  NQ = NQCUR, L = NQ+1, N, TN, and H are
    !  communicated by COMMON.  The above sum is done in reverse order.
    !  IFLAG is returned negative if either K or T is out of bounds.
    !
    !  Discussion above and comments in driver explain all variables.
    ! -----------------------------------------------------------------------
    !
    ! Note: the following variable replacements have been made ---
    !  yh => vstate % yh
    !  t => vstate % tout
    !  dky => vstate % y

    implicit none

    ! Declare arguments
    type(dvode_t), intent(inout) :: vstate
    integer,       intent(  out) :: IFLAG

    ! Declare local variables
    real(rt) :: C, R, S, TFUZZ, TN1, TP
    integer    :: I, IC, J, JB, JB2, JJ, JJ1, JP1

    !$gpu

    IFLAG = 0

    TFUZZ = 100.0_rt * vstate % UROUND * (vstate % TN + vstate % HU)
    TP = vstate % TN - vstate % HU - TFUZZ
    TN1 = vstate % TN + TFUZZ

    S = (vstate % TOUT - vstate % TN)/vstate % H
    IC = 1
    C = REAL(IC)
    do I = 1, VODE_NEQS
       vstate % Y(I) = C * vstate % YH(I,vstate % L)
    end do

    JB2 = vstate % NQ
    do JB = 1, JB2
       J = vstate % NQ - JB
       JP1 = J + 1
       IC = 1
       C = REAL(IC)
       do I = 1, VODE_NEQS
          vstate % Y(I) = C * vstate % YH(I,JP1) + S*vstate % Y(I)
       end do
    end do

  end subroutine dvindy

end module cuvode_dvindy_module
