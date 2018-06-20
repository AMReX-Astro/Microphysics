module dvode_dvindy_module

  use vode_type_module, only: rwork_t
  use vode_parameters_module, only: VODE_LMAX, VODE_NEQS, VODE_LIW,   &
                                    VODE_LENWM, VODE_MAXORD, VODE_ITOL
  use dvode_type_module, only: dvode_t
  use amrex_fort_module, only: rt => amrex_real
#ifndef CUDA  
  use dvode_output_module, only: xerrwd
#endif
  use blas_module

  use dvode_constants_module

  implicit none

contains

  AMREX_DEVICE subroutine dvindy(vstate, rwork, IFLAG)
  
    !$acc routine seq
    
    ! -----------------------------------------------------------------------
    !  Call sequence input -- T, K, YH, LDYH
    !  Call sequence output -- DKY, IFLAG
    !  COMMON block variables accessed:
    !      /DVOD01/ --  H, TN, UROUND, L, N, NQ
    !      /DVOD02/ --  HU
    ! 
    !  Subroutines called by DVINDY: DSCAL, XERRWD
    !  Function routines called by DVINDY: None
    ! -----------------------------------------------------------------------
    !  DVINDY computes interpolated values of the K-th derivative of the
    !  dependent variable vector y, and stores it in DKY.  This routine
    !  is called within the package with K = 0 and T = TOUT, but may
    !  also be called by the user for any K up to the current order.
    !  (See detailed instructions in the usage documentation.)
    ! -----------------------------------------------------------------------
    !  The computed values in DKY are gotten by interpolation using the
    !  Nordsieck history array YH.  This array corresponds uniquely to a
    !  vector-valued polynomial of degree NQCUR or less, and DKY is set
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
    !  yh => rwork % yh
    !  t => vstate % tout
    !  dky => vstate % y

    implicit none

    ! Declare arguments
    type(dvode_t), intent(inout) :: vstate
    type(rwork_t), intent(in   ) :: rwork
    integer,       intent(  out) :: IFLAG

    ! Declare local variables
    real(rt) :: C, R, S, TFUZZ, TN1, TP
    integer    :: I, IC, J, JB, JB2, JJ, JJ1, JP1
#ifndef CUDA
    character (len=80) :: MSG
#endif

    integer, parameter :: K = ZERO
    !don -- remove logic for K != 0 (we only use K = 0).
    
    IFLAG = 0
    IF (K .LT. 0 .OR. K .GT. vstate % NQ) GO TO 80
    TFUZZ = HUN * vstate % UROUND * (vstate % TN + vstate % HU)
    TP = vstate % TN - vstate % HU - TFUZZ
    TN1 = vstate % TN + TFUZZ
    IF ((vstate % TOUT-TP)*(vstate % TOUT-TN1) .GT. ZERO) GO TO 90

    S = (vstate % TOUT - vstate % TN)/vstate % H
    IC = 1
    IF (K .EQ. 0) GO TO 15
    JJ1 = vstate % L - K
    do JJ = JJ1, vstate % NQ
       IC = IC*JJ
    end do
15  continue
    C = REAL(IC)
    do I = 1, VODE_NEQS
       vstate % Y(I) = C * rwork % YH(I,vstate % L)
    end do
    IF (K .EQ. vstate % NQ) GO TO 55
    JB2 = vstate % NQ - K
    do JB = 1, JB2
       J = vstate % NQ - JB
       JP1 = J + 1
       IC = 1
       IF (K .EQ. 0) GO TO 35
       JJ1 = JP1 - K
       do JJ = JJ1, J
          IC = IC*JJ
       end do
35     continue
       C = REAL(IC)
       do I = 1, VODE_NEQS
          vstate % Y(I) = C * rwork % YH(I,JP1) + S*vstate % Y(I)
       end do
    end do
    IF (K .EQ. 0) RETURN
55  continue
    R = vstate % H**(-K)

    CALL DSCALN (VODE_NEQS, R, vstate % Y, 1)
    RETURN

80  continue
#ifndef CUDA    
    MSG = 'DVINDY-- K (=I1) illegal      '
    CALL XERRWD (MSG, 30, 51, 1, 1, K, 0, 0, ZERO, ZERO)
#endif    
    IFLAG = -1
    RETURN
90  continue
#ifndef CUDA    
    MSG = 'DVINDY-- vstate % TOUT (=R1) illegal      '
    CALL XERRWD (MSG, 30, 52, 1, 0, 0, 0, 1, vstate % TOUT, ZERO)
    MSG='      vstate % TOUT not in interval TCUR - HU (= R1) to TCUR (=R2)      '
    CALL XERRWD (MSG, 60, 52, 1, 0, 0, 0, 2, TP, vstate % TN)
#endif    
     IFLAG = -2
    RETURN
  end subroutine dvindy

end module dvode_dvindy_module
