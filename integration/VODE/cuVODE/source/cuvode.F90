module cuvode_module

  use cuvode_parameters_module, only: VODE_LMAX, VODE_NEQS, VODE_LIW,   &
                                      VODE_LENWM, VODE_MAXORD
  use cuvode_types_module, only: dvode_t, rwork_t
  use vode_rpar_indices
  use amrex_fort_module, only: rt => amrex_real
  use linpack_module
#ifdef AMREX_USE_CUDA
  use cudafor
#endif

  use cuvode_constants_module
  use cuvode_dvhin_module
  use cuvode_dvindy_module
  use cuvode_dvstep_module
  
  implicit none

  public :: dvode
  
contains

#if defined(AMREX_USE_CUDA) && !defined(AMREX_USE_GPU_PRAGMA)
  attributes(device) &
#endif
  subroutine dvode(vstate, rwork, IWORK, MF)

    !$acc routine seq

#ifndef AMREX_USE_CUDA
    use cuvode_output_module, only: xerrwd
#endif
#ifdef TRUE_SDC
    use sdc_vode_rhs_module, only: f_rhs, jac
#else
    use vode_rhs_module, only: f_rhs, jac
#endif
    use cuvode_dvnorm_module, only: dvnorm ! function

    implicit none

    ! Declare arguments
    type(dvode_t), intent(inout) :: vstate
    type(rwork_t), intent(inout) :: rwork
    integer,       intent(inout) :: IWORK(VODE_LIW)
    integer,       intent(in   ) :: MF

    ! Declare local variables
    logical    :: IHIT
    real(rt) :: ATOLI, EWTI, H0, HMAX, HMX
    real(rt) :: RH, RTOLI, SIZE, TCRIT, TNEXT, TOLSF, TP
    integer    :: I, IER, IFLAG, KGO, LENJ, LENP
    integer    :: MBAND, MFA, ML, MU, NITER
    integer    :: NSLAST
    integer    :: pivot(VODE_NEQS)
#ifndef AMREX_USE_CUDA
    character (len=80) :: MSG
#endif

    ! Parameter declarations
    integer, parameter :: MXSTP0 = 500
    integer, parameter :: MXHNL0 = 10
    real(rt), parameter :: PT2 = 0.2e0_rt

    logical :: skip_loop_start

    !$gpu

    vstate % INIT = 0
    IF (vstate % TOUT .EQ. vstate % T) RETURN

    ! Compute Jacobian choices based on MF.

    vstate % JSV = SIGN(1,MF)
    MFA = ABS(MF)
    vstate % METH = MFA/10
    vstate % MITER = MFA - 10*vstate % METH

    ! Set some defaults that the user can override.

    if (IWORK(6) <= 0) then
       vstate % MXSTEP = MXSTP0
    else
       vstate % MXSTEP = IWORK(6)
    end if

    if (IWORK(7) <= 0) then
       vstate % MXHNIL = MXHNL0
    else
       vstate % MXHNIL = IWORK(7)
    end if

    H0 = ZERO

    vstate % HMIN = ZERO
    vstate % HMXI = ZERO

    ! -----------------------------------------------------------------------
    !  Arrays stored in RWORK are denoted  CONDOPT, YH, WM, EWT, SAVF, ACOR.
    !  Within WM, LOCJS is the location of the saved Jacobian (JSV .gt. 0).
    ! -----------------------------------------------------------------------

    IF (vstate % MITER .EQ. 1 .OR. vstate % MITER .EQ. 2) THEN
       vstate % LOCJS = VODE_NEQS*VODE_NEQS + 3
    ENDIF

    ! Check RTOL and ATOL for legality. ------------------------------------
    RTOLI = vstate % RTOL(1)
    ATOLI = vstate % ATOL(1)
    do I = 1,VODE_NEQS

       RTOLI = vstate % RTOL(I)
       ATOLI = vstate % ATOL(I)

       if (RTOLI .LT. ZERO) then
#ifndef AMREX_USE_CUDA
          MSG = 'DVODE--  RTOL(I1) is R1 .lt. 0.0_rt        '
          CALL XERRWD (MSG, 40, 19, 1, 1, I, 0, 1, RTOLI, ZERO)
#endif
          vstate % ISTATE = -3
          return
       end if

       if (ATOLI .LT. ZERO) then
#ifndef AMREX_USE_CUDA
          MSG = 'DVODE--  ATOL(I1) is R1 .lt. 0.0_rt        '
          CALL XERRWD (MSG, 40, 20, 1, 1, I, 0, 1, ATOLI, ZERO)
#endif
          vstate % ISTATE = -3
          return
       end if

    end do

    ! -----------------------------------------------------------------------
    !  All remaining initializations, the initial call to F,
    !  and the calculation of the initial step size.
    !  The error weights in EWT are inverted after being loaded.
    ! -----------------------------------------------------------------------

    vstate % UROUND = epsilon(1.0_rt)
    vstate % TN = vstate % T

    vstate % JSTART = 0
    IF (vstate % MITER .GT. 0) RWORK % wm(1) = SQRT(vstate % UROUND)
    vstate % CCMXJ = PT2
    vstate % MSBJ = 50
    vstate % NHNIL = 0
    vstate % NST = 0
    vstate % NJE = 0
    vstate % NNI = 0
    vstate % NCFN = 0
    vstate % NETF = 0
    vstate % NLU = 0
    vstate % NSLJ = 0
    NSLAST = 0
    vstate % HU = ZERO
    vstate % NQU = 0

    ! Initial call to F.  -------------------------

    CALL f_rhs (vstate % T, vstate, rwork % yh(:,2))
    vstate % NFE = 1
    ! Load the initial value array in YH. ---------------------------------
    rwork % YH(1:VODE_NEQS,1) = vstate % Y(1:VODE_NEQS)

    ! Load and invert the EWT array.  (H is temporarily set to 1.0.) -------
    vstate % NQ = 1
    vstate % H = ONE

    do I = 1,VODE_NEQS
       rwork % EWT(I) = vstate % RTOL(I) * abs(rwork % YH(I,1)) + vstate % ATOL(I)
       if (rwork % ewt(I) .LE. ZERO) then
#ifndef AMREX_USE_CUDA
          EWTI = rwork % ewt(I)
          MSG = 'DVODE--  EWT(I1) is R1 .le. 0.0_rt         '
          CALL XERRWD (MSG, 40, 21, 1, 1, I, 0, 1, EWTI, ZERO)
#endif
          vstate % ISTATE = -3
          return
       end if

       rwork % ewt(I) = ONE/rwork % ewt(I)
    end do

    ! Call DVHIN to set initial step size H0 to be attempted. --------------
    CALL DVHIN (vstate, rwork, H0, NITER, IER)
    vstate % NFE = vstate % NFE + NITER

    if (IER .NE. 0) then
#ifndef AMREX_USE_CUDA
       MSG='DVODE--  TOUT (=R1) too close to T(=R2) to start integration'
       CALL XERRWD (MSG, 60, 22, 1, 0, 0, 0, 2, vstate % TOUT, vstate % T)
#endif
       vstate % ISTATE = -3
       return
    end if

    ! Load H with H0 and scale YH(*,2) by H0. ------------------------------
    vstate % H = H0
    rwork % YH(:,2) = rwork % YH(:,2) * H0

    skip_loop_start = .true.

    ! -----------------------------------------------------------------------
    !  Block E.
    !  The next block is normally executed for all calls and contains
    !  the call to the one-step core integrator DVSTEP.
    !
    !  This is a looping point for the integration steps.
    !
    !  First check for too many steps being taken, update EWT (if not at
    !  start of problem), check for too much accuracy being requested, and
    !  check for H below the roundoff level in T.
    ! -----------------------------------------------------------------------

    do while (.true.)

       if (.not. skip_loop_start) then

          IF ((vstate % NST-NSLAST) .GE. vstate % MXSTEP) then
             ! The maximum number of steps was taken before reaching TOUT. ----------
#ifndef AMREX_USE_CUDA
             MSG = 'DVODE--  At current T (=R1), MXSTEP (=I1) steps   '
             CALL XERRWD (MSG, 50, 201, 1, 0, 0, 0, 0, ZERO, ZERO)
             MSG = '      taken on this call before reaching TOUT     '
             CALL XERRWD (MSG, 50, 201, 1, 1, vstate % MXSTEP, 0, 1, vstate % TN, ZERO)
#endif
             vstate % ISTATE = -1

             vstate % Y(1:VODE_NEQS) = rwork % YH(1:VODE_NEQS,1)

             vstate % T = vstate % TN
             IWORK(11) = vstate % NST
             IWORK(12) = vstate % NFE
             IWORK(13) = vstate % NJE
             IWORK(14) = vstate % NQU
             IWORK(15) = vstate % NQ
             IWORK(19) = vstate % NLU
             IWORK(20) = vstate % NNI
             IWORK(21) = vstate % NCFN
             IWORK(22) = vstate % NETF
             return

          end IF

          do I = 1,VODE_NEQS

             rwork % EWT(I) = vstate % RTOL(I) * abs(rwork % YH(I,1)) + vstate % ATOL(I)

             IF (rwork % ewt(I) .LE. ZERO) then
                ! EWT(i) .le. 0.0 for some i (not at start of problem). ----------------
#ifndef AMREX_USE_CUDA
                EWTI = rwork % ewt(I)
                MSG = 'DVODE--  At T (=R1), EWT(I1) has become R2 .le. 0.'
                CALL XERRWD (MSG, 50, 202, 1, 1, I, 0, 2, vstate % TN, EWTI)
#endif
                vstate % ISTATE = -6

                vstate % Y(1:VODE_NEQS) = rwork % YH(1:VODE_NEQS,1)

                vstate % T = vstate % TN
                IWORK(11) = vstate % NST
                IWORK(12) = vstate % NFE
                IWORK(13) = vstate % NJE
                IWORK(14) = vstate % NQU
                IWORK(15) = vstate % NQ
                IWORK(19) = vstate % NLU
                IWORK(20) = vstate % NNI
                IWORK(21) = vstate % NCFN
                IWORK(22) = vstate % NETF

                return
             end IF
             rwork % ewt(I) = ONE/rwork % ewt(I)
          end do

       else
          skip_loop_start = .false.
       end if

       TOLSF = vstate % UROUND * DVNORM (rwork % YH(:,1), rwork % EWT)
       IF (TOLSF > ONE) then
          TOLSF = TOLSF*TWO

          if (vstate % NST .EQ. 0) then
#ifndef AMREX_USE_CUDA
             MSG = 'DVODE--  At start of problem, too much accuracy   '
             CALL XERRWD (MSG, 50, 26, 1, 0, 0, 0, 0, ZERO, ZERO)
             MSG='      requested for precision of machine:   see TOLSF (=R1) '
             CALL XERRWD (MSG, 60, 26, 1, 0, 0, 0, 1, TOLSF, ZERO)
#endif
             vstate % ISTATE = -3
             return
          end if

          ! Too much accuracy requested for machine precision. -------------------
#ifndef AMREX_USE_CUDA
          MSG = 'DVODE--  At T (=R1), too much accuracy requested  '
          CALL XERRWD (MSG, 50, 203, 1, 0, 0, 0, 0, ZERO, ZERO)
          MSG = '      for precision of machine:   see TOLSF (=R2) '
          CALL XERRWD (MSG, 50, 203, 1, 0, 0, 0, 2, vstate % TN, TOLSF)
#endif
          vstate % ISTATE = -2

          vstate % Y(1:VODE_NEQS) = rwork % YH(1:VODE_NEQS,1)

          vstate % T = vstate % TN
          IWORK(11) = vstate % NST
          IWORK(12) = vstate % NFE
          IWORK(13) = vstate % NJE
          IWORK(14) = vstate % NQU
          IWORK(15) = vstate % NQ
          IWORK(19) = vstate % NLU
          IWORK(20) = vstate % NNI
          IWORK(21) = vstate % NCFN
          IWORK(22) = vstate % NETF

          return

       end IF

       IF ((vstate % TN + vstate % H) == vstate % TN) then
          vstate % NHNIL = vstate % NHNIL + 1
          IF (vstate % NHNIL <= vstate % MXHNIL) then
#ifndef AMREX_USE_CUDA
             MSG = 'DVODE--  Warning: internal T (=R1) and H (=R2) are'
             CALL XERRWD (MSG, 50, 101, 1, 0, 0, 0, 0, ZERO, ZERO)
             MSG='      such that in the machine, T + H = T on the next step  '
             CALL XERRWD (MSG, 60, 101, 1, 0, 0, 0, 0, ZERO, ZERO)
             MSG = '      (H = step size). solver will continue anyway'
             CALL XERRWD (MSG, 50, 101, 1, 0, 0, 0, 2, vstate % TN, vstate % H)
#endif
             IF (vstate % NHNIL == vstate % MXHNIL) then
#ifndef AMREX_USE_CUDA
                MSG = 'DVODE--  Above warning has been issued I1 times.  '
                CALL XERRWD (MSG, 50, 102, 1, 0, 0, 0, 0, ZERO, ZERO)
                MSG = '      it will not be issued again for this problem'
                CALL XERRWD (MSG, 50, 102, 1, 1, vstate % MXHNIL, 0, 0, ZERO, ZERO)
#endif
             end IF
          end IF
       end IF

       CALL DVSTEP(pivot, rwork, vstate)

       ! Branch on KFLAG.  Note: In this version, KFLAG can not be set to -3.
       !  KFLAG .eq. 0,   -1,  -2

       if (vstate % kflag == -1) then
          ! KFLAG = -1.  Error test failed repeatedly or with ABS(H) = HMIN. -----
#ifndef AMREX_USE_CUDA
          MSG = 'DVODE--  At T(=R1) and step size H(=R2), the error'
          CALL XERRWD (MSG, 50, 204, 1, 0, 0, 0, 0, ZERO, ZERO)
          MSG = '      test failed repeatedly or with abs(H) = HMIN'
          CALL XERRWD (MSG, 50, 204, 1, 0, 0, 0, 2, vstate % TN, vstate % H)
#endif
          vstate % ISTATE = -4

          ! Set Y array, T, and optional output. --------------------------------
          vstate % Y(1:VODE_NEQS) = rwork % YH(1:VODE_NEQS,1)

          vstate % T = vstate % TN
          IWORK(11) = vstate % NST
          IWORK(12) = vstate % NFE
          IWORK(13) = vstate % NJE
          IWORK(14) = vstate % NQU
          IWORK(15) = vstate % NQ
          IWORK(19) = vstate % NLU
          IWORK(20) = vstate % NNI
          IWORK(21) = vstate % NCFN
          IWORK(22) = vstate % NETF

          return

       else if (vstate % kflag == -2) then
          ! KFLAG = -2.  Convergence failed repeatedly or with ABS(H) = HMIN. ----
#ifndef AMREX_USE_CUDA
          MSG = 'DVODE--  At T (=R1) and step size H (=R2), the    '
          CALL XERRWD (MSG, 50, 205, 1, 0, 0, 0, 0, ZERO, ZERO)
          MSG = '      corrector convergence failed repeatedly     '
          CALL XERRWD (MSG, 50, 205, 1, 0, 0, 0, 0, ZERO, ZERO)
          MSG = '      or with abs(H) = HMIN   '
          CALL XERRWD (MSG, 30, 205, 1, 0, 0, 0, 2, vstate % TN, vstate % H)
#endif
          vstate % ISTATE = -5

          ! Set Y array, T, and optional output. --------------------------------
          vstate % Y(1:VODE_NEQS) = rwork % YH(1:VODE_NEQS,1)

          vstate % T = vstate % TN
          IWORK(11) = vstate % NST
          IWORK(12) = vstate % NFE
          IWORK(13) = vstate % NJE
          IWORK(14) = vstate % NQU
          IWORK(15) = vstate % NQ
          IWORK(19) = vstate % NLU
          IWORK(20) = vstate % NNI
          IWORK(21) = vstate % NCFN
          IWORK(22) = vstate % NETF

          return
       end if

       ! -----------------------------------------------------------------------
       !  Block F.
       !  The following block handles the case of a successful return from the
       !  core integrator (KFLAG = 0).  Test for stop conditions.
       ! -----------------------------------------------------------------------

       vstate % INIT = 1
       vstate % KUTH = 0

       ! If TOUT has been reached, interpolate. -------------------
       IF ((vstate % TN - vstate % TOUT) * vstate % H .LT. ZERO) cycle
       CALL DVINDY (vstate, rwork, IFLAG)
       vstate % T = vstate % TOUT

       vstate % ISTATE = 2
       IWORK(11) = vstate % NST
       IWORK(12) = vstate % NFE
       IWORK(13) = vstate % NJE
       IWORK(14) = vstate % NQU
       IWORK(15) = vstate % NEWQ
       IWORK(19) = vstate % NLU
       IWORK(20) = vstate % NNI
       IWORK(21) = vstate % NCFN
       IWORK(22) = vstate % NETF

       return

    end do

    ! -----------------------------------------------------------------------
    !  Block G.
    !  The following block handles all successful returns from DVODE.
    !  vstate % ISTATE is set to 2, and the optional output is loaded into the work
    !  arrays before returning.
    ! -----------------------------------------------------------------------

    vstate % Y(1:VODE_NEQS) = rwork % YH(1:VODE_NEQS,1)

    vstate % T = vstate % TN

    vstate % ISTATE = 2
    IWORK(11) = vstate % NST
    IWORK(12) = vstate % NFE
    IWORK(13) = vstate % NJE
    IWORK(14) = vstate % NQU
    IWORK(15) = vstate % NEWQ
    IWORK(19) = vstate % NLU
    IWORK(20) = vstate % NNI
    IWORK(21) = vstate % NCFN
    IWORK(22) = vstate % NETF

    return

  end subroutine dvode
      
end module cuvode_module
