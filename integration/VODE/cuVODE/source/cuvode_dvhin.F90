module cuvode_dvhin_module

  use cuvode_parameters_module, only: VODE_LMAX, VODE_NEQS, VODE_MAXORD
  use cuvode_types_module, only: dvode_t
  use amrex_fort_module, only: rt => amrex_real

  implicit none

contains

#if defined(AMREX_USE_CUDA) && !defined(AMREX_USE_GPU_PRAGMA)
  attributes(device) &
#endif
  subroutine dvhin(vstate, H0, NITER, IER)
  
    !$acc routine seq
    
    ! -----------------------------------------------------------------------
    !  This routine computes the step size, H0, to be attempted on the
    !  first step, when the user has not supplied a value for this.
    ! 
    !  First we check that TOUT - T0 differs significantly from zero.  Then
    !  an iteration is done to approximate the initial second derivative
    !  and this is used to define h from w.r.m.s.norm(h**2 * yddot / 2) = 1.
    !  A bias factor of 1/2 is applied to the resulting h.
    !  The sign of H0 is inferred from the initial values of TOUT and T0.
    ! -----------------------------------------------------------------------
  
#ifdef TRUE_SDC
    use sdc_vode_rhs_module, only : f_rhs
#else
    use vode_rhs_module, only: f_rhs
#endif

    implicit none

    ! Declare arguments
    type(dvode_t), intent(inout) :: vstate
    real(rt),      intent(inout) :: H0
    integer,       intent(  out) :: NITER, IER

    ! Declare local variables
    real(rt) :: AFI, ATOLI, DELYI, H, HG, HLB, HNEW, HRAT
    real(rt) :: HUB, T1, TDIST, TROUND, YDDNRM
    integer  :: I, ITER

    real(rt), parameter :: PT1 = 0.1e0_rt

    logical :: do_iterations

    !$gpu

    NITER = 0
    TDIST = abs(vstate % TOUT - vstate % T)
    TROUND = vstate % UROUND * max(abs(vstate % T), abs(vstate % TOUT))

    if (TDIST < 2.0_rt*TROUND) then
       ! Error return for vstate % TOUT - vstate % T too small. --------------------------------
       IER = -1
       return
    end if

    ! Set a lower bound on h based on the roundoff level in vstate % T and vstate % TOUT. ---
    HLB = 100.0_rt * TROUND
    ! Set an upper bound on h based on vstate % TOUT-vstate % T and the initial Y and YDOT. -
    HUB = PT1*TDIST
    ATOLI = vstate % ATOL(1)

    do I = 1, VODE_NEQS
       ATOLI = vstate % ATOL(I)
       DELYI = PT1 * abs(vstate % YH(I,1)) + ATOLI
       AFI = abs(vstate % YH(I,2))
       if (AFI*HUB > DELYI) HUB = DELYI/AFI
    end do

    ! Set initial guess for h as geometric mean of upper and lower bounds. -
    ITER = 0
    HG = sqrt(HLB*HUB)

    ! If the bounds have crossed, exit with the mean value. ----------------
    do_iterations = .true.
    if (HUB < HLB) then
       H0 = HG
       do_iterations = .false.
    end if

    if (do_iterations) then

       ! Looping point for iteration. -----------------------------------------
       do while (.true.)

          ! Estimate the second derivative as a difference quotient in f. --------
          H = sign(HG, vstate % TOUT - vstate % T)
          T1 = vstate % T + H
          do I = 1, VODE_NEQS
             vstate % Y(I) = vstate % YH(I,1) + H*vstate % YH(I,2)
          end do

          call f_rhs(T1, vstate, vstate % ACOR)

          do I = 1, VODE_NEQS
             vstate % ACOR(I) = (vstate % ACOR(I) - vstate % YH(I,2))/H
          end do
          YDDNRM = sqrt(sum((vstate % ACOR * vstate % EWT)**2) / VODE_NEQS)

          ! Get the corresponding new value of h. --------------------------------
          if (YDDNRM*HUB*HUB > 2.0_rt) then
             HNEW = sqrt(2.0_rt/YDDNRM)
          else
             HNEW = sqrt(HG*HUB)
          end if
          ITER = ITER + 1

          ! -----------------------------------------------------------------------
          !  Test the stopping conditions.
          !  Stop if the new and previous h values differ by a factor of .lt. 2.
          !  Stop if four iterations have been done.  Also, stop with previous h
          !  if HNEW/HG .gt. 2 after first iteration, as this probably means that
          !  the second derivative value is bad because of cancellation error.
          ! -----------------------------------------------------------------------
          if (iter >= 4) exit

          HRAT = HNEW/HG
          if ( (HRAT > 0.5_rt) .and. (HRAT < 2.0_rt) ) exit
          if ( (ITER >= 2) .and. (HNEW > 2.0_rt*HG) ) then
             HNEW = HG
             exit
          end if
          HG = HNEW
       end do

       ! Iteration done.  Apply bounds, bias factor, and sign.  Then exit. ----
       H0 = HNEW*0.5_rt
       if (H0 < HLB) H0 = HLB
       if (H0 > HUB) H0 = HUB

    end if

    H0 = sign(H0, vstate % TOUT - vstate % T)
    NITER = ITER
    IER = 0

  end subroutine dvhin

end module cuvode_dvhin_module
