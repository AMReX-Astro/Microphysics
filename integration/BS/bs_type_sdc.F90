module bs_type_module

  use bl_types, only: dp_t
  use sdc_type_module, only: SVAR

  implicit none

  ! BS parameters -- see the discussion in 16.4
  integer, parameter :: KMAXX = 7
  integer :: nseq(KMAXX+1)

  integer, parameter :: bs_neqs = SVAR

  ! error codes
  integer, parameter :: IERR_NONE = 0
  integer, parameter :: IERR_DT_TOO_SMALL = -100
  integer, parameter :: IERR_TOO_MANY_STEPS = -101
  integer, parameter :: IERR_DT_UNDERFLOW = -102
  integer, parameter :: IERR_NO_CONVERGENCE = -103

  integer, parameter :: IERR_LU_DECOMPOSITION_ERROR = -200

  real(kind=dp_t), parameter :: S1 = 0.25_dp_t
  real(kind=dp_t), parameter :: S2 = 0.7_dp_t

  real(kind=dp_t), parameter :: RED_BIG_FACTOR = 0.7_dp_t
  real(kind=dp_t), parameter :: RED_SMALL_FACTOR = 1.e-5_dp_t
  real(kind=dp_t), parameter :: SCALMX = 0.1_dp_t

  type bs_t

     logical :: first
     real(kind=dp_t) :: eps_old
     real(kind=dp_t) :: dt_did
     real(kind=dp_t) :: dt_next
     real(kind=dp_t) :: a(KMAXX+1)
     real(kind=dp_t) :: alpha(KMAXX, KMAXX)
     real(kind=dp_t) :: t_new
     integer :: kmax
     integer :: kopt

     real(kind=dp_t) :: y(SVAR), ydot(SVAR), jac(SVAR, SVAR)
     real(kind=dp_t) :: ydot_a(SVAR)
     real(kind=dp_t) :: atol(SVAR), rtol(SVAR)
     real(kind=dp_t) :: t, dt, tmax
     integer         :: n
     integer         :: n_rhs, n_jac

     integer :: i, j, k
     logical :: T_from_eden

  end type bs_t

  !$acc declare create(nseq)

contains

  subroutine clean_state(state)

    !$acc routine seq

    use actual_network, only: nspec
    use sdc_type_module, only: SRHO, SFS

    implicit none

    ! this is the absolute cutoff for species -- note that this might
    ! be larger than small_x that the user set, but the issue is that
    ! we can have underflow issues if the integrator has to keep track
    ! of species mass fractions much smaller than this.
    real (kind=dp_t), parameter :: SMALL_X_SAFE = 1.0d-30

    type (bs_t) :: state

    ! Ensure that mass fractions always stay positive and less than one.

    state % y(SFS:SFS+nspec-1) = max(min(state % y(SFS:SFS+nspec-1), state % y(SRHO)), &
                                     state % y(SRHO) * SMALL_X_SAFE)

  end subroutine clean_state



  subroutine renormalize_species(state)

    !$acc routine seq

    use sdc_type_module, only: SRHO, SFS
    use actual_network, only: nspec

    implicit none

    type (bs_t) :: state

    real(dp_t) :: nspec_sum

    nspec_sum = sum(state % y(SFS:SFS+nspec-1)) / state % y(SRHO)

    state % y(SFS:SFS+nspec-1) = state % y(SFS:SFS+nspec-1) / nspec_sum

  end subroutine renormalize_species



  subroutine sdc_to_bs(sdc, bs)

    !$acc routine seq

    use sdc_type_module, only: sdc_t

    implicit none

    type (sdc_t) :: sdc
    type (bs_t) :: bs

    bs % y = sdc % y

    bs % ydot_a = sdc % ydot_a

    bs % i = sdc % i
    bs % j = sdc % j
    bs % k = sdc % k

    bs % T_from_eden = sdc % T_from_eden

  end subroutine sdc_to_bs



  subroutine bs_to_sdc(sdc, bs)

    !$acc routine seq

    use sdc_type_module, only: sdc_t

    implicit none

    type (sdc_t) :: sdc
    type (bs_t) :: bs

    sdc % y = bs % y

  end subroutine bs_to_sdc



  subroutine rhs_to_bs(bs, burn)

    !$acc routine seq

    use actual_network, only: nspec_evolve, aion
    use burn_type_module, only: burn_t, net_ienuc
    use sdc_type_module, only: SRHO, SEDEN, SEINT, SFS

    implicit none

    type (bs_t) :: bs
    type (burn_t) :: burn

    integer :: n

    ! Start with the contribution from the non-reacting sources

    bs % ydot = bs % ydot_a

    ! Add in the reacting terms from the burn_t

    bs % ydot(SFS:SFS+nspec_evolve-1) = bs % ydot(SFS:SFS+nspec_evolve-1) + &
                                        bs % y(SRHO) * burn % ydot(1:nspec_evolve) * aion(1:nspec_evolve)

    bs % ydot(SEINT) = bs % ydot(SEINT) + bs % y(SRHO) * burn % ydot(net_ienuc)
    bs % ydot(SEDEN) = bs % ydot(SEDEN) + bs % y(SRHO) * burn % ydot(net_ienuc)

  end subroutine rhs_to_bs



  subroutine jac_to_bs(bs, burn)

    !$acc routine seq

    use actual_network, only: nspec_evolve, aion
    use burn_type_module, only: burn_t, net_ienuc
    use sdc_type_module, only: SRHO, SEDEN, SEINT, SFS

    implicit none

    type (bs_t) :: bs
    type (burn_t) :: burn

    integer :: n

    ! Copy the Jacobian from the burn

    bs % jac(SFS:SFS+nspec_evolve-1,SFS:SFS+nspec_evolve-1) = burn % jac(1:nspec_evolve,1:nspec_evolve)
    bs % jac(SFS:SFS+nspec_evolve-1,SEDEN) = burn % jac(1:nspec_evolve,net_ienuc)
    bs % jac(SFS:SFS+nspec_evolve-1,SEINT) = burn % jac(1:nspec_evolve,net_ienuc)

    bs % jac(SEDEN,SFS:SFS+nspec_evolve-1) = burn % jac(net_ienuc,1:nspec_evolve)
    bs % jac(SEDEN,SEDEN) = burn % jac(net_ienuc,net_ienuc)
    bs % jac(SEDEN,SEINT) = burn % Jac(net_ienuc,net_ienuc)

    bs % jac(SEINT,SFS:SFS+nspec_evolve-1) = burn % jac(net_ienuc,1:nspec_evolve)
    bs % jac(SEINT,SEDEN) = burn % jac(net_ienuc,net_ienuc)
    bs % jac(SEINT,SEINT) = burn % jac(net_ienuc,net_ienuc)

    ! Scale it to match our variables. We don't need to worry about the rho
    ! dependence, since every one of the SDC variables is linear in rho, so
    ! we just need to focus on the Y --> X conversion.

    do n = 1, nspec_evolve
       bs % jac(SFS+n-1,:) = bs % jac(SFS+n-1,:) * aion(n)
       bs % jac(:,SFS+n-1) = bs % jac(:,SFS+n-1) / aion(n)
    enddo

  end subroutine jac_to_bs



  subroutine bs_to_burn(bs, burn)
    
    !$acc routine seq

    use actual_network, only: nspec
    use burn_type_module, only: burn_t, burn_to_eos, eos_to_burn
    use bl_constants_module, only: HALF, ONE
    use bl_types, only: dp_t
    use eos_type_module, only: eos_input_re, eos_t
    use eos_module, only: eos
    use sdc_type_module, only: SRHO, SMX, SMZ, SEDEN, SEINT, SFS

    implicit none

    type (bs_t) :: bs
    type (burn_t) :: burn
    type (eos_t) :: eos_state

    real(kind=dp_t) :: e, rhoInv

    rhoInv = ONE / bs % y(SRHO)

    burn % rho = bs % y(SRHO)
    burn % xn  = bs % y(SFS:SFS+nspec-1) * rhoInv

    if (bs % T_from_eden) then

       e = (bs % y(SEDEN) - HALF * rhoInv * sum(bs % y(SMX:SMZ)**2))

    else

       e = bs % y(SEINT) * rhoInv

    endif

    ! If the temperature is smaller than the EOS can handle, allow it to
    ! reset the temperature accordingly.

    eos_state % reset = .true.

    ! Do not check the validity of inputs going into the EOS call.
    ! Sometimes we may stray into a meaningless state and we want
    ! to be able to get through the EOS call without a failure so
    ! that we can return to a meaningful state in the convergence.

    eos_state % check_inputs = .false.

    call burn_to_eos(burn, eos_state)
    call eos(eos_input_re, eos_state)
    call eos_to_burn(eos_state, burn)

    burn % time = bs % t

  end subroutine bs_to_burn

end module bs_type_module
