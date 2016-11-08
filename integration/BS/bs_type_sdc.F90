module bs_type_module

  use bl_constants_module, only: HALF, ONE
  use bl_types, only: dp_t
  use sdc_type_module, only: SVAR, SVAR_EVOLVE
  use rpar_indices, only: n_rpar_comps

  implicit none

  ! BS parameters -- see the discussion in 16.4
  integer, parameter :: KMAXX = 7
  integer :: nseq(KMAXX+1)

  integer, parameter :: bs_neqs = SVAR_EVOLVE

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

     real(kind=dp_t) :: y(SVAR_EVOLVE), ydot(SVAR_EVOLVE), jac(SVAR_EVOLVE, SVAR_EVOLVE)
     real(kind=dp_t) :: ydot_a(SVAR_EVOLVE)
     real(kind=dp_t) :: atol(SVAR_EVOLVE), rtol(SVAR_EVOLVE)
     real(kind=dp_t) :: y_init(SVAR_EVOLVE)
     real(kind=dp_t) :: u(n_rpar_comps), u_init(n_rpar_comps), udot_a(n_rpar_comps)
     real(kind=dp_t) :: t, dt, tmax
     integer         :: n

     integer         :: n_rhs, n_jac
     
     integer :: i, j, k
     logical :: T_from_eden
     logical :: self_heat  ! needed only for compatibilty with no-SDC

  end type bs_t

  !$acc declare create(nseq)

contains

  subroutine clean_state(state)

    !$acc routine seq

    use extern_probin_module, only: SMALL_X_SAFE, renormalize_abundances
    use actual_network, only: nspec
    use sdc_type_module, only: SFS, SEDEN, SEINT
    use rpar_indices, only: irp_SRHO, irp_SMX, irp_SMZ
    use eos_module, only: eos_get_small_dens, eos_get_max_dens, eos
    use eos_type_module, only: eos_input_rt, eos_t


    implicit none

    ! this should be larger than any reasonable temperature we will encounter
    real (kind=dp_t), parameter :: MAX_TEMP = 1.0d11

    real (kind=dp_t) :: max_e, ke

    type (bs_t) :: state

    type (eos_t) :: eos_state

    ! Update rho, rho*u, etc.
    call fill_unevolved_variables(state)

    ! Ensure that mass fractions always stay positive and less than one.
    state % y(SFS:SFS+nspec-1) = max(min(state % y(SFS:SFS+nspec-1), state % u(irp_SRHO)), &
                                     state % u(irp_SRHO) * SMALL_X_SAFE)

    ! Renormalize abundances as necessary.
    if (renormalize_abundances) then
       call renormalize_species(state)
    endif

    ! Ensure that internal energy never goes above the maximum limit
    ! provided by the EOS. Same for the internal energy implied by the
    ! total energy (which we get by subtracting kinetic energy).

    eos_state % rho = state % u(irp_SRHO)
    eos_state % T = MAX_TEMP
    eos_state % xn = state % y(SFS:SFS+nspec-1) / state % u(irp_SRHO)

    call eos(eos_input_rt, eos_state)

    max_e = eos_state % e

    state % y(SEINT) = min(state % u(irp_SRHO) * max_e, state % y(SEINT))

    ke = state % y(SEDEN) - HALF * sum(state % u(irp_SMX:irp_SMZ)**2) / state % u(irp_SRHO)

    state % y(SEDEN) = min(state % u(irp_SRHO) * max_e + ke, state % y(SEDEN))

  end subroutine clean_state



  subroutine fill_unevolved_variables(state)

    !$acc routine seq

    use sdc_type_module, only: SRHO, SMX, SMZ
    use rpar_indices, only: irp_SRHO, irp_SMX, irp_SMZ

    implicit none

    type (bs_t) :: state

    ! we are always integrating from t = 0, so there is no offset time
    ! needed here

    state % u(irp_SRHO) = state % u_init(irp_SRHO) + &
         state % udot_a(irp_SRHO) * state % t
    state % u(irp_SMX:irp_SMZ) = state % u_init(irp_SMX:irp_SMZ) + &
         state % udot_a(irp_SMX:irp_SMZ) * state % t

  end subroutine fill_unevolved_variables



  subroutine renormalize_species(state)

    !$acc routine seq

    use sdc_type_module, only: SFS
    use actual_network, only: nspec
    use rpar_indices, only: irp_SRHO

    implicit none

    type (bs_t) :: state

    real(dp_t) :: nspec_sum

    ! Update rho, rho*u, etc.

    call fill_unevolved_variables(state)

    nspec_sum = sum(state % y(SFS:SFS+nspec-1)) / state % u(irp_SRHO)

    state % y(SFS:SFS+nspec-1) = state % y(SFS:SFS+nspec-1) / nspec_sum

  end subroutine renormalize_species



  subroutine sdc_to_bs(sdc, bs)

    !$acc routine seq

    use sdc_type_module, only: sdc_t, SVAR_EVOLVE, SRHO, SMX, SMZ
    use rpar_indices, only: irp_SRHO, irp_SMX, irp_SMZ

    implicit none

    type (sdc_t) :: sdc
    type (bs_t) :: bs

    bs % y = sdc % y(1:SVAR_EVOLVE)
    bs % y_init = bs % y
    bs % ydot_a = sdc % ydot_a(1:SVAR_EVOLVE)

    bs % u(irp_SRHO) = sdc % y(SRHO)
    bs % u(irp_SMX:irp_SMZ) = sdc % y(SMX:SMZ)

    bs % udot_a(irp_SRHO) = sdc % ydot_a(SRHO)
    bs % udot_a(irp_SMX:irp_SMZ) = sdc % ydot_a(SMX:SMZ)

    bs % u_init = bs % u

    bs % i = sdc % i
    bs % j = sdc % j
    bs % k = sdc % k

    bs % n_rhs = sdc % n_rhs
    bs % n_rhs = sdc % n_rhs

    bs % T_from_eden = sdc % T_from_eden

  end subroutine sdc_to_bs



  subroutine bs_to_sdc(sdc, bs)

    !$acc routine seq

    use sdc_type_module, only: sdc_t, SRHO, SMX, SMZ
    use rpar_indices, only: irp_SRHO, irp_SMX, irp_SMZ

    implicit none

    type (sdc_t) :: sdc
    type (bs_t) :: bs

    sdc % y(1:SVAR_EVOLVE) = bs % y

    sdc % y(SRHO) = bs % u(irp_SRHO)
    sdc % y(SMX:SMZ) = bs % u(irp_SMX:irp_SMZ)

    sdc % n_rhs = bs % n_rhs
    sdc % n_jac = bs % n_jac

    sdc % self_heat = bs % self_heat

  end subroutine bs_to_sdc



  subroutine rhs_to_bs(bs, burn)

    !$acc routine seq

    use actual_network, only: nspec_evolve, aion
    use burn_type_module, only: burn_t, net_ienuc
    use sdc_type_module, only: SVAR_EVOLVE, SEDEN, SEINT, SFS
    use rpar_indices, only: irp_SRHO

    implicit none

    type (bs_t) :: bs
    type (burn_t) :: burn

    integer :: n

    call fill_unevolved_variables(bs)

    ! Start with the contribution from the non-reacting sources

    bs % ydot = bs % ydot_a(1:SVAR_EVOLVE)

    ! Add in the reacting terms from the burn_t

    bs % ydot(SFS:SFS+nspec_evolve-1) = bs % ydot(SFS:SFS+nspec_evolve-1) + &
                                        bs % u(irp_SRHO) * burn % ydot(1:nspec_evolve) * aion(1:nspec_evolve)

    bs % ydot(SEINT) = bs % ydot(SEINT) + bs % u(irp_SRHO) * burn % ydot(net_ienuc)
    bs % ydot(SEDEN) = bs % ydot(SEDEN) + bs % u(irp_SRHO) * burn % ydot(net_ienuc)

  end subroutine rhs_to_bs



  subroutine jac_to_bs(bs, burn)

    !$acc routine seq

    use actual_network, only: nspec_evolve, aion
    use burn_type_module, only: burn_t, net_ienuc
    use sdc_type_module, only: SEDEN, SEINT, SFS

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
    bs % jac(SEDEN,SEINT) = burn % jac(net_ienuc,net_ienuc)

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
    use bl_types, only: dp_t
    use eos_type_module, only: eos_input_re, eos_t
    use eos_module, only: eos, eos_get_small_temp, eos_get_max_temp
    use sdc_type_module, only: SEDEN, SEINT, SFS
    use rpar_indices, only: irp_SRHO, irp_SMX, irp_SMZ

    implicit none

    type (bs_t) :: bs
    type (burn_t) :: burn
    type (eos_t) :: eos_state

    real(kind=dp_t) :: rhoInv, min_temp, max_temp

    ! Update rho, rho*u, etc.

    call fill_unevolved_variables(bs)

    rhoInv = ONE / bs % u(irp_SRHO)

    eos_state % rho = bs % u(irp_SRHO)
    eos_state % xn  = bs % y(SFS:SFS+nspec-1) * rhoInv

    if (bs % T_from_eden) then
       eos_state % e = (bs % y(SEDEN) - HALF * rhoInv * sum(bs % u(irp_SMX:irp_SMZ)**2))
    else
       eos_state % e = bs % y(SEINT) * rhoInv
    endif

    ! Give the temperature an initial guess -- use the geometric mean
    ! of the minimum and maximum temperatures.

    call eos_get_small_temp(min_temp)
    call eos_get_max_temp(max_temp)

    eos_state % T = sqrt(min_temp * max_temp)

    call eos(eos_input_re, eos_state)
    call eos_to_burn(eos_state, burn)

    burn % time = bs % t

    burn % n_rhs = bs % n_rhs
    burn % n_jac = bs % n_jac

    burn % self_heat = bs % self_heat

  end subroutine bs_to_burn

  subroutine dump_bs_state(bs)
    
    use eos_type_module, only: eos_input_re, eos_t
    use sdc_type_module, only: SEDEN, SEINT, SFS
    use rpar_indices, only: irp_SRHO, irp_SMX, irp_SMZ
    use actual_network, only: nspec
    use eos_module, only: eos, eos_get_small_temp, eos_get_max_temp

    type (bs_t) :: bs
    type (eos_t) :: eos_state

    real (kind=dp_t) :: rhoInv, min_temp, max_temp

    call fill_unevolved_variables(bs)

    rhoInv = ONE / bs % u(irp_SRHO)

    eos_state % rho = bs % u(irp_SRHO)
    eos_state % xn  = bs % y(SFS:SFS+nspec-1) * rhoInv

    if (bs % T_from_eden) then
       eos_state % e = (bs % y(SEDEN) - HALF * rhoInv * sum(bs % u(irp_SMX:irp_SMZ)**2))
    else
       eos_state % e = bs % y(SEINT) * rhoInv
    endif

    ! Give the temperature an initial guess -- use the geometric mean
    ! of the minimum and maximum temperatures.

    call eos_get_small_temp(min_temp)
    call eos_get_max_temp(max_temp)

    eos_state % T = sqrt(min_temp * max_temp)

    call eos(eos_input_re, eos_state)

    print *, "time: ", bs % t
    print *, "T:    ", eos_state % T
    print *, "rho:  ", eos_state % rho
    print *, "X:    ", eos_state % xn(:)

  end subroutine dump_bs_state

end module bs_type_module
