module bs_type_module

  use amrex_fort_module, only : rt => amrex_real
  use burn_type_module, only: neqs, burn_t
  use bs_rpar_indices, only: n_rpar_comps

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  ! BS parameters -- see the discussion in 16.4
  integer, parameter :: KMAXX = 7
  integer :: nseq(KMAXX+1)

  integer, parameter :: bs_neqs = neqs

  type bs_t
     logical :: first
     real(rt) :: eps_old
     real(rt) :: dt_did
     real(rt) :: dt_next
     real(rt) :: a(KMAXX+1)
     real(rt) :: alpha(KMAXX, KMAXX)
     real(rt) :: t_new
     integer :: kmax
     integer :: kopt

     real(rt) :: y(neqs)
     real(rt) :: atol(neqs), rtol(neqs)
     real(rt) :: upar(n_rpar_comps)
     real(rt) :: t, dt, tmax
     integer         :: n

     type(burn_t) :: burn_s

  end type bs_t

  !$acc declare create(nseq)

contains

  subroutine clean_state(state)

    !$acc routine seq

    use amrex_constants_module, only: ONE
    use extern_probin_module, only: SMALL_X_SAFE, renormalize_abundances, MAX_TEMP
    use actual_network, only: nspec, nspec_evolve
    use burn_type_module, only: net_itemp
    use eos_type_module, only : eos_get_small_temp

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    type (bs_t), intent(inout) :: state

    real(rt)  :: small_temp

    ! Ensure that mass fractions always stay positive and sum to 1.
    state % y(1:nspec_evolve) = &
         max(min(state % y(1:nspec_evolve), ONE), SMALL_X_SAFE)

    ! Renormalize abundances as necessary.
    if (renormalize_abundances) then
       call renormalize_species(state)
    endif

    ! Ensure that the temperature always stays within reasonable limits.
    call eos_get_small_temp(small_temp)

    state % y(net_itemp) = min(MAX_TEMP, max(state % y(net_itemp), small_temp))

  end subroutine clean_state


  subroutine renormalize_species(state)

    !$acc routine seq

    use amrex_fort_module, only : rt => amrex_real
    use actual_network, only: nspec, nspec_evolve
    use bs_rpar_indices, only: irp_nspec, n_not_evolved

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    type (bs_t) :: state

    real(rt) :: nspec_sum

    nspec_sum = &
         sum(state % y(1:nspec_evolve)) + &
         sum(state % upar(irp_nspec:irp_nspec+n_not_evolved-1))

    state % y(1:nspec_evolve) = state % y(1:nspec_evolve) / nspec_sum
    state % upar(irp_nspec:irp_nspec+n_not_evolved-1) = &
         state % upar(irp_nspec:irp_nspec+n_not_evolved-1) / nspec_sum

  end subroutine renormalize_species


  subroutine update_thermodynamics(state)

    !$acc routine seq

    use amrex_constants_module, only: ZERO
    use eos_type_module, only: eos_t, eos_input_rt
    use eos_composition_module, only : composition
    use eos_module, only: eos
    use extern_probin_module, only: call_eos_in_rhs, dT_crit
    ! these shouldn't be needed
    use bs_rpar_indices, only: irp_nspec, n_not_evolved
    use actual_network, only : nspec, nspec_evolve

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    type (bs_t) :: state
    type (eos_t) :: eos_state

    ! Several thermodynamic quantities come in via bs % upar -- note: these
    ! are evaluated at the start of the integration, so if things change
    ! dramatically, they will fall out of sync with the current
    ! thermodynamics.

    call bs_to_eos(eos_state, state)

    ! Evaluate the thermodynamics -- if desired. Note that
    ! even if this option is selected, we don't need to do it
    ! for non-self-heating integrations because the temperature
    ! isn't being updated. Also, if it is, we can optionally
    ! set a fraction dT_crit such that we don't call the EOS
    ! if the last temperature we evaluated the EOS at is relatively
    ! close to the current temperature.

    ! Otherwise just do the composition calculations since
    ! that's needed to construct dY/dt. Then make sure
    ! the abundances are safe.

    if (call_eos_in_rhs .and. state % burn_s % self_heat) then

       call eos(eos_input_rt, eos_state)
       call eos_to_bs(eos_state, state)

    else if (abs(eos_state % T - state % burn_s % T_old) > &
         dT_crit * eos_state % T .and. state % burn_s % self_heat) then

       call eos(eos_input_rt, eos_state)

       state % burn_s % dcvdt = (eos_state % cv - state % burn_s % cv) / &
            (eos_state % T - state % burn_s % T_old)
       state % burn_s % dcpdt = (eos_state % cp - state % burn_s % cp) / &
            (eos_state % T - state % burn_s % T_old)
       state % burn_s % T_old  = eos_state % T

       ! note: the update to state % upar(irp_cv) and irp_cp is done
       ! in the call to eos_to_bs that follows 
       call eos_to_bs(eos_state, state)

    else

       call composition(eos_state)

       ! just update what is needed here
       state % burn_s % y_e = eos_state % y_e
       state % burn_s % abar = eos_state % abar
       state % burn_s % zbar = eos_state % zbar

       ! this shouldn't actually be necessary -- we haven't changed X at all,
       ! but roundoff in the multiply / divide change answers slightly.  Leaving
       ! this in for now for the test suite
       state % y(1:nspec_evolve) = eos_state % xn(1:nspec_evolve)
       state % upar(irp_nspec:irp_nspec+n_not_evolved-1) = &
            eos_state % xn(nspec_evolve+1:nspec) 

    endif



  end subroutine update_thermodynamics



  ! Given a bs_t, set up an EOS state.
  ! We could fill the energy component by storing the initial energy in
  ! ts % upar if we wanted, but we won't do that because (1) if we do call the EOS,
  ! it is always in (rho, T) mode and (2) converting back would imply subtracting
  ! off the nuclear energy from the zone's internal energy, which could lead to
  ! issues from roundoff error if the energy released from burning is small.

  subroutine bs_to_eos(state, bs)

    !$acc routine seq

    use actual_network, only: nspec, nspec_evolve
    use eos_type_module, only: eos_t
    use bs_rpar_indices, only: irp_nspec, n_not_evolved
    use burn_type_module, only: net_itemp

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    type (eos_t) :: state
    type (bs_t) :: bs

    state % rho     = bs % burn_s % rho
    state % T       = bs % y(net_itemp)

    state % xn(1:nspec_evolve) = bs % y(1:nspec_evolve)
    state % xn(nspec_evolve+1:nspec) = &
         bs % upar(irp_nspec:irp_nspec+n_not_evolved-1)

    ! we don't copy any of the other quantities, since we can always
    ! access them through the original bs type

  end subroutine bs_to_eos



  ! Given an EOS state, fill the rpar and integration state data.

  subroutine eos_to_bs(state, bs)

    !$acc routine seq

    use network, only: nspec, nspec_evolve
    use eos_type_module, only: eos_t
    use bs_rpar_indices, only: irp_nspec, n_not_evolved
    use burn_type_module, only: net_itemp

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    type (eos_t) :: state
    type (bs_t) :: bs

    bs % burn_s % rho = state % rho

    ! T is funny -- it is both an integration variable and a member of burn_t
    bs % y(net_itemp) = state % T
    bs % burn_s % T = state % T 

    bs % y(1:nspec_evolve) = state % xn(1:nspec_evolve)
    bs % upar(irp_nspec:irp_nspec+n_not_evolved-1) = &
         state % xn(nspec_evolve+1:nspec) 

    bs % burn_s % cp = state % cp
    bs % burn_s % cv = state % cv
    bs % burn_s % abar = state % abar
    bs % burn_s % zbar = state % zbar
    bs % burn_s % eta = state % eta
    bs % burn_s % y_e = state % y_e
    bs % burn_s % cs = state % cs

  end subroutine eos_to_bs



  subroutine burn_to_bs(bs)

    ! Given a burn state, fill the bs_t.  For this implementation, we just
    ! modify the burn_t that is part of the bs_t, in place

    !$acc routine seq

    use network, only: nspec, nspec_evolve
    use bs_rpar_indices, only: irp_nspec, n_not_evolved
    use burn_type_module, only: burn_t, net_itemp, net_ienuc
    use amrex_constants_module, only: ONE

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    type (bs_t) :: bs

    bs % burn_s % rho = bs % burn_s % rho
    bs % y(net_itemp) = bs % burn_s % T

    bs % y(1:nspec_evolve) = bs % burn_s % xn(1:nspec_evolve)
    bs % upar(irp_nspec:irp_nspec+n_not_evolved-1) = &
         bs % burn_s % xn(nspec_evolve+1:nspec) 

    bs % y(net_ienuc) = bs % burn_s % e

    ! we don't need to do anything to thermodynamic quantities, cp, cv, ...

    ! we don't need to do anything to self_heat

  end subroutine burn_to_bs



  subroutine bs_to_burn(bs)
    ! Given a bs_t, set up a burn state.  For this implementation, we just
    ! modify the burn_t that is part of the bs_t, in place
    
    !$acc routine seq

    use actual_network, only: nspec, nspec_evolve
    use bs_rpar_indices, only: irp_nspec, n_not_evolved
    use burn_type_module, only: burn_t, net_itemp, net_ienuc
    use amrex_constants_module, only: ZERO, ONE

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    type (bs_t) :: bs

    bs % burn_s % rho = bs % burn_s % rho
    bs % burn_s % T = bs % y(net_itemp)
    bs % burn_s % e = bs % y(net_ienuc)

    bs % burn_s % xn(1:nspec_evolve) = bs % y(1:nspec_evolve)
    bs % burn_s % xn(nspec_evolve+1:nspec) = &
         bs % upar(irp_nspec:irp_nspec+n_not_evolved-1)

    ! all the other thermodynamic quantities (cp, cv, ...) are already
    ! in the burn_t

    ! self_heat is already set

    bs % burn_s % time = bs % t

  end subroutine bs_to_burn

  subroutine dump_bs_state(bs)

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    type (bs_t) :: bs

    print *, "time: ", bs % t
    print *, "T:    ", bs % burn_s % T 
    print *, "rho:  ", bs % burn_s % rho 
    print *, "X:    ", bs % burn_s % xn(:)

  end subroutine dump_bs_state

end module bs_type_module
