module sdc_type_module

  use amrex_fort_module, only : rt => amrex_real
  use burn_type_module, only: neqs, burn_t
  use sdc_rpar_indices, only: n_rpar_comps
  use sdc_quadrature_module, only : SDC_NODES
  implicit none

  integer, parameter :: sdc_neqs = neqs

  type sdc_t
     logical :: first
     real(rt) :: eps_old
     real(rt) :: dt_did
     real(rt) :: dt_next
     real(rt) :: t_new

     ! storage for the solution
     real(rt) :: y(neqs)

     ! storage for the solution on each time node
     real(rt) :: y_node(0:SDC_NODES-1, neqs)

     ! storage for the previous iteration's RHS
     real(rt) :: f_old(0:SDC_NODES-1, neqs)

     real(rt) :: atol(neqs), rtol(neqs)
     real(rt) :: upar(n_rpar_comps)
     real(rt) :: t, dt, tmax
     integer :: n

     type(burn_t) :: burn_s

  end type sdc_t

  !$acc declare create(nseq)

contains

  subroutine clean_state(state)

    !$acc routine seq

    use amrex_constants_module, only: ONE
    use extern_probin_module, only: SMALL_X_SAFE, renormalize_abundances, MAX_TEMP
    use actual_network, only: nspec, nspec_evolve
    use burn_type_module, only: net_itemp
    use eos_type_module, only : eos_get_small_temp

    implicit none

    type (sdc_t), intent(inout) :: state

    real (rt) :: small_temp

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
    use sdc_rpar_indices, only: irp_nspec, n_not_evolved

    implicit none

    type (sdc_t) :: state

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
    use eos_type_module, only: eos_t, eos_input_rt, composition
    use eos_module, only: eos
    use extern_probin_module, only: call_eos_in_rhs, dT_crit
    ! these shouldn't be needed
    use sdc_rpar_indices, only: irp_nspec, n_not_evolved
    use actual_network, only : nspec, nspec_evolve

    implicit none

    type (sdc_t) :: state
    type (eos_t) :: eos_state

    ! Several thermodynamic quantities come in via bs % upar -- note: these
    ! are evaluated at the start of the integration, so if things change
    ! dramatically, they will fall out of sync with the current
    ! thermodynamics.

    call sdc_to_eos(eos_state, state)

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



  ! Given a sdc_t, set up an EOS state.
  ! We could fill the energy component by storing the initial energy in
  ! ts % upar if we wanted, but we won't do that because (1) if we do call the EOS,
  ! it is always in (rho, T) mode and (2) converting back would imply subtracting
  ! off the nuclear energy from the zone's internal energy, which could lead to
  ! issues from roundoff error if the energy released from burning is small.

  subroutine sdc_to_eos(state, sdc)

    !$acc routine seq

    use actual_network, only: nspec, nspec_evolve
    use eos_type_module, only: eos_t
    use bs_rpar_indices, only: irp_nspec, n_not_evolved
    use burn_type_module, only: net_itemp

    implicit none

    type (eos_t) :: state
    type (sdc_t) :: sdc

    state % rho     = sdc % burn_s % rho
    state % T       = sdc % y(net_itemp)

    state % xn(1:nspec_evolve) = sdc % y(1:nspec_evolve)
    state % xn(nspec_evolve+1:nspec) = &
         sdc % upar(irp_nspec:irp_nspec+n_not_evolved-1)

    ! we don't copy any of the other quantities, since we can always
    ! access them through the original sdc type

  end subroutine sdc_to_eos



  ! Given an EOS state, fill the rpar and integration state data.

  subroutine eos_to_sdc(state, sdc)

    !$acc routine seq

    use network, only: nspec, nspec_evolve
    use eos_type_module, only: eos_t
    use sdc_rpar_indices, only: irp_nspec, n_not_evolved
    use burn_type_module, only: net_itemp

    implicit none

    type (eos_t) :: state
    type (sdc_t) :: sdc

    sdc % burn_s % rho = state % rho

    ! T is funny -- it is both an integration variable and a member of burn_t
    sdc % y(net_itemp) = state % T
    sdc % burn_s % T = state % T 

    sdc % y(1:nspec_evolve) = state % xn(1:nspec_evolve)
    sdc % upar(irp_nspec:irp_nspec+n_not_evolved-1) = &
         state % xn(nspec_evolve+1:nspec) 

    sdc % burn_s % cp = state % cp
    sdc % burn_s % cv = state % cv
    sdc % burn_s % abar = state % abar
    sdc % burn_s % zbar = state % zbar
    sdc % burn_s % eta = state % eta
    sdc % burn_s % y_e = state % y_e
    sdc % burn_s % cs = state % cs

  end subroutine eos_to_sdc


  subroutine burn_to_sdc(sdc)

    ! Given a burn state, fill the sdc_t.  For this implementation, we just
    ! modify the burn_t that is part of the sdc_t, in place

    !$acc routine seq

    use network, only: nspec, nspec_evolve
    use sdc_rpar_indices, only: irp_nspec, n_not_evolved
    use burn_type_module, only: burn_t, net_itemp, net_ienuc
    use amrex_constants_module, only: ONE

    implicit none

    type (sdc_t) :: sdc

    sdc % burn_s % rho = sdc % burn_s % rho
    sdc % y(net_itemp) = sdc % burn_s % T

    sdc % y(1:nspec_evolve) = sdc % burn_s % xn(1:nspec_evolve)
    sdc % upar(irp_nspec:irp_nspec+n_not_evolved-1) = &
         sdc % burn_s % xn(nspec_evolve+1:nspec) 

    sdc % y(net_ienuc) = sdc % burn_s % e

    ! we don't need to do anything to thermodynamic quantities, cp, cv, ...

    ! we don't need to do anything to self_heat

  end subroutine burn_to_sdc


  subroutine sdc_to_burn(sdc)
    ! Given a sdc_t, set up a burn state.  For this implementation, we just
    ! modify the burn_t that is part of the sdc_t, in place
    
    !$acc routine seq

    use actual_network, only: nspec, nspec_evolve
    use sdc_rpar_indices, only: irp_nspec, n_not_evolved
    use burn_type_module, only: burn_t, net_itemp, net_ienuc
    use amrex_constants_module, only: ZERO, ONE

    implicit none

    type (sdc_t) :: sdc

    sdc % burn_s % rho = sdc % burn_s % rho
    sdc % burn_s % T = sdc % y(net_itemp)
    sdc % burn_s % e = sdc % y(net_ienuc)

    sdc % burn_s % xn(1:nspec_evolve) = sdc % y(1:nspec_evolve)
    sdc % burn_s % xn(nspec_evolve+1:nspec) = &
         sdc % upar(irp_nspec:irp_nspec+n_not_evolved-1)

    ! all the other thermodynamic quantities (cp, cv, ...) are already
    ! in the burn_t

    ! self_heat is already set

    sdc % burn_s % time = sdc % t

  end subroutine sdc_to_burn

  subroutine dump_sdc_state(sdc)

    implicit none

    type (sdc_t) :: sdc

    print *, "time: ", sdc % t
    print *, "T:    ", sdc % burn_s % T 
    print *, "rho:  ", sdc % burn_s % rho 
    print *, "X:    ", sdc % burn_s % xn(:)

  end subroutine dump_sdc_state

end module sdc_type_module
