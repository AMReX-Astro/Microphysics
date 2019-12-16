module vode_type_module

  use amrex_fort_module, only: rt => amrex_real
  use amrex_constants_module
  use cuvode_parameters_module, only : VODE_NEQS

  use network, only : nspec, nspec_evolve, aion, aion_inv

  use vode_rpar_indices
  use sdc_type_module

  use extern_probin_module, only: renormalize_abundances

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  ! this should be larger than any reasonable temperature we will encounter   
  real (kind=rt), parameter :: MAX_TEMP = 1.0e11_rt          

  public

contains

  subroutine clean_state(time, y, rpar)

    use eos_type_module, only : eos_t, eos_input_re, eos_input_rt
    use eos_module, only : eos

    use amrex_fort_module, only : rt => amrex_real
    real(rt), intent(in) :: time
    real(rt) :: y(SVAR_EVOLVE), rpar(n_rpar_comps)

    real(rt) :: max_e, ke

    type (eos_t) :: eos_state

    !$gpu

    ! update rho, rho*u, etc.
    call fill_unevolved_variables(time, y, rpar)

    ! Ensure that mass fractions always stay positive.
    y(SFS:SFS+nspec-1) = &
         max(min(y(SFS:SFS+nspec-1), rpar(irp_SRHO)), &
             rpar(irp_SRHO) * 1.e-200_rt)

    ! renormalize abundances as necessary
    if (renormalize_abundances) then
       call renormalize_species(time, y, rpar)
    endif

#ifdef SDC_EVOLVE_ENERGY

    ! Ensure that internal energy never goes above the maximum limit
    ! provided by the EOS. Same for the internal energy implied by the
    ! total energy (which we get by subtracting kinetic energy).
    eos_state % rho = rpar(irp_SRHO)
    eos_state % T = MAX_TEMP
    eos_state % xn = y(SFS:SFS+nspec-1) / rpar(irp_SRHO)

    call eos(eos_input_rt, eos_state)

    max_e = eos_state % e

    y(SEINT) = min(rpar(irp_SRHO) * max_e, y(SEINT))

    ke = y(SEDEN) - HALF * sum(rpar(irp_SMX:irp_SMZ)**2) / rpar(irp_SRHO)

    y(SEDEN) = min(rpar(irp_SRHO) * max_e + ke, y(SEDEN))

#endif

  end subroutine clean_state


  subroutine fill_unevolved_variables(time, y, rpar)

    use amrex_fort_module, only : rt => amrex_real
    real(rt), intent(in) :: time
    real(rt) :: y(SVAR_EVOLVE), rpar(n_rpar_comps)

    !$gpu

#if defined(SDC_EVOLVE_ENERGY)

    ! we are always integrating from t = 0, so there is no offset
    ! time needed here.  The indexing of irp_ydot_a is based on
    ! the indices in sdc_type_module
    rpar(irp_SRHO) = rpar(irp_u_init-1+irp_SRHO) + &
         rpar(irp_ydot_a-1+SRHO) * time

    rpar(irp_SMX) = rpar(irp_u_init-1+irp_SMX) + rpar(irp_ydot_a-1+SMX) * time
    rpar(irp_SMY) = rpar(irp_u_init-1+irp_SMY) + rpar(irp_ydot_a-1+SMY) * time
    rpar(irp_SMZ) = rpar(irp_u_init-1+irp_SMZ) + rpar(irp_ydot_a-1+SMZ) * time

#elif defined(SDC_EVOLVE_ENTHALPY)

    ! Keep density consistent with the partial densities.
    rpar(irp_SRHO) = sum(y(SFS:SFS - 1 + nspec))

#endif

  end subroutine fill_unevolved_variables

  subroutine renormalize_species(time, y, rpar)

    use amrex_fort_module, only : rt => amrex_real
    real(rt), intent(in) :: time
    real(rt) :: y(SVAR_EVOLVE), rpar(n_rpar_comps)

    real(rt) :: nspec_sum

    !$gpu

    ! We only renormalize species when evolving energy because
    ! when we evolve enthalpy, we define the density as
    ! the sum of the partial densities rho*X for each species.
#ifdef SDC_EVOLVE_ENERGY

    ! update rho, rho*u, etc.
    call fill_unevolved_variables(time, y, rpar)

    nspec_sum = sum(y(SFS:SFS-1+nspec)) / rpar(irp_SRHO)

    y(SFS:SFS-1+nspec) = y(SFS:SFS-1+nspec) / nspec_sum

#endif

  end subroutine renormalize_species


  ! Given a burn state, fill the rpar and integration state data.

  subroutine sdc_to_vode(sdc, y, rpar)

    use amrex_fort_module, only : rt => amrex_real
    type (sdc_t) :: sdc
    real(rt)   :: rpar(n_rpar_comps)
    real(rt)   :: y(SVAR_EVOLVE)

    !$gpu

    y(:) = sdc % y(1:SVAR_EVOLVE)

    ! advective sources
    rpar(irp_ydot_a:irp_ydot_a-1+SVAR) = sdc % ydot_a(:)

#if defined(SDC_EVOLVE_ENERGY)

    ! unevolved state variables
    rpar(irp_SRHO) = sdc % y(SRHO)
    rpar(irp_SMX:irp_SMZ) = sdc % y(SMX:SMZ)

    ! initial state for unevolved variables
    rpar(irp_u_init-1+irp_SRHO) = sdc % y(SRHO)
    rpar(irp_u_init-1+irp_SMX:irp_u_init-1+irp_SMZ) = sdc % y(SMX:SMZ)

    ! other parameters
    if (sdc % T_from_eden) then
       rpar(irp_T_from_eden) = ONE
    else
       rpar(irp_T_from_eden) = -ONE
    endif

#elif defined(SDC_EVOLVE_ENTHALPY)

    rpar(irp_p0)   = sdc % p0
    rpar(irp_SRHO) = sdc % rho

#endif

#ifdef NONAKA_PLOT

    ! bookkeeping information
    rpar(irp_i) = sdc % i
    rpar(irp_j) = sdc % j
    rpar(irp_k) = sdc % k
    rpar(irp_iter) = sdc % sdc_iter

#endif

  end subroutine sdc_to_vode

  subroutine vode_to_sdc(time, y, rpar, sdc)

    use amrex_fort_module, only : rt => amrex_real
    real(rt), intent(in) :: time
    type (sdc_t) :: sdc
    real(rt)    :: rpar(n_rpar_comps)
    real(rt)    :: y(SVAR_EVOLVE)

    !$gpu

    sdc % y(1:SVAR_EVOLVE) = y(:)

    ! unevolved state variables
    call fill_unevolved_variables(time, y, rpar)

#if defined(SDC_EVOLVE_ENERGY)

    sdc % y(SRHO) = rpar(irp_SRHO)
    sdc % y(SMX:SMZ) = rpar(irp_SMX:irp_SMZ)

#elif defined(SDC_EVOLVE_ENTHALPY)

    sdc % p0  = rpar(irp_p0)
    sdc % rho = rpar(irp_SRHO)

#endif

#ifdef NONAKA_PLOT

    sdc % i = rpar(irp_i)
    sdc % j = rpar(irp_j)
    sdc % k = rpar(irp_k)
    sdc % sdc_iter = rpar(irp_iter)

#endif

  end subroutine vode_to_sdc


  subroutine rhs_to_vode(time, burn_state, ydot_react, y, ydot, rpar)

    use burn_type_module, only : burn_t, net_ienuc, neqs

    use amrex_fort_module, only : rt => amrex_real
    real(rt), intent(in) :: time
    real(rt), intent(in) :: rpar(n_rpar_comps)
    real(rt), intent(in) :: y(SVAR_EVOLVE)
    type(burn_t), intent(in) :: burn_state
    real(rt), intent(in) :: ydot_react(neqs)
    real(rt), intent(out) :: ydot(SVAR_EVOLVE)

    !$gpu

    call fill_unevolved_variables(time, y, rpar)

    ! ydot_react has just the contribution to the RHS from the
    ! reaction network.  Note that these are in terms of dY/dt

    ! start with the contribution from the non-reacting sources
    ydot(:) = rpar(irp_ydot_a:irp_ydot_a-1+SVAR_EVOLVE)

    ! add in the reacting terms -- here we convert from dY/dt to dX/dt
    ydot(SFS:SFS-1+nspec_evolve) = ydot(SFS:SFS-1+nspec_evolve) + &
         rpar(irp_SRHO) * aion(1:nspec_evolve) * ydot_react(1:nspec_evolve)

#if defined(SDC_EVOLVE_ENERGY)

    ydot(SEINT) = ydot(SEINT) + rpar(irp_SRHO) * ydot_react(net_ienuc)
    ydot(SEDEN) = ydot(SEDEN) + rpar(irp_SRHO) * ydot_react(net_ienuc)

#elif defined(SDC_EVOLVE_ENTHALPY)

    ydot(SENTH) = ydot(SENTH) + rpar(irp_SRHO) * ydot_react(net_ienuc)

#endif

  end subroutine rhs_to_vode


  subroutine jac_to_vode(time, jac_react, y, jac, rpar)

    ! this is only used with an analytic Jacobian

    use burn_type_module, only : net_ienuc, neqs

    use amrex_fort_module, only : rt => amrex_real
    real(rt), intent(in) :: time
    real(rt)    :: rpar(n_rpar_comps)
    real(rt)    :: y(SVAR_EVOLVE)
    real(rt), intent(in) :: jac_react(neqs, neqs)
    real(rt)    :: jac(SVAR_EVOLVE,SVAR_EVOLVE)

    integer :: n

    !$gpu

#if defined(SDC_EVOLVE_ENERGY)

    jac(SFS:SFS+nspec_evolve-1,SFS:SFS+nspec_evolve-1) = jac_react(1:nspec_evolve,1:nspec_evolve)
    jac(SFS:SFS+nspec_evolve-1,SEDEN) = jac_react(1:nspec_evolve,net_ienuc)
    jac(SFS:SFS+nspec_evolve-1,SEINT) = jac_react(1:nspec_evolve,net_ienuc)

    jac(SEDEN,SFS:SFS+nspec_evolve-1) = jac_react(net_ienuc,1:nspec_evolve)
    jac(SEDEN,SEDEN) = jac_react(net_ienuc,net_ienuc)
    jac(SEDEN,SEINT) = jac_react(net_ienuc,net_ienuc)

    jac(SEINT,SFS:SFS+nspec_evolve-1) = jac_react(net_ienuc,1:nspec_evolve)
    jac(SEINT,SEDEN) = jac_react(net_ienuc,net_ienuc)
    jac(SEINT,SEINT) = jac_react(net_ienuc,net_ienuc)

#elif defined(SDC_EVOLVE_ENTHALPY)

    jac(SFS:SFS+nspec_evolve-1,SFS:SFS+nspec_evolve-1) = jac_react(1:nspec_evolve,1:nspec_evolve)
    jac(SFS:SFS+nspec_evolve-1,SENTH) = jac_react(1:nspec_evolve,net_ienuc)

    jac(SENTH,SFS:SFS+nspec_evolve-1) = jac_react(net_ienuc,1:nspec_evolve)
    jac(SENTH,SENTH) = jac_react(net_ienuc,net_ienuc)

#endif

    ! Scale it to match our variables. We don't need to worry about
    ! the rho dependence, since every one of the SDC variables is
    ! linear in rho, so we just need to focus on the Y --> X
    ! conversion.
    do n = 1, nspec
       jac(SFS+n-1,:) = jac(SFS+n-1,:) * aion(n)
       jac(:,SFS+n-1) = jac(:,SFS+n-1) * aion_inv(n)
    enddo

  end subroutine jac_to_vode


  subroutine vode_to_burn(time, y, rpar, burn_state)

    use eos_type_module, only : eos_t, eos_input_re, eos_input_rt, eos_input_rp, eos_input_rh
    use eos_type_module, only : eos_get_small_temp, eos_get_max_temp
    use eos_module, only : eos
    use burn_type_module, only : eos_to_burn, burn_t

#ifdef SDC_EVOLVE_ENTHALPY
    use meth_params_module, only: use_tfromp
#endif

    use amrex_fort_module, only : rt => amrex_real
    type (burn_t) :: burn_state
    real(rt), intent(in) :: time
    real(rt)    :: rpar(n_rpar_comps)
    real(rt)    :: y(SVAR_EVOLVE)

    type(eos_t) :: eos_state

    real(rt) :: rhoInv, min_temp, max_temp                               

    !$gpu

    ! update rho, rho*u, etc.
    call fill_unevolved_variables(time, y, rpar)

    rhoInv = ONE / rpar(irp_SRHO)

    eos_state % rho = rpar(irp_SRHO)
    eos_state % xn  = y(SFS:SFS+nspec-1) * rhoInv

#if defined(SDC_EVOLVE_ENERGY)

    if (rpar(irp_T_from_eden) > ZERO) then
       eos_state % e = (y(SEDEN) - HALF*rhoInv*sum(rpar(irp_SMX:irp_SMZ)**2)) * rhoInv
    else
       eos_state % e = y(SEINT) * rhoInv
    endif

#elif defined(SDC_EVOLVE_ENTHALPY)

    if (use_tfromp) then
       ! NOT SURE IF THIS IS VALID
       eos_state % p = rpar(irp_p0)
    else
       eos_state % h = y(SENTH) * rhoInv
    endif

#endif

    ! Give the temperature an initial guess -- use the geometric mean
    ! of the minimum and maximum temperatures.

    call eos_get_small_temp(min_temp)
    call eos_get_max_temp(max_temp)

    eos_state % T = sqrt(min_temp * max_temp)

#if defined(SDC_EVOLVE_ENERGY)

    call eos(eos_input_re, eos_state)

#elif defined(SDC_EVOLVE_ENTHALPY)

    if (use_tfromp) then
       ! NOT SURE IF THIS IS VALID
       ! used to be an Abort statement
       call eos(eos_input_rp, eos_state)
    else
       call eos(eos_input_rh, eos_state)
    endif

#endif

    call eos_to_burn(eos_state, burn_state)

    burn_state % time = time

    if (rpar(irp_self_heat) > ZERO) then
       burn_state % self_heat = .true.
    else
       burn_state % self_heat = .false.       
    endif

#ifdef SDC_EVOLVE_ENTHALPY

    burn_state % p0 = rpar(irp_p0)

#endif

#ifdef NONAKA_PLOT

    burn_state % i = rpar(irp_i)
    burn_state % j = rpar(irp_j)
    burn_state % k = rpar(irp_k)

#endif

  end subroutine vode_to_burn

end module vode_type_module
