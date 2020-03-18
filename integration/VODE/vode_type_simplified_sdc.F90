module vode_type_module

  use amrex_fort_module, only: rt => amrex_real
  use amrex_constants_module
  use cuvode_parameters_module, only : VODE_NEQS
  use cuvode_types_module, only : dvode_t
  use network, only : nspec, aion, aion_inv

  use vode_rpar_indices
  use sdc_type_module

  use extern_probin_module, only: renormalize_abundances

  implicit none

  ! this should be larger than any reasonable temperature we will encounter   
  real (kind=rt), parameter :: MAX_TEMP = 1.0e11_rt          

  public

contains

  subroutine clean_state(time, vode_state)

    use eos_type_module, only : eos_t, eos_input_re, eos_input_rt
    use eos_module, only : eos

    real(rt), intent(in) :: time
    type(dvode_t), intent(inout) :: vode_state

    real(rt) :: max_e, ke

    type (eos_t) :: eos_state

    !$gpu

    ! update rho, rho*u, etc.
    call fill_unevolved_variables(time, vode_state)

    ! Ensure that mass fractions always stay positive.
    vode_state % y(SFS:SFS+nspec-1) = &
         max(min(vode_state % y(SFS:SFS+nspec-1), vode_State % rpar(irp_SRHO)), &
             vode_state % rpar(irp_SRHO) * 1.e-200_rt)

    ! renormalize abundances as necessary
    if (renormalize_abundances) then
       call renormalize_species(time, vode_state)
    endif

#ifdef SDC_EVOLVE_ENERGY

    ! Ensure that internal energy never goes above the maximum limit
    ! provided by the EOS. Same for the internal energy implied by the
    ! total energy (which we get by subtracting kinetic energy).
    eos_state % rho = vode_state % rpar(irp_SRHO)
    eos_state % T = MAX_TEMP
    eos_state % xn = vode_state % y(SFS:SFS+nspec-1) / vode_state % rpar(irp_SRHO)

    call eos(eos_input_rt, eos_state)

    max_e = eos_state % e

    vode_state % y(SEINT) = min(vode_state % rpar(irp_SRHO) * max_e, vode_state % y(SEINT))

    ke = vode_state % y(SEDEN) - &
         HALF * sum(vode_state % rpar(irp_SMX:irp_SMZ)**2) / vode_state % rpar(irp_SRHO)

    vode_state % y(SEDEN) = min(vode_state % rpar(irp_SRHO) * max_e + ke, vode_state % y(SEDEN))

#endif

  end subroutine clean_state


  subroutine fill_unevolved_variables(time, vode_state)

    real(rt), intent(in) :: time
    type(dvode_t), intent(inout) :: vode_state

    !$gpu

#if defined(SDC_EVOLVE_ENERGY)

    ! we are always integrating from t = 0, so there is no offset
    ! time needed here.  The indexing of irp_ydot_a is based on
    ! the indices in sdc_type_module
    vode_state % rpar(irp_SRHO) = vode_state % rpar(irp_u_init-1+irp_SRHO) + &
         vode_state % rpar(irp_ydot_a-1+SRHO) * time

    vode_state % rpar(irp_SMX) = vode_state % rpar(irp_u_init-1+irp_SMX) + &
         vode_state % rpar(irp_ydot_a-1+SMX) * time
    vode_state % rpar(irp_SMY) = vode_state % rpar(irp_u_init-1+irp_SMY) + &
         vode_state % rpar(irp_ydot_a-1+SMY) * time
    vode_state % rpar(irp_SMZ) = vode_state % rpar(irp_u_init-1+irp_SMZ) + &
         vode_state % rpar(irp_ydot_a-1+SMZ) * time

#elif defined(SDC_EVOLVE_ENTHALPY)

    ! Keep density consistent with the partial densities.
    vode_state % rpar(irp_SRHO) = sum(vode_state % y(SFS:SFS - 1 + nspec))

#endif

  end subroutine fill_unevolved_variables

  subroutine renormalize_species(time, vode_state)

    real(rt), intent(in) :: time
    type(dvode_t), intent(inout) :: vode_state

    real(rt) :: nspec_sum

    !$gpu

    ! We only renormalize species when evolving energy because
    ! when we evolve enthalpy, we define the density as
    ! the sum of the partial densities rho*X for each species.
#ifdef SDC_EVOLVE_ENERGY

    ! update rho, rho*u, etc.
    call fill_unevolved_variables(time, vode_state)

    nspec_sum = sum(vode_state % y(SFS:SFS-1+nspec)) / vode_state % rpar(irp_SRHO)

    vode_state % y(SFS:SFS-1+nspec) = vode_state % y(SFS:SFS-1+nspec) / nspec_sum

#endif

  end subroutine renormalize_species


  ! Given a burn state, fill the rpar and integration state data.

  subroutine sdc_to_vode(sdc, vode_state)

    type (sdc_t) :: sdc
    type(dvode_t), intent(inout) :: vode_state

    !$gpu

    vode_state % y(:) = sdc % y(1:SVAR_EVOLVE)

    ! advective sources
    vode_state % rpar(irp_ydot_a:irp_ydot_a-1+SVAR) = sdc % ydot_a(:)

#if defined(SDC_EVOLVE_ENERGY)

    ! unevolved state variables
    vode_state % rpar(irp_SRHO) = sdc % y(SRHO)
    vode_state % rpar(irp_SMX:irp_SMZ) = sdc % y(SMX:SMZ)

    ! initial state for unevolved variables
    vode_state % rpar(irp_u_init-1+irp_SRHO) = sdc % y(SRHO)
    vode_state % rpar(irp_u_init-1+irp_SMX:irp_u_init-1+irp_SMZ) = sdc % y(SMX:SMZ)

    ! other parameters
    if (sdc % T_from_eden) then
       vode_state % rpar(irp_T_from_eden) = ONE
    else
       vode_state % rpar(irp_T_from_eden) = -ONE
    endif

#elif defined(SDC_EVOLVE_ENTHALPY)

    vode_state % rpar(irp_p0)   = sdc % p0
    vode_state % rpar(irp_SRHO) = sdc % rho

#endif

#ifdef NONAKA_PLOT

    ! bookkeeping information
    vode_state % rpar(irp_i) = sdc % i
    vode_state % rpar(irp_j) = sdc % j
    vode_state % rpar(irp_k) = sdc % k
    vode_state % rpar(irp_iter) = sdc % sdc_iter

#endif

  end subroutine sdc_to_vode

  subroutine vode_to_sdc(time, vode_state, sdc)

    real(rt), intent(in) :: time
    type (sdc_t) :: sdc
    type(dvode_t), intent(inout) :: vode_state

    !$gpu

    sdc % y(1:SVAR_EVOLVE) = vode_state % y(:)

    ! unevolved state variables
    call fill_unevolved_variables(time, vode_state)

#if defined(SDC_EVOLVE_ENERGY)

    sdc % y(SRHO) = vode_state % rpar(irp_SRHO)
    sdc % y(SMX:SMZ) = vode_state % rpar(irp_SMX:irp_SMZ)

#elif defined(SDC_EVOLVE_ENTHALPY)

    sdc % p0  = vode_state % rpar(irp_p0)
    sdc % rho = vode_state % rpar(irp_SRHO)

#endif

#ifdef NONAKA_PLOT

    sdc % i = vode_state % rpar(irp_i)
    sdc % j = vode_state % rpar(irp_j)
    sdc % k = vode_state % rpar(irp_k)
    sdc % sdc_iter = vode_state % rpar(irp_iter)

#endif

  end subroutine vode_to_sdc


  subroutine rhs_to_vode(time, burn_state, ydot_react, vode_state, ydot)

    use burn_type_module, only : burn_t, net_ienuc, neqs

    real(rt), intent(in) :: time
    type(dvode_t), intent(inout) :: vode_state
    type(burn_t), intent(in) :: burn_state
    real(rt), intent(in) :: ydot_react(neqs)
    real(rt), intent(out) :: ydot(SVAR_EVOLVE)

    !$gpu

    call fill_unevolved_variables(time, vode_state)

    ! ydot_react has just the contribution to the RHS from the
    ! reaction network.  Note that these are in terms of dY/dt

    ! start with the contribution from the non-reacting sources
    ydot(:) = vode_state % rpar(irp_ydot_a:irp_ydot_a-1+SVAR_EVOLVE)

    ! add in the reacting terms -- here we convert from dY/dt to dX/dt
    ydot(SFS:SFS-1+nspec) = ydot(SFS:SFS-1+nspec) + &
         vode_state % rpar(irp_SRHO) * aion(1:nspec) * ydot_react(1:nspec)

#if defined(SDC_EVOLVE_ENERGY)

    ydot(SEINT) = ydot(SEINT) + vode_state % rpar(irp_SRHO) * ydot_react(net_ienuc)
    ydot(SEDEN) = ydot(SEDEN) + vode_state % rpar(irp_SRHO) * ydot_react(net_ienuc)

#elif defined(SDC_EVOLVE_ENTHALPY)

    ydot(SENTH) = ydot(SENTH) + vode_state % rpar(irp_SRHO) * ydot_react(net_ienuc)

#endif

  end subroutine rhs_to_vode


  subroutine jac_to_vode(time, jac_react, vode_state, jac)

    ! this is only used with an analytic Jacobian

    use burn_type_module, only : net_ienuc, neqs

    real(rt), intent(in) :: time
    type(dvode_t), intent(inout) :: vode_state
    real(rt), intent(in) :: jac_react(neqs, neqs)
    real(rt)    :: jac(SVAR_EVOLVE,SVAR_EVOLVE)

    integer :: n

    !$gpu

#if defined(SDC_EVOLVE_ENERGY)

    jac(SFS:SFS+nspec-1,SFS:SFS+nspec-1) = jac_react(1:nspec,1:nspec)
    jac(SFS:SFS+nspec-1,SEDEN) = jac_react(1:nspec,net_ienuc)
    jac(SFS:SFS+nspec-1,SEINT) = jac_react(1:nspec,net_ienuc)

    jac(SEDEN,SFS:SFS+nspec-1) = jac_react(net_ienuc,1:nspec)
    jac(SEDEN,SEDEN) = jac_react(net_ienuc,net_ienuc)
    jac(SEDEN,SEINT) = jac_react(net_ienuc,net_ienuc)

    jac(SEINT,SFS:SFS+nspec-1) = jac_react(net_ienuc,1:nspec)
    jac(SEINT,SEDEN) = jac_react(net_ienuc,net_ienuc)
    jac(SEINT,SEINT) = jac_react(net_ienuc,net_ienuc)

#elif defined(SDC_EVOLVE_ENTHALPY)

    jac(SFS:SFS+nspec-1,SFS:SFS+nspec-1) = jac_react(1:nspec,1:nspec)
    jac(SFS:SFS+nspec-1,SENTH) = jac_react(1:nspec,net_ienuc)

    jac(SENTH,SFS:SFS+nspec-1) = jac_react(net_ienuc,1:nspec)
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


  subroutine vode_to_burn(time, vode_state, burn_state)

    use eos_type_module, only : eos_t, eos_input_re, eos_input_rt, eos_input_rp, eos_input_rh
    use eos_type_module, only : eos_get_small_temp, eos_get_max_temp
    use eos_module, only : eos
    use burn_type_module, only : eos_to_burn, burn_t

#ifdef SDC_EVOLVE_ENTHALPY
    use meth_params_module, only: use_tfromp
#endif

    type (burn_t) :: burn_state
    real(rt), intent(in) :: time
    type(dvode_t), intent(inout) :: vode_state

    type(eos_t) :: eos_state

    real(rt) :: rhoInv, min_temp, max_temp

    !$gpu

    ! update rho, rho*u, etc.
    call fill_unevolved_variables(time, vode_state)

    rhoInv = ONE / vode_state % rpar(irp_SRHO)

    eos_state % rho = vode_state % rpar(irp_SRHO)
    eos_state % xn  = vode_state % y(SFS:SFS+nspec-1) * rhoInv

#if defined(SDC_EVOLVE_ENERGY)

    if (vode_state % rpar(irp_T_from_eden) > ZERO) then
       eos_state % e = (vode_state % y(SEDEN) - &
            HALF*rhoInv*sum(vode_state % rpar(irp_SMX:irp_SMZ)**2)) * rhoInv
    else
       eos_state % e = vode_state % y(SEINT) * rhoInv
    endif

#elif defined(SDC_EVOLVE_ENTHALPY)

    if (use_tfromp) then
       ! NOT SURE IF THIS IS VALID
       eos_state % p = vode_state % rpar(irp_p0)
    else
       eos_state % h = vode_state % y(SENTH) * rhoInv
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

    if (vode_state % rpar(irp_self_heat) > ZERO) then
       burn_state % self_heat = .true.
    else
       burn_state % self_heat = .false.
    endif

#ifdef SDC_EVOLVE_ENTHALPY

    burn_state % p0 = vode_state % rpar(irp_p0)

#endif

#ifdef NONAKA_PLOT

    burn_state % i = vode_state % rpar(irp_i)
    burn_state % j = vode_state % rpar(irp_j)
    burn_state % k = vode_state % rpar(irp_k)

#endif

  end subroutine vode_to_burn

end module vode_type_module
