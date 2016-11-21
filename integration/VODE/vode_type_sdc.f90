module vode_type_module

  use bl_types, only: dp_t
  use rpar_indices, only: n_rpar_comps
  use sdc_type_module, only: SVAR, SVAR_EVOLVE

  implicit none

  private

  integer, parameter :: VODE_NEQS = SVAR_EVOLVE

  public :: VODE_NEQS

contains

  subroutine clean_state(time, y, rpar)

    real(dp_t) :: y(SVAR_EVOLVE), rpar(n_rpar_comps)

    real (kind=dp_t) :: max_e, ke

    type (eos_t) :: eos_state

    ! update rho, rho*u, etc.
    call fill_unevolved_variables(time, y, rpar)

    ! Ensure that mass fractions always stay positive.
    y(SFS:SFS+nspec-1) = &
         max(min(y(SFS:SFS+nspec-1), rpar(irp_SRHO)), &
             rpar(irp_SRHO) * 1.d-200)

    ! renormalize abundances as necessary
    if (renormalize_abundances) then
       call renormalize_species(y, rpar)
    endif

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

  end subroutine clean_state


  subroutine fill_unevolved_variables(time, y, rpar)

    real(dp_t), intent(in) :: time
    real(dp_t) :: y(SVAR_EVOLVE), rpar(n_rpar_comps)

    ! we are always integrating from t = 0, so there is no offset
    ! time needed here
    rpar(irp_SRHO) = rpar(irp_u_init-1+irp_SRHO) + &
         rpar(irp_ydot_a-1+irp_SRHO) * time

    do m = irp_SMX, irp_SMZ
       rpar(m) = rpar(irp_u_init-1+m) + rpar(irp_ydot_a-1+m) * time
    enddo

  end subroutine fill_unevolved_variables

  subroutine renormalize_species(time, y, rpar)

    real(dp_t), intent(in) :: time
    real(dp_t) :: y(SVAR_EVOLVE), rpar(n_rpar_comps)

    real(dp_t) :: nspec_sum

    ! update rho, rho*u, etc.
    call fill_unevolved_variables(time, y, rpar)

    nspec_sum = sum(y(SFS:SFS-1+nspec)) / rpar(irp_SRHO)

    y(SFS:SFS-1+nspec) = y(SFS:SFS-1+nspec) / nspec_sum

  end subroutine renormalize_species


  ! Given a burn state, fill the rpar and integration state data.

  subroutine sdc_to_vode(sdc, y, rpar)

    type (sdc_t) :: sdc
    real(dp_t    :: rpar(n_rpar_comps)
    real(dp_t    :: y(SVAR_EVOLVE)

    y(:) = sdc % y(1:SVAR_EVOLVE)

    ! unevolved state variables
    rpar(irp_dens) = sdc % y(SRHO)
    rpar(irp_SMX:irp_SMZ) = sdc % y(SMX:SMZ)

    ! advective sources
    rpar(irp_ydot_a:irp_ydot_a-1+SVAR) = sdc % ydot_a(:)

    ! initial state for unevolved variables
    rpar(irp_u_init-1+irp_SRHO) = sdc % y(SRHO)
    rpar(irp_u_init-1+irp_SMX:irp_u_init-1+irp_SMZ) = sdc % y(SMX:SMZ)

    ! other parameters
    if (sdc % T_from_eden) then
       rpar(irp_T_from_eden) = ONE
    else
       rpar(irp_T_from_eden) = -ONE
    endif

  end subroutine sdc_to_vode

  subroutine vode_to_sdc(sdc, y, rpar)

    type (sdc_t) :: sdc
    real(dp_t)    :: rpar(n_rpar_comps)
    real(dp_t)    :: y(SVAR_EVOLVE)

    sdc % y(1:SVAR_EVOLVE) = y(:)

    ! unevolved state variables
    call fill_unevolved_variables(time, y, rpar)

    sdc % y(SRHO) = rpar(irp_dens)
    sdc % y(SMX:SMZ) = rpar(irp_SMX:irp_SMZ)

  end subroutine vode_to_sdc


  subroutine rhs_to_vode(time, burn_state, y, ydot, rpar)

    real(dp_t), intent(in) :: time
    real(dp_t)    :: rpar(n_rpar_comps)
    real(dp_t)    :: y(SVAR_EVOLVE), ydot(SVAR_EVOLVE)
    type(burn_t), intent(in) :: burn_state

    call fill_unevolved_variables(time, y, rpar)

    ! burn_t % ydot has just the contribution to the RHS from the
    ! reaction network.  Note that these are in terms of dY/dt

    ! start with the contribution from the non-reacting sources
    ydot(:) = rpar(irp_ydot_a:irp_ydot_a-1+SVAR_EVOLVE)

    ! add in the reacting terms
    ydot(SFS:SFS-1+nspec) = ydot(SFS:SFS-1+nspec) + &
         rpar(irp_SRHO) * aion(1:nspec) * burn_state % ydot(1:nspec)

    ydot(SEINT) = ydot(SEINT) + rpar(irp_SRHO) * burn_state % ydot(net_ienuc)
    ydot(SEDEN) = ydot(SEDEN) + rpar(irp_SRHO) * burn_state % ydot(net_ienuc)

  end subroutine rhs_to_vode


  subroutine jac_to_vode(time, burn_state, y, jac, rpar)

    ! this is only used with an analytic Jacobian

    real(dp_t), intent(in) :: time
    real(dp_t)    :: rpar(n_rpar_comps)
    real(dp_t)    :: y(SVAR_EVOLVE)
    type(burn_t), intent(in) :: burn_state
    real(dp_t)    :: jac(SVAR_EVOLVE,SVAR_EVOLVE)

    jac(SFS:SFS+nspec-1,SFS:SFS+nspec-1) = burn_state % jac(1:nspec,1:nspec)
    jac(SFS:SFS+nspec-1,SEDEN) = burn_state % jac(1:nspec,net_ienuc)
    jac(SFS:SFS+nspec-1,SEINT) = burn_state % jac(1:nspec,net_ienuc)

    jac(SEDEN,SFS:SFS+nspec-1) = burn_state % jac(net_ienuc,1:nspec)
    jac(SEDEN,SEDEN) = burn_state % jac(net_ienuc,net_ienuc)
    jac(SEDEN,SEINT) = burn_state % jac(net_ienuc,net_ienuc)

    jac(SEINT,SFS:SFS+nspec-1) = burn_state % jac(net_ienuc,1:nspec)
    jac(SEINT,SEDEN) = burn_state % jac(net_ienuc,net_ienuc)
    jac(SEINT,SEINT) = burn_state % jac(net_ienuc,net_ienuc)

    ! Scale it to match our variables. We don't need to worry about
    ! the rho dependence, since every one of the SDC variables is
    ! linear in rho, so we just need to focus on the Y --> X
    ! conversion.
    do n = 1, nspec
       jac(SFS+n-1,:) = jac(SFS+n-1,:) * aion(n)
       jac(:,SFS+n-1) = jac(:,SFS+n-1) / aion(n)
    enddo

  end subroutine jac_to_vode


  subroutine vode_to_burn(time, y, rpar, state)

    type (burn_t) :: state
    real(dp_t), intent(in) :: time
    real(dp_t)    :: rpar(n_rpar_comps)
    real(dp_t)    :: y(SVAR_EVOLVE)

    integer :: n

    ! update rho, rho*u, etc.
    call fill_unevolved_variables(time, y, rpar)

    rhoInv = ONE / rpar(irp_SRHO)

    eos_state % rho = rpar(irp_SRHO)
    eos_state % xn  = y(SFS:SFS+nspec-1) * rhoInv

    if (rpar(irp_T_from_eden > ZERO) then
       eos_state % e = (y(SEDEN) - HALF*rhoInv*sum(rpar(irp_SMX:irp_SMZ)**2)) * rhoInv
    else
       eos_state % e = y(SEINT) * rhoInv
    endif

    ! Give the temperature an initial guess -- use the geometric mean
    ! of the minimum and maximum temperatures.

    call eos_get_small_temp(min_temp)
    call eos_get_max_temp(max_temp)

    eos_state % T = sqrt(min_temp * max_temp)

    call eos(eos_input_re, eos_state)
    call eos_to_burn(eos_state, burn)

    burn % time = time

    if (rpar(irp_self_heat) > ZERO) then
       burn % self_heat = .true.
    else
       burn % self_heat = .false.       
    endif

  end subroutine vode_to_burn

end module vode_type_module
