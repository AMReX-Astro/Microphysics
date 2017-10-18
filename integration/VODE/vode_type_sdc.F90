module vode_type_module

  use bl_types, only: dp_t
  use bl_constants_module

  use burn_type_module, only : burn_t, net_ienuc, eos_to_burn
  use eos_type_module, only : eos_t, eos_input_re, eos_input_rh, eos_input_rt, eos_get_small_temp, eos_get_max_temp
  use eos_module, only : eos

  use network, only : nspec, aion, aion_inv

  use rpar_indices

#if (SDC_METHOD == 1)

  use sdc_type_module, only: SEDEN, SEINT, SFS, SRHO, SMX, SMY, SMZ, SVAR, SVAR_EVOLVE, sdc_t

#elif (SDC_METHOD == 2)

  use sdc_type_module, only: SFS, SENTH, SVAR, SVAR_EVOLVE, sdc_t

#endif

  use extern_probin_module, only: renormalize_abundances

  implicit none

  private

  ! this should be larger than any reasonable temperature we will encounter
  real (kind=dp_t), parameter :: MAX_TEMP = 1.0d11

  integer, parameter :: VODE_NEQS = SVAR_EVOLVE

  public :: VODE_NEQS
  public :: clean_state, fill_unevolved_variables, &
            renormalize_species, sdc_to_vode, vode_to_sdc, &
            rhs_to_vode, jac_to_vode, vode_to_burn

contains

  subroutine clean_state(time, y, rpar)

    real(dp_t), intent(in) :: time
    real(dp_t) :: y(SVAR_EVOLVE), rpar(n_rpar_comps)

    real (kind=dp_t) :: max_e, ke

    type (eos_t) :: eos_state

    ! update rho, rho*u, etc., depending on SDC_METHOD
    call fill_unevolved_variables(time, y, rpar)

    ! Ensure that mass fractions always stay positive.
    y(SFS:SFS+nspec-1) = &
         max(min(y(SFS:SFS+nspec-1), rpar(irp_SRHO)), &
             rpar(irp_SRHO) * 1.d-200)

    ! renormalize abundances as necessary
    if (renormalize_abundances) then
       call renormalize_species(time, y, rpar)
    endif

#if (SDC_METHOD == 1)

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

    real(dp_t), intent(in) :: time
    real(dp_t) :: y(SVAR_EVOLVE), rpar(n_rpar_comps)

#if (SDC_METHOD == 1)

    ! we are always integrating from t = 0, so there is no offset
    ! time needed here.  The indexing of irp_ydot_a is based on
    ! the indices in sdc_type_module
    rpar(irp_SRHO) = rpar(irp_u_init-1+irp_SRHO) + &
         rpar(irp_ydot_a-1+SRHO) * time

    rpar(irp_SMX) = rpar(irp_u_init-1+irp_SMX) + rpar(irp_ydot_a-1+SMX) * time
    rpar(irp_SMY) = rpar(irp_u_init-1+irp_SMY) + rpar(irp_ydot_a-1+SMY) * time
    rpar(irp_SMZ) = rpar(irp_u_init-1+irp_SMZ) + rpar(irp_ydot_a-1+SMZ) * time

#elif (SDC_METHOD == 2)

    ! SDC_METHOD = 2 doesn't update unevolved variables
    ! with advective source terms, but the density is
    ! kept consistent with the partial densities.
    rpar(irp_SRHO) = sum(y(SFS:SFS - 1 + nspec))

#endif

  end subroutine fill_unevolved_variables

  subroutine renormalize_species(time, y, rpar)

    real(dp_t), intent(in) :: time
    real(dp_t) :: y(SVAR_EVOLVE), rpar(n_rpar_comps)

    real(dp_t) :: nspec_sum

    ! We only renormalize below for SDC_METHOD = 1 because
    ! in SDC_METHOD = 2, we define the density as
    ! the sum of the partial densities so nspec_sum = 1

#if (SDC_METHOD == 1)

    ! update rho, rho*u, etc., depending on SDC_METHOD
    call fill_unevolved_variables(time, y, rpar)

    nspec_sum = sum(y(SFS:SFS-1+nspec)) / rpar(irp_SRHO)

    y(SFS:SFS-1+nspec) = y(SFS:SFS-1+nspec) / nspec_sum

#endif

  end subroutine renormalize_species


  ! Given a burn state, fill the rpar and integration state data.

  subroutine sdc_to_vode(sdc, y, rpar)

    type (sdc_t) :: sdc
    real(dp_t)   :: rpar(n_rpar_comps)
    real(dp_t)   :: y(SVAR_EVOLVE)

    y(:) = sdc % y(1:SVAR_EVOLVE)

    ! advective sources
    rpar(irp_ydot_a:irp_ydot_a-1+SVAR) = sdc % ydot_a(:)

#if (SDC_METHOD == 1)

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

#elif (SDC_METHOD == 2)

    rpar(irp_p0)   = sdc % p0
    rpar(irp_SRHO) = sdc % rho

#endif

  end subroutine sdc_to_vode

  subroutine vode_to_sdc(time, y, rpar, sdc)

    real(dp_t), intent(in) :: time
    type (sdc_t) :: sdc
    real(dp_t)    :: rpar(n_rpar_comps)
    real(dp_t)    :: y(SVAR_EVOLVE)

    sdc % y(1:SVAR_EVOLVE) = y(:)

    ! unevolved state variables
    call fill_unevolved_variables(time, y, rpar)

#if (SDC_METHOD == 1)

    sdc % y(SRHO) = rpar(irp_SRHO)
    sdc % y(SMX:SMZ) = rpar(irp_SMX:irp_SMZ)

#elif (SDC_METHOD == 2)

    sdc % p0  = rpar(irp_p0)
    sdc % rho = rpar(irp_SRHO)

#endif

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

    ! add in the reacting terms -- here we convert from dY/dt to dX/dt
    ydot(SFS:SFS-1+nspec) = ydot(SFS:SFS-1+nspec) + &
         rpar(irp_SRHO) * aion(1:nspec) * burn_state % ydot(1:nspec)

#if (SDC_METHOD == 1)

    ydot(SEINT) = ydot(SEINT) + rpar(irp_SRHO) * burn_state % ydot(net_ienuc)
    ydot(SEDEN) = ydot(SEDEN) + rpar(irp_SRHO) * burn_state % ydot(net_ienuc)

#elif (SDC_METHOD == 2)

    ydot(SENTH) = ydot(SENTH) + rpar(irp_SRHO) * burn_state % ydot(net_ienuc)

#endif

  end subroutine rhs_to_vode


  subroutine jac_to_vode(time, burn_state, y, jac, rpar)

    ! this is only used with an analytic Jacobian

    real(dp_t), intent(in) :: time
    real(dp_t)    :: rpar(n_rpar_comps)
    real(dp_t)    :: y(SVAR_EVOLVE)
    type(burn_t), intent(in) :: burn_state
    real(dp_t)    :: jac(SVAR_EVOLVE,SVAR_EVOLVE)

    integer :: n

#if (SDC_METHOD == 1)

    jac(SFS:SFS+nspec-1,SFS:SFS+nspec-1) = burn_state % jac(1:nspec,1:nspec)
    jac(SFS:SFS+nspec-1,SEDEN) = burn_state % jac(1:nspec,net_ienuc)
    jac(SFS:SFS+nspec-1,SEINT) = burn_state % jac(1:nspec,net_ienuc)

    jac(SEDEN,SFS:SFS+nspec-1) = burn_state % jac(net_ienuc,1:nspec)
    jac(SEDEN,SEDEN) = burn_state % jac(net_ienuc,net_ienuc)
    jac(SEDEN,SEINT) = burn_state % jac(net_ienuc,net_ienuc)

    jac(SEINT,SFS:SFS+nspec-1) = burn_state % jac(net_ienuc,1:nspec)
    jac(SEINT,SEDEN) = burn_state % jac(net_ienuc,net_ienuc)
    jac(SEINT,SEINT) = burn_state % jac(net_ienuc,net_ienuc)

#elif (SDC_METHOD == 2)

    jac(SFS:SFS+nspec-1,SFS:SFS+nspec-1) = burn_state % jac(1:nspec,1:nspec)
    jac(SFS:SFS+nspec-1,SENTH) = burn_state % jac(1:nspec,net_ienuc)

    jac(SENTH,SFS:SFS+nspec-1) = burn_state % jac(net_ienuc,1:nspec)
    jac(SENTH,SENTH) = burn_state % jac(net_ienuc,net_ienuc)

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

#if (SDC_METHOD == 2)

    use probin_module, only: use_tfromp

#endif

    type (burn_t) :: burn_state
    real(dp_t), intent(in) :: time
    real(dp_t)    :: rpar(n_rpar_comps)
    real(dp_t)    :: y(SVAR_EVOLVE)

    type(eos_t) :: eos_state

    real(kind=dp_t) :: rhoInv, min_temp, max_temp

    ! update rho, rho*u, etc.
    call fill_unevolved_variables(time, y, rpar)

    rhoInv = ONE / rpar(irp_SRHO)

    eos_state % rho = rpar(irp_SRHO)
    eos_state % xn  = y(SFS:SFS+nspec-1) * rhoInv

#if (SDC_METHOD == 1)

    if (rpar(irp_T_from_eden) > ZERO) then
       eos_state % e = (y(SEDEN) - HALF*rhoInv*sum(rpar(irp_SMX:irp_SMZ)**2)) * rhoInv
    else
       eos_state % e = y(SEINT) * rhoInv
    endif

#elif (SDC_METHOD == 2)

    eos_state % h = y(SENTH) * rhoInv

#endif

    ! Give the temperature an initial guess -- use the geometric mean
    ! of the minimum and maximum temperatures.

    call eos_get_small_temp(min_temp)
    call eos_get_max_temp(max_temp)

    eos_state % T = sqrt(min_temp * max_temp)

#if (SDC_METHOD == 1)

    call eos(eos_input_re, eos_state)

#elif (SDC_METHOD == 2)

    if (use_tfromp) then

       call bl_error("Error: SDC_METHOD = 2 requires use_tfromp = F")

    else

       call eos(eos_input_rh, eos_state)

    endif

#endif

    call eos_to_burn(eos_state, burn_state)

    if (rpar(irp_self_heat) > ZERO) then
       burn_state % self_heat = .true.
    else
       burn_state % self_heat = .false.
    endif

#if (SDC_METHOD == 2)

    burn_state % p0 = rpar(irp_p0)

#endif

  end subroutine vode_to_burn

end module vode_type_module
