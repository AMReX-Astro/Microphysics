! Common variables and routines for burners
! that use SDC for their integration.

module actual_integrator_module

  use eos_module
  use eos_type_module
  use network
  use sdc_rpar_indices
  use burn_type_module
  use sdc_ode_module
  use sdc_type_module

  implicit none

contains

  subroutine actual_integrator_init()

    implicit none

  end subroutine actual_integrator_init


  ! Main interface

  subroutine actual_integrator(state_in, state_out, dt, time)

    !$acc routine seq

    use sdc_rpar_indices
    use amrex_fort_module, only : rt => amrex_real
    use actual_rhs_module, only : update_unevolved_species
    use temperature_integration_module, only: self_heat
    use extern_probin_module

    implicit none

    ! Input arguments

    type (burn_t), intent(in   ) :: state_in
    type (burn_t), intent(inout) :: state_out
    real(rt)     , intent(in   ) :: dt, time

    ! Local variables
    integer :: ierr

    real(rt) :: t0, t1

    type (eos_t) :: eos_state_in, eos_state_temp
    type (sdc_t) :: sdc

    real(rt) :: ener_offset
    real(rt) :: t_enuc, t_sound, limit_factor

    logical :: success

    type(sdc_t) :: sdc_state

    ! Set the tolerances.  We will be more relaxed on the temperature
    ! since it is only used in evaluating the rates.
    !
    ! **NOTE** if you reduce these tolerances, you probably will need
    ! to (a) decrease dT_crit, (b) increase the maximum number of
    ! steps allowed.

    sdc_state % atol(1:nspec_evolve) = atol_spec ! mass fractions
    sdc_state % atol(net_itemp)      = atol_temp ! temperature
    sdc_state % atol(net_ienuc)      = atol_enuc ! energy generated

    sdc_state % rtol(1:nspec_evolve) = rtol_spec ! mass fractions
    sdc_state % rtol(net_itemp)      = rtol_temp ! temperature
    sdc_state % rtol(net_ienuc)      = rtol_enuc ! energy generated


    ! Initialize the integration time.
    sdc_state % t = ZERO
    sdc_state % tmax = t0 + dt

    ! Convert our input burn state into an EOS type.

    call burn_to_eos(state_in, eos_state_in)

    ! We assume that (rho, T) coming in are valid, do an EOS call
    ! to fill the rest of the thermodynamic variables.

    call eos(eos_input_rt, eos_state_in)

    ! Convert the EOS state data into the form SDC expects.

    ! at this point, the burn_t that we care about is part of SDC
    call eos_to_sdc(eos_state_in, sdc_state % y, sdc_state % rpar)


    ! Get the internal energy e that is consistent with this T.
    ! We will start the zone with this energy and subtract it at
    ! the end. This helps a lot with convergence, rather than
    ! starting e off at zero.

    ener_offset = eos_state_in % e

    sdc_state % y(net_ienuc) = ener_offset

    ! Pass through whether we are doing self-heating.
    if (self_heat) then
       sdc_state % rpar(irp_self_heat) = ONE
    else
       sdc_state % rpar(irp_self_heat) = -ONE
    end if

    ! Copy in the zone size.

    sdc_state % rpar(irp_dx) = state_in % dx

    ! Set the sound crossing time.

    sdc_state % rpar(irp_t_sound) = state_in % dx / eos_state_in % cs

    ! set the time offset -- we integrate from 0 to dt, so this
    ! is the offset to simulation time
    
    sdc_state % rpar(irp_t0) = time


    ! If we are using the dT_crit functionality and therefore doing a
    ! linear interpolation of the specific heat in between EOS calls,
    ! do a second EOS call here to establish an initial slope.

    sdc_state % rpar(irp_Told) = eos_state_in % T

    if (dT_crit < 1.0d19) then

       call copy_eos_t(eos_state_temp, eos_state_in)
       eos_state_temp % T = eos_state_in % T * (ONE + sqrt(epsilon(ONE)))

       call eos(eos_input_rt, eos_state_temp)

       sdc_state % rpar(irp_dcvdt) = (eos_state_temp % cv - eos_state_in % cv) / &
                                     (eos_state_temp % T - eos_state_in % T)
       sdc_state % rpar(irp_dcpdt) = (eos_state_temp % cp - eos_state_in % cp) / &
                                     (eos_state_temp % T - eos_state_in % T)
    endif

    ! Save the initial state.

    sdc_state % rpar(irp_y_init:irp_y_init + neqs - 1) = sdc_state % y

    ! Call the integration routine.
    call ode(sdc, ierr)

    ! If we are using hybrid burning and the energy release was
    ! negative (or we failed), re-run this in self-heating mode.

    if ( burning_mode == 2 .and. &
         (sdc_state % y(net_ienuc) - ener_offset < ZERO .or. &
         ierr /= IERR_NONE) ) then

       sdc_state % rpar(irp_self_heat) = ONE
       sdc_state % t = ZERO
       sdc_state % tmax = dt

       call eos_to_sdc(eos_state_in, sdc_state % y, sdc_state % rpar)

       sdc_state % y(net_ienuc) = ener_offset

       ! redo the T_old, cv / cp extrapolation
       sdc_state % rpar(irp_Told) = eos_state_in % T

       if (dT_crit < 1.0d19) then
          sdc_state % rpar(irp_dcvdt) = (eos_state_temp % cv - eos_state_in % cv) / &
                                        (eos_state_temp % T - eos_state_in % T)
          sdc_state % rpar(irp_dcpdt) = (eos_state_temp % cp - eos_state_in % cp) / &
                                        (eos_state_temp % T - eos_state_in % T)
       endif

       call ode(sdc, ierr)

    endif

    ! If we still failed, print out the current state of the integration.

    if (ierr /= IERR_NONE) then

#ifndef ACC
       print *, 'ERROR: integration failed in net'
       print *, 'ierr = ', ierr
       print *, 'dt = ', dt
       print *, 'time start = ', time
       print *, 'time current = ', sdc % t
       print *, 'dens = ', state_in % rho
       print *, 'temp start = ', state_in % T
       print *, 'xn start = ', state_in % xn
       print *, 'temp current = ', sdc % y(net_itemp)
       print *, 'xn current = ', sdc % y(1:nspec_evolve), &
            sdc_state % rpar(irp_nspec:irp_nspec+n_not_evolved-1)
       print *, 'energy generated = ', sdc % y(net_ienuc) - ener_offset
#endif

       state_out % success = .false.
       return

    endif

    ! Subtract the energy offset.

    sdc_state % y(net_ienuc) = sdc_state % y(net_ienuc) - ener_offset

    ! Store the final data, and then normalize abundances.
    call sdc_to_burn(sdc_state % y, sdc_state % rpar, state_out)

    ! cache the success
    success = state_out % success

    state_out % success = success

    if (nspec_evolve < nspec) then
       call update_unevolved_species(state_out)
    endif

    ! For burning_mode == 3, limit the burning.

    if (burning_mode == 3) then

       t_enuc = eos_state_in % e / max(abs(state_out % e - state_in % e) / max(dt, 1.d-50), 1.d-50)
       t_sound = state_in % dx / eos_state_in % cs

       limit_factor = min(1.0d0, burning_mode_factor * t_enuc / t_sound)

       state_out % e = state_in % e + limit_factor * (state_out % e - state_in % e)
       state_out % xn(:) = state_in % xn(:) + limit_factor * (state_out % xn(:) - state_in % xn(:))

    endif

    call normalize_abundances_burn(state_out)

#ifndef ACC
    if (burner_verbose) then

       ! Print out some integration statistics, if desired.

       print *, 'integration summary: '
       print *, 'dt = ', dt
       print *, 'time start = ', time
       print *, 'time final = ', sdc_state % t
       print *, 'dens = ', state_out % rho
       print *, 'temp start = ', state_in % T
       print *, 'xn start = ', state_in % xn
       print *, 'temp final = ', state_out % T
       print *, 'xn final = ', state_out % xn
       print *, 'energy generated = ', state_out % e - state_in % e
       print *, 'number of steps taken: ', sdc_state % n
       print *, 'number of RHS evaluations: ', state_out % n_rhs
       print *, 'number of Jacobian evaluations: ', state_out % n_jac

    endif
#endif

  end subroutine actual_integrator

end module actual_integrator_module
