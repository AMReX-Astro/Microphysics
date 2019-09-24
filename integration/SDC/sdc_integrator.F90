! Common variables and routines for burners
! that use SDC for their integration.

module sdc_integrator_module

  use eos_module
  use eos_type_module
  use network
  use sdc_rpar_indices
  use burn_type_module
  use stiff_ode
  use sdc_type_module

  implicit none

contains

  subroutine sdc_integrator_init()

    implicit none

  end subroutine sdc_integrator_init


  ! Main interface

  subroutine sdc_integrator(state_in, state_out, dt, time, status)

    !$acc routine seq

    use sdc_rpar_indices
    use amrex_fort_module, only : rt => amrex_real
    use extern_probin_module, only: burner_verbose, burning_mode, burning_mode_factor, dT_crit
    use actual_rhs_module, only : update_unevolved_species
    use integration_data, only: integration_status_t
    use temperature_integration_module, only: self_heat

    implicit none

    ! Input arguments

    type (burn_t), intent(in   ) :: state_in
    type (burn_t), intent(inout) :: state_out
    real(rt)     , intent(in   ) :: dt, time
    type (integration_status_t), intent(inout) :: status

    ! Local variables
    integer :: ierr

    real(rt) :: atol(neqs), rtol(neqs)   ! input state, abs and rel tolerances
    real(rt) :: t0, t1

    type (eos_t) :: eos_state_in, eos_state_temp
    type (sdc_t) :: sdc

    real(rt) :: ener_offset
    real(rt) :: t_enuc, t_sound, limit_factor

    logical :: success

    ! Set the tolerances.  We will be more relaxed on the temperature
    ! since it is only used in evaluating the rates.
    !
    ! **NOTE** if you reduce these tolerances, you probably will need
    ! to (a) decrease dT_crit, (b) increase the maximum number of
    ! steps allowed.

    atol(1:nspec_evolve) = status % atol_spec ! mass fractions
    atol(net_itemp)      = status % atol_temp ! temperature
    atol(net_ienuc)      = status % atol_enuc ! energy generated

    rtol(1:nspec_evolve) = status % rtol_spec ! mass fractions
    rtol(net_itemp)      = status % rtol_temp ! temperature
    rtol(net_ienuc)      = status % rtol_enuc ! energy generated

    ! Note that at present, we use a uniform error tolerance chosen
    ! to be the largest of the relative error tolerances for any
    ! equation. We may expand this capability in the future.

    sdc % atol = atol
    sdc % rtol = rtol

    ! Start out by assuming a successful burn.

    state_out % success = .true.

    ! Initialize the integration time.
    t0 = ZERO
    t1 = t0 + dt

    ! Convert our input burn state into an EOS type.

    call burn_to_eos(state_in, eos_state_in)

    ! We assume that (rho, T) coming in are valid, do an EOS call
    ! to fill the rest of the thermodynamic variables.

    call eos(eos_input_rt, eos_state_in)

    ! Convert the EOS state data into the form SDC expects.

    ! at this point, the burn_t that we care about is part of SDC
    call eos_to_sdc(eos_state_in, sdc)

    sdc % burn_s % i = state_in % i
    sdc % burn_s % j = state_in % j
    sdc % burn_s % k = state_in % k

    sdc % burn_s % n_rhs = 0
    sdc % burn_s % n_jac = 0

    ! Get the internal energy e that is consistent with this T.
    ! We will start the zone with this energy and subtract it at
    ! the end. This helps a lot with convergence, rather than
    ! starting e off at zero.

    ener_offset = eos_state_in % e

    sdc % y(net_ienuc) = ener_offset

    ! Pass through whether we are doing self-heating.

    sdc % burn_s % self_heat = self_heat

    ! Copy in the zone size.

    sdc % burn_s % dx = state_in % dx

    ! Set the sound crossing time.

    sdc % upar(irp_t_sound) = state_in % dx / eos_state_in % cs

    ! set the time offset -- we integrate from 0 to dt, so this
    ! is the offset to simulation time
    
    sdc % upar(irp_t0) = time


    ! If we are using the dT_crit functionality and therefore doing a
    ! linear interpolation of the specific heat in between EOS calls,
    ! do a second EOS call here to establish an initial slope.

    sdc % burn_s % T_old = eos_state_in % T

    if (dT_crit < 1.0d19) then

       eos_state_temp = eos_state_in
       eos_state_temp % T = eos_state_in % T * (ONE + sqrt(epsilon(ONE)))

       call eos(eos_input_rt, eos_state_temp)

       sdc % burn_s % dcvdt = (eos_state_temp % cv - eos_state_in % cv) / &
            (eos_state_temp % T - eos_state_in % T)

       sdc % burn_s % dcpdt = (eos_state_temp % cp - eos_state_in % cp) / &
            (eos_state_temp % T - eos_state_in % T)
    endif

    ! Save the initial state.

    sdc % upar(irp_y_init:irp_y_init + neqs - 1) = sdc % y

    ! Call the integration routine.
    call ode(sdc, t0, t1, maxval(rtol), ierr)

    ! If we are using hybrid burning and the energy release was
    ! negative (or we failed), re-run this in self-heating mode.

    if ( burning_mode == 2 .and. &
         (sdc % y(net_ienuc) - ener_offset < ZERO .or. &
         ierr /= IERR_NONE) ) then

       sdc % burn_s % self_heat = .true.

       call eos_to_sdc(eos_state_in, sdc)

       sdc % y(net_ienuc) = ener_offset

       ! redo the T_old, cv / cp extrapolation
       sdc % burn_s % T_old = eos_state_in % T

       if (dT_crit < 1.0d19) then
          sdc % burn_s % dcvdt = (eos_state_temp % cv - eos_state_in % cv) / &
               (eos_state_temp % T - eos_state_in % T)

          sdc % burn_s % dcpdt = (eos_state_temp % cp - eos_state_in % cp) / &
               (eos_state_temp % T - eos_state_in % T)
       endif

       call ode(sdc, t0, t1, maxval(rtol), ierr)

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
            sdc % upar(irp_nspec:irp_nspec+n_not_evolved-1)
       print *, 'energy generated = ', sdc % y(net_ienuc) - ener_offset
#endif

       state_out % success = .false.
       return

    endif

    ! Subtract the energy offset.

    sdc % y(net_ienuc) = sdc % y(net_ienuc) - ener_offset

    ! Store the final data, and then normalize abundances.
    call sdc_to_burn(sdc)

    ! cache the success
    success = state_out % success

    state_out = sdc % burn_s

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
       print *, 'time final = ', sdc % t
       print *, 'dens = ', state_out % rho
       print *, 'temp start = ', state_in % T
       print *, 'xn start = ', state_in % xn
       print *, 'temp final = ', state_out % T
       print *, 'xn final = ', state_out % xn
       print *, 'energy generated = ', state_out % e - state_in % e
       print *, 'number of steps taken: ', sdc % n
       print *, 'number of RHS evaluations: ', state_out % n_rhs
       print *, 'number of Jacobian evaluations: ', state_out % n_jac

    endif
#endif

  end subroutine sdc_integrator

end module sdc_integrator_module
