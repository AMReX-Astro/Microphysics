! Common variables and routines for burners
! that use BS for their integration.

module bs_integrator_module

  use eos_module
  use eos_type_module
  use network
  use bs_rpar_indices
  use burn_type_module
  use stiff_ode
  use bs_type_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

contains

  subroutine bs_integrator_init()

    implicit none

    nseq = [2, 6, 10, 14, 22, 34, 50, 70]

    !$acc update device(nseq)

  end subroutine bs_integrator_init



  ! Main interface

  subroutine bs_integrator(state_in, state_out, dt, time, status)

    !$acc routine seq

    use bs_rpar_indices
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
    type (bs_t) :: bs

    real(rt) :: ener_offset
    real(rt) :: t_enuc, t_sound, limit_factor

    logical :: success

    ! Set the tolerances.  We will be more relaxed on the temperature
    ! since it is only used in evaluating the rates.
    !
    ! **NOTE** if you reduce these tolerances, you probably will need
    ! to (a) decrease dT_crit, (b) increase the maximum number of
    ! steps allowed.

    atol(1:nspec)   = status % atol_spec ! mass fractions
    atol(net_itemp) = status % atol_temp ! temperature
    atol(net_ienuc) = status % atol_enuc ! energy generated

    rtol(1:nspec)   = status % rtol_spec ! mass fractions
    rtol(net_itemp) = status % rtol_temp ! temperature
    rtol(net_ienuc) = status % rtol_enuc ! energy generated

    ! Note that at present, we use a uniform error tolerance chosen
    ! to be the largest of the relative error tolerances for any
    ! equation. We may expand this capability in the future.

    bs % atol = atol
    bs % rtol = rtol

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

    ! Convert the EOS state data into the form BS expects.

    ! at this point, the burn_t that we care about is part of bs
    call eos_to_bs(eos_state_in, bs)

    bs % burn_s % i = state_in % i
    bs % burn_s % j = state_in % j
    bs % burn_s % k = state_in % k

    bs % burn_s % n_rhs = 0
    bs % burn_s % n_jac = 0

    ! Get the internal energy e that is consistent with this T.
    ! We will start the zone with this energy and subtract it at
    ! the end. This helps a lot with convergence, rather than
    ! starting e off at zero.

    ener_offset = eos_state_in % e

    bs % y(net_ienuc) = ener_offset

    ! Pass through whether we are doing self-heating.

    bs % burn_s % self_heat = self_heat

    ! Copy in the zone size.

    bs % burn_s % dx = state_in % dx

    ! Set the sound crossing time.

    bs % upar(irp_t_sound) = state_in % dx / eos_state_in % cs

    ! set the time offset -- we integrate from 0 to dt, so this
    ! is the offset to simulation time
    
    bs % upar(irp_t0) = time


    ! If we are using the dT_crit functionality and therefore doing a
    ! linear interpolation of the specific heat in between EOS calls,
    ! do a second EOS call here to establish an initial slope.

    bs % burn_s % T_old = eos_state_in % T

    if (dT_crit < 1.0e19_rt) then

       eos_state_temp = eos_state_in
       eos_state_temp % T = eos_state_in % T * (ONE + sqrt(epsilon(ONE)))

       call eos(eos_input_rt, eos_state_temp)

       bs % burn_s % dcvdt = (eos_state_temp % cv - eos_state_in % cv) / &
            (eos_state_temp % T - eos_state_in % T)

       bs % burn_s % dcpdt = (eos_state_temp % cp - eos_state_in % cp) / &
            (eos_state_temp % T - eos_state_in % T)
    endif

    ! Call the integration routine.
    call ode(bs, t0, t1, maxval(rtol), ierr)

    ! If we are using hybrid burning and the energy release was
    ! negative (or we failed), re-run this in self-heating mode.

    if ( burning_mode == 2 .and. &
         (bs % y(net_ienuc) - ener_offset < ZERO .or. &
         ierr /= IERR_NONE) ) then

       bs % burn_s % self_heat = .true.

       call eos_to_bs(eos_state_in, bs)

       bs % y(net_ienuc) = ener_offset

       ! redo the T_old, cv / cp extrapolation
       bs % burn_s % T_old = eos_state_in % T

       if (dT_crit < 1.0e19_rt) then
          bs % burn_s % dcvdt = (eos_state_temp % cv - eos_state_in % cv) / &
               (eos_state_temp % T - eos_state_in % T)

          bs % burn_s % dcpdt = (eos_state_temp % cp - eos_state_in % cp) / &
               (eos_state_temp % T - eos_state_in % T)
       endif

       call ode(bs, t0, t1, maxval(rtol), ierr)

    endif

    ! If we still failed, print out the current state of the integration.

    if (ierr /= IERR_NONE) then

#ifndef CUDA
       print *, 'ERROR: integration failed in net'
       print *, 'ierr = ', ierr
       print *, 'dt = ', dt
       print *, 'time start = ', time
       print *, 'time current = ', bs % t
       print *, 'dens = ', state_in % rho
       print *, 'temp start = ', state_in % T
       print *, 'xn start = ', state_in % xn
       print *, 'temp current = ', bs % y(net_itemp)
       print *, 'xn current = ', bs % y(1:nspec)
       print *, 'energy generated = ', bs % y(net_ienuc) - ener_offset
#endif

       state_out % success = .false.
       return

    endif

    ! Subtract the energy offset.

    bs % y(net_ienuc) = bs % y(net_ienuc) - ener_offset

    ! Store the final data, and then normalize abundances.
    call bs_to_burn(bs)

    ! cache the success
    success = state_out % success

    state_out = bs % burn_s

    state_out % success = success

    ! For burning_mode == 3, limit the burning.

    if (burning_mode == 3) then

       t_enuc = eos_state_in % e / max(abs(state_out % e - state_in % e) / max(dt, 1.e-50_rt), 1.e-50_rt)
       t_sound = state_in % dx / eos_state_in % cs

       limit_factor = min(1.0e0_rt, burning_mode_factor * t_enuc / t_sound)

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
       print *, 'time final = ', bs % t
       print *, 'dens = ', state_out % rho
       print *, 'temp start = ', state_in % T
       print *, 'xn start = ', state_in % xn
       print *, 'temp final = ', state_out % T
       print *, 'xn final = ', state_out % xn
       print *, 'energy generated = ', state_out % e - state_in % e
       print *, 'number of steps taken: ', bs % n
       print *, 'number of RHS evaluations: ', state_out % n_rhs
       print *, 'number of Jacobian evaluations: ', state_out % n_jac

    endif
#endif

  end subroutine bs_integrator

end module bs_integrator_module
