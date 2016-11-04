! Common variables and routines for burners
! that use BS for their integration.

module actual_integrator_module

  use eos_module
  use network
  use rpar_indices
  use burn_type_module
  use bl_types
  use stiff_ode
  use bs_type_module

  implicit none

contains

  subroutine actual_integrator_init()

    implicit none

    nseq = [2, 6, 10, 14, 22, 34, 50, 70]

    !$acc update device(nseq)

  end subroutine actual_integrator_init



  ! Main interface

  subroutine actual_integrator(state_in, state_out, dt, time)

    !$acc routine seq

    use rpar_indices
    use extern_probin_module, only: burner_verbose, &
                                    rtol_spec, rtol_temp, rtol_enuc, &
                                    atol_spec, atol_temp, atol_enuc, &
                                    burning_mode, retry_burn, &
                                    retry_burn_factor, retry_burn_max_change, &
                                    dT_crit
    use actual_rhs_module, only : update_unevolved_species

    implicit none

    ! Input arguments

    type (burn_t), intent(in   ) :: state_in
    type (burn_t), intent(inout) :: state_out
    real(dp_t),    intent(in   ) :: dt, time

    ! Local variables
    integer :: ierr

    real(kind=dp_t) :: atol(neqs), rtol(neqs)   ! input state, abs and rel tolerances
    real(kind=dp_t) :: t0, t1

    type (eos_t) :: eos_state_in, eos_state_temp
    type (bs_t) :: bs

    real(dp_t) :: retry_change_factor

    real(dp_t) :: ener_offset

    ! Set the tolerances.  We will be more relaxed on the temperature
    ! since it is only used in evaluating the rates.
    !
    ! **NOTE** if you reduce these tolerances, you probably will need
    ! to (a) decrease dT_crit, (b) increase the maximum number of
    ! steps allowed.


    atol(1:nspec_evolve) = atol_spec ! mass fractions
    atol(net_itemp)      = atol_temp ! temperature
    atol(net_ienuc)      = atol_enuc ! energy generated

    rtol(1:nspec_evolve) = rtol_spec ! mass fractions
    rtol(net_itemp)      = rtol_temp ! temperature
    rtol(net_ienuc)      = rtol_enuc ! energy generated

    ! Note that at present, we use a uniform error tolerance chosen
    ! to be the largest of the relative error tolerances for any
    ! equation. We may expand this capability in the future.

    bs % atol = atol
    bs % rtol = rtol

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

    if (burning_mode == 0 .or. burning_mode == 2) then
       bs % burn_s % self_heat = .false.

    else if (burning_mode == 1 .or. burning_mode == 3) then
       bs % burn_s % self_heat = .true.
    endif

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

    if (dT_crit < 1.0d19) then

       eos_state_temp = eos_state_in
       eos_state_temp % T = eos_state_in % T * (ONE + sqrt(epsilon(ONE)))

       call eos(eos_input_rt, eos_state_temp)

       bs % burn_s % dcvdt = (eos_state_temp % cv - eos_state_in % cv) / &
            (eos_state_temp % T - eos_state_in % T)

       bs % burn_s % dcpdt = (eos_state_temp % cp - eos_state_in % cp) / &
            (eos_state_temp % T - eos_state_in % T)
    endif

    ! Save the initial state.

    bs % upar(irp_y_init:irp_y_init + neqs - 1) = bs % y

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

       if (dT_crit < 1.0d19) then
          bs % burn_s % dcvdt = (eos_state_temp % cv - eos_state_in % cv) / &
               (eos_state_temp % T - eos_state_in % T)

          bs % burn_s % dcpdt = (eos_state_temp % cp - eos_state_in % cp) / &
               (eos_state_temp % T - eos_state_in % T)
       endif

       call ode(bs, t0, t1, maxval(rtol), ierr)

    endif

    ! If we still failed, print out the current state of the integration.

    if (ierr /= IERR_NONE) then

#ifndef ACC
       print *, 'ERROR: integration failed in net'
       print *, 'ierr = ', ierr
       print *, 'dt = ', dt
       print *, 'time start = ', time
       print *, 'time current = ', bs % t
       print *, 'dens = ', state_in % rho
       print *, 'temp start = ', state_in % T
       print *, 'xn start = ', state_in % xn
       print *, 'temp current = ', bs % y(net_itemp)
       print *, 'xn current = ', bs % y(1:nspec_evolve), &
            bs % upar(irp_nspec:irp_nspec+n_not_evolved-1)
       print *, 'energy generated = ', bs % y(net_ienuc) - ener_offset
#endif

       if (.not. retry_burn) then

#ifndef ACC
          call bl_error("ERROR in burner: integration failed")
#endif

       else

          print *, 'Retrying burn with looser tolerances'

          retry_change_factor = ONE

          do while (ierr /= IERR_NONE .and. &
                    retry_change_factor <= retry_burn_max_change)

             retry_change_factor = retry_change_factor * retry_burn_factor

             bs % atol = atol * retry_burn_factor
             bs % rtol = rtol * retry_burn_factor

             call eos_to_bs(eos_state_in, bs)

             bs % y(net_ienuc) = ener_offset

             ! redo the T_old, cv / cp extrapolation
             bs % burn_s % T_old = eos_state_in % T

             if (dT_crit < 1.0d19) then
                bs % burn_s % dcvdt = (eos_state_temp % cv - eos_state_in % cv) / &
                     (eos_state_temp % T - eos_state_in % T)
                
                bs % burn_s % dcpdt = (eos_state_temp % cp - eos_state_in % cp) / &
                     (eos_state_temp % T - eos_state_in % T)
             endif

             call ode(bs, t0, t1, maxval(rtol), ierr)

          enddo

          if (retry_change_factor > retry_burn_max_change .and. &
               ierr /= IERR_NONE) then

#ifndef ACC
             call bl_error("ERROR in burner: integration failed")
#endif

          endif

       endif

    endif

    ! Subtract the energy offset.

    bs % y(net_ienuc) = bs % y(net_ienuc) - ener_offset

    ! Store the final data, and then normalize abundances.
    call bs_to_burn(bs)

    state_out = bs % burn_s

    if (nspec_evolve < nspec) then
       call update_unevolved_species(state_out)
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

  end subroutine actual_integrator

end module actual_integrator_module
