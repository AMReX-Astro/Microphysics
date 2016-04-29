! Common variables and routines for burners
! that use BS for their integration.

module actual_integrator_module

  use eos_module
  use network
  use rpar_indices
  use bs_convert_module
  use burn_type_module
  use bl_types
  use stiff_ode
  use bs_type_module

  implicit none

contains

  subroutine actual_integrator_init()

    implicit none

  end subroutine actual_integrator_init



  ! Main interface

  subroutine actual_integrator(state_in, state_out, dt, time)

    use rpar_indices
    use extern_probin_module, only: jacobian, burner_verbose, &
                                    rtol_spec, rtol_temp, rtol_enuc, &
                                    atol_spec, atol_temp, atol_enuc, &
                                    burning_mode, retry_burn, &
                                    retry_burn_factor, retry_burn_max_change, &
                                    call_eos_in_rhs, dT_crit
    use integration_data, only: ener_scale

    implicit none

    ! Input arguments

    type (burn_t),       intent(in   ) :: state_in
    type (burn_t),       intent(inout) :: state_out
    double precision,    intent(in   ) :: dt, time

    logical, parameter :: RESET = .true.  !.true. means we want to initialize the bdf_ts object
    logical, parameter :: REUSE = .false. !.false. means don't reuse the Jacobian
    real(kind=dp_t), parameter :: DT0 = 1.0d-9 !Initial dt to be used in getting from 
                                               !t to tout.  Also arbitrary,
                                               !multiple values should be
                                               !explored.

    ! Local variables
    integer :: n, i, j, ierr

    real(kind=dp_t) :: atol(neqs), rtol(neqs)   ! input state, abs and rel tolerances
    real(kind=dp_t) :: t0, t1

    type (eos_t) :: eos_state_in, eos_state_out, eos_state_temp
    type (bs_t) :: bs

    double precision :: sum
    double precision :: retry_change_factor

    ! Set the tolerances.  We will be more relaxed on the temperature
    ! since it is only used in evaluating the rates.
    !
    ! **NOTE** if you reduce these tolerances, you probably will need
    ! to (a) decrease dT_crit, (b) increase the maximum number of
    ! steps allowed.

    atol(1:nspec)   = atol_spec ! mass fractions
    atol(net_itemp) = atol_temp ! temperature
    atol(net_ienuc) = atol_enuc ! energy generated

    rtol(1:nspec)   = rtol_spec ! mass fractions
    rtol(net_itemp) = rtol_temp ! temperature
    rtol(net_ienuc) = rtol_enuc ! energy generated

    bs % atol = atol
    bs % rtol = rtol

    ! Initialize the integration time.

    t0 = ZERO
    t1 = t0 + dt

    ! Convert our input burn state into an EOS type.

    call burn_to_eos(state_in, eos_state_in)

    ! We assume that the valid quantities coming in are (rho, e); do an EOS call
    ! to make sure all other variables are consistent.

    call eos(eos_input_burn, eos_state_in)

    ! Send this data back to the burn state in case the energy changed
    ! due to a reset/flooring.

    call eos_to_burn(eos_state_in, state_in)

    ! Convert the EOS state data into the form BS expects.

    call eos_to_bs(eos_state_in, bs)

    bs % y(net_ienuc) = ZERO

    ! Pass through whether we are doing self-heating.

    if (burning_mode == 0 .or. burning_mode == 2) then
       bs % upar(irp_self_heat) = -ONE
    else if (burning_mode == 1) then
       bs % upar(irp_self_heat) = ONE
    else
       call bl_error("Error: unknown burning_mode in actual_integrator.f90.")
    endif

    ! If we are using the dT_crit functionality and therefore doing a linear
    ! interpolation of the specific heat in between EOS calls, do a second
    ! EOS call here to establish an initial slope.

    bs % upar(irp_Told) = eos_state_in % T

    if (dT_crit < 1.0d19) then

       eos_state_temp = eos_state_in
       eos_state_temp % T = eos_state_in % T * (ONE + sqrt(epsilon(ONE)))

       call eos(eos_input_rt, eos_state_temp)

       bs % upar(irp_dcvdt) = (eos_state_temp % cv - eos_state_in % cv) / (eos_state_temp % T - eos_state_in % T)
       bs % upar(irp_dcpdt) = (eos_state_temp % cp - eos_state_in % cp) / (eos_state_temp % T - eos_state_in % T)

    endif

    ! Call the integration routine.

    call ode(bs, t0, t1, maxval(rtol), ierr)

    ! If we are using hybrid burning and the energy release was negative (or we failed),
    ! re-run this in self-heating mode.

    if ( burning_mode == 2 .and. (bs % y(net_ienuc) < ZERO .or. ierr /= IERR_NONE) ) then

       bs % upar(irp_self_heat) = ONE

       call eos_to_bs(eos_state_in, bs)

       bs % y(net_ienuc) = ZERO

       call ode(bs, t0, t1, maxval(rtol), ierr)

    endif

    ! If we still failed, print out the current state of the integration.

    if (ierr /= IERR_NONE) then
       print *, 'ERROR: integration failed in net'
       print *, 'ierr = ', ierr
       print *, 'time = ', bs % t
       print *, 'dens = ', state_in % rho
       print *, 'temp start = ', state_in % T
       print *, 'xn start = ', state_in % xn
       print *, 'temp current = ', bs % y(net_itemp) * temp_scale
       print *, 'xn current = ', bs % y(1:nspec_evolve) * aion(1:nspec_evolve), &
            bs % upar(irp_nspec:irp_nspec+nspec-nspec_evolve-1) * aion(nspec_evolve+1:)
       print *, 'energy generated = ', bs % y(net_ienuc) * ener_scale

       if (.not. retry_burn) then

          call bl_error("ERROR in burner: integration failed")

       else

          print *, 'Retrying burn with looser tolerances'

          retry_change_factor = ONE

          do while (ierr /= IERR_NONE .and. retry_change_factor <= retry_burn_max_change)

             retry_change_factor = retry_change_factor * retry_burn_factor

             bs % atol = atol * retry_burn_factor
             bs % rtol = rtol * retry_burn_factor

             call eos_to_bs(eos_state_in, bs)

             bs % y(net_ienuc) = ZERO

             call ode(bs, t0, t1, maxval(rtol), ierr)

          enddo

          if (retry_change_factor > retry_burn_max_change .and. ierr /= IERR_NONE) then

             call bl_error("ERROR in burner: integration failed")

          endif

       endif

    endif

    ! Store the final data.

    call bs_to_eos(eos_state_out, bs)

    call normalize_abundances(eos_state_out)

    ! Energy was integrated in the system -- we use this integrated
    ! energy which contains both the reaction energy release and
    ! neutrino losses. The final energy is the initial energy
    ! plus this energy release. Note that we get a new temperature too,
    ! but we will discard it and call the EOS to get a final temperature
    ! consistent with this new energy.

    eos_state_out % e = eos_state_in % e + bs % y(net_ienuc) * ener_scale

    eos_state_out % reset = .true.

    call eos(eos_input_re, eos_state_out)

    call eos_to_burn(eos_state_out, state_out)

    if (burner_verbose) then

       ! Print out some integration statistics, if desired.

       print *, 'integration summary: '
       print *, 'dens: ', state_out % rho
       print *, ' temp: ', state_out % T
       print *, ' energy released: ', bs % y(net_ienuc) * ener_scale
       print *, 'number of steps taken: ', bs % n
!       print *, 'number of f evaluations: ', bs % nfe

    endif

  end subroutine actual_integrator

end module actual_integrator_module
