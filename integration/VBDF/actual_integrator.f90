! Common variables and routines for burners
! that use VBDF for their integration.

module actual_integrator_module

  use eos_module
  use network
  use rpar_indices
  use vbdf_convert_module
  use burn_type_module
  use bl_types
  use bdf_type_module
  use bdf

  implicit none

contains

  subroutine actual_integrator_init()

    implicit none

    call init_rpar_indices()

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

    implicit none

    ! Input arguments

    type (burn_t),       intent(in   ) :: state_in
    type (burn_t),       intent(inout) :: state_out
    double precision,    intent(in   ) :: dt, time

    integer, parameter :: MAX_ORDER = 5   !This is arbitrary, should investigate other values
    logical, parameter :: RESET = .true.  !.true. means we want to initialize the bdf_ts object
    logical, parameter :: REUSE = .false. !.false. means don't reuse the Jacobian
    real(kind=dp_t), parameter :: DT0 = 1.0d-9 !Initial dt to be used in getting from 
                                               !t to tout.  Also arbitrary,
                                               !multiple values should be
                                               !explored.

    !!Local variables
    integer, parameter :: bdf_npt = 1 ! Number of points for BDF to operate on
    integer :: n, i, j, ierr

    real(kind=dp_t), dimension(neqs) :: atol, rtol   ! input state, abs and rel tolerances
    real(kind=dp_t), allocatable :: upar(:,:), y0(:,:), y1(:,:)
    real(kind=dp_t) :: t0, t1, enuc, dX

    type (eos_t)  :: eos_state_in, eos_state_out, eos_state_temp
    type (bdf_ts) :: ts

    double precision :: sum
    double precision :: retry_change_factor

    allocate(upar(n_rpar_comps, bdf_npt))
    allocate(y0(neqs, bdf_npt), y1(neqs, bdf_npt))

    ! Build the bdf_ts time-stepper object
    call bdf_ts_build(ts, neqs, bdf_npt, rtol, atol, MAX_ORDER, upar)

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

    ts % atol = atol
    ts % rtol = rtol

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

    ! Convert the EOS state data into the form VBDF expects.

    call eos_to_vbdf(eos_state_in, ts)

    ts % y(net_ienuc,1) = ZERO

    ! Pass through whether we are doing self-heating.

    if (burning_mode == 0 .or. burning_mode == 2) then
       ts % upar(irp_self_heat,1) = -ONE
    else if (burning_mode == 1) then
       ts % upar(irp_self_heat,1) = ONE
    else
       call bl_error("Error: unknown burning_mode in actual_integrator.f90.")
    endif

    ! If we are using the dT_crit functionality and therefore doing a linear
    ! interpolation of the specific heat in between EOS calls, do a second
    ! EOS call here to establish an initial slope.

    ts % upar(irp_Told,1) = eos_state_in % T

    if (dT_crit < 1.0d19) then

       eos_state_temp = eos_state_in
       eos_state_temp % T = eos_state_in % T * (ONE + sqrt(epsilon(ONE)))

       call eos(eos_input_rt, eos_state_temp)

       ts % upar(irp_dcvdt,1) = (eos_state_temp % cv - eos_state_in % cv) / (eos_state_temp % T - eos_state_in % T)
       ts % upar(irp_dcpdt,1) = (eos_state_temp % cp - eos_state_in % cp) / (eos_state_temp % T - eos_state_in % T)

    endif

    ! Call the integration routine.

    do n = 1, neqs
       y0(n,1) = ts % y(n,1)
    end do
    call bdf_advance(ts, neqs, bdf_npt, y0, t0, y1, t1, &
                     DT0, RESET, REUSE, ierr, .true.)
    do n = 1, neqs
       ts % y(n,1) = y1(n,1)
    end do

    ! If we are using hybrid burning and the energy release was negative (or we failed),
    ! re-run this in self-heating mode.

    if ( burning_mode == 2 .and. (ts % y(net_ienuc,1) < ZERO .or. ierr /= BDF_ERR_SUCCESS) ) then

       ts % upar(irp_self_heat,1) = ONE

       call eos_to_vbdf(eos_state_in, ts)

       ts % y(net_ienuc,1) = ZERO

       do n = 1, neqs
          y0(n,1) = ts % y(n,1)
       end do
       call bdf_advance(ts, neqs, bdf_npt, y0, t0, y1, t1, &
                        DT0, RESET, REUSE, ierr, .true.)
       do n = 1, neqs
          ts % y(n,1) = y1(n,1)
       end do

    endif

    ! If we still failed, print out the current state of the integration.

    if (ierr /= BDF_ERR_SUCCESS) then
       print *, 'ERROR: integration failed in net'
       print *, 'ierr = ', ierr
       print *, 'time = ', ts % t
       print *, 'dens = ', state_in % rho
       print *, 'temp start = ', state_in % T
       print *, 'xn start = ', state_in % xn
       print *, 'temp current = ', ts % y(net_itemp,1) * temp_scale
       print *, 'xn current = ', ts % y(1:nspec,1)
       print *, 'energy generated = ', ts % y(net_ienuc,1)

       if (.not. retry_burn) then

          call bl_error("ERROR in burner: integration failed")

       else

          print *, 'Retrying burn with looser tolerances'

          retry_change_factor = ONE

          do while (ierr /= BDF_ERR_SUCCESS .and. retry_change_factor <= retry_burn_max_change)

             retry_change_factor = retry_change_factor * retry_burn_factor

             ts % atol = atol * retry_burn_factor
             ts % rtol = rtol * retry_burn_factor

             call eos_to_vbdf(eos_state_in, ts)

             ts % y(net_ienuc,1) = ZERO

             do n = 1, neqs
                y0(n,1) = ts % y(n,1)
             end do
             call bdf_advance(ts, neqs, bdf_npt, y0, t0, y1, t1, &
                              DT0, RESET, REUSE, ierr, .true.)
             do n = 1, neqs
                ts % y(n,1) = y1(n,1)
             end do

          enddo

          if (retry_change_factor > retry_burn_max_change .and. ierr /= BDF_ERR_SUCCESS) then

             call bl_error("ERROR in burner: integration failed")

          endif

       endif

    endif

    ! Store the final data.

    call vbdf_to_eos(eos_state_out, ts)

    call normalize_abundances(eos_state_out)

    ! Energy was integrated in the system -- we use this integrated
    ! energy which contains both the reaction energy release and
    ! neutrino losses. The final energy is the initial energy
    ! plus this energy release. Note that we get a new temperature too,
    ! but we will discard it and call the EOS to get a final temperature
    ! consistent with this new energy.

    eos_state_out % e = eos_state_in % e + ts % y(net_ienuc,1)

    eos_state_out % reset = .true.

    call eos(eos_input_re, eos_state_out)

    call eos_to_burn(eos_state_out, state_out)

    call bdf_ts_destroy(ts)

    if (burner_verbose) then

       ! Print out some integration statistics, if desired.

       print *, 'integration summary: '
       print *, 'dens: ', state_out % rho, ' temp: ', state_out % T, ' energy released: ', ts % y(net_ienuc,1)
       print *, 'number of steps taken: ', ts % n
       print *, 'number of f evaluations: ', ts % nfe

    endif

  end subroutine actual_integrator

end module actual_integrator_module
