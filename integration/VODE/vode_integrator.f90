! Common variables and routines for burners
! that use VODE for their integration.

module vode_integrator_module

  use eos_module, only : eos
  use eos_type_module, only: eos_t, eos_input_rt
  use network
  use rpar_indices
  use vode_type_module
  use burn_type_module
  use amrex_error_module
  use amrex_fort_module, only : rt => amrex_real

  implicit none

  ! Our problem is stiff, so tell ODEPACK that. 21 means stiff, jacobian
  ! function is supplied; 22 means stiff, figure out my jacobian through
  ! differencing.

  integer, parameter :: MF_ANALYTIC_JAC = 21, MF_NUMERICAL_JAC = 22

  ! Tolerance parameters:
  !
  !  itol specifies whether to use an single absolute tolerance for
  !  all variables (1), or to pass an array of absolute tolerances, one
  !  for each variable with a scalar relative tol (2), a scalar absolute
  !  and array of relative tolerances (3), or arrays for both (4).
  !
  !  The error is determined as e(i) = rtol*abs(y(i)) + atol, and must
  !  be > 0.  Since we have some compositions that may be 0 initially,
  !  we will specify both an absolute and a relative tolerance.
  !
  ! We will use arrays for both the absolute and relative tolerances,
  ! since we want to be easier on the temperature than the species.

  integer, parameter :: ITOL = 4

  ! We want to do a normal computation, and get the output values of y(t)
  ! after stepping though dt.

  integer, PARAMETER :: ITASK = 1

  ! We will override the maximum number of steps, so turn on the
  ! optional arguments flag.

  integer, parameter :: IOPT = 1

  ! Declare a real work array of size 22 + 9*NEQ + 2*NEQ**2 and an
  ! integer work array of size 30 + NEQ. These are VODE constants
  ! that depend on the integration mode we're using -- see dvode.f.

  integer, parameter :: LRW = 22 + 9*neqs + 2*neqs**2
  integer, parameter :: LIW = 30 + neqs

  ! Maximum number of timesteps to take.

  integer, save :: max_steps

contains

  subroutine vode_integrator_init()

    use extern_probin_module, only: ode_max_steps

    implicit none

    max_steps = ode_max_steps

  end subroutine vode_integrator_init



  ! Main interface

  subroutine vode_integrator(state_in, state_out, dt, time, status)

    use rpar_indices
    use extern_probin_module, only: jacobian, burner_verbose, &
                                    burning_mode, burning_mode_factor, dT_crit
    use actual_rhs_module, only : update_unevolved_species
    use amrex_constants_module, only : ZERO, ONE
    use integration_data, only: integration_status_t

    implicit none

    ! Input arguments

    type (burn_t), intent(in   ) :: state_in
    type (burn_t), intent(inout) :: state_out
    real(rt),    intent(in   ) :: dt, time
    type (integration_status_t), intent(inout) :: status

    ! Local variables

    real(rt) :: local_time
    type (eos_t) :: eos_state_in, eos_state_temp

    ! Work arrays

    real(rt) :: y(neqs)
    real(rt) :: atol(neqs), rtol(neqs)
    real(rt) :: rwork(LRW)
    integer    :: iwork(LIW)
    real(rt) :: rpar(n_rpar_comps)

    integer :: MF_JAC

    ! istate determines the state of the calculation.  A value of 1 meeans
    ! this is the first call to the problem -- this is what we will want.

    integer :: istate

    integer :: ipar

    real(rt) :: ener_offset
    real(rt) :: edot, t_enuc, t_sound, limit_factor

    logical :: integration_failed
    real(rt), parameter :: failure_tolerance = 1.d-2

    EXTERNAL jac, f_rhs

    if (jacobian == 1) then ! Analytical
       MF_JAC = MF_ANALYTIC_JAC
    else if (jacobian == 2) then ! Numerical
       MF_JAC = MF_NUMERICAL_JAC
    else
       call amrex_error("Error: unknown Jacobian mode in vode_integrator.f90.")
    endif

    integration_failed = .false.

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

    ! We want VODE to re-initialize each time we call it.

    istate = 1

    ! Initialize work arrays to zero.

    rwork(:) = ZERO
    iwork(:) = 0

    ! Set the maximum number of steps allowed (the VODE default is 500).

    iwork(6) = max_steps

    ! Disable printing of messages about T + H == T unless we are in verbose mode.

    if (burner_verbose) then
       iwork(7) = 1
    else
       iwork(7) = 0
    endif

    ! Start off by assuming a successful burn.

    state_out % success = .true.

    ! Initialize the integration time.

    local_time = ZERO

    ! Convert our input burn state into an EOS type.

    call burn_to_eos(state_in, eos_state_in)

    ! We assume that (rho, T) coming in are valid, do an EOS call
    ! to fill the rest of the thermodynamic variables.

    call eos(eos_input_rt, eos_state_in)

    ! Convert the EOS state data into the form VODE expects.

    call eos_to_vode(eos_state_in, y, rpar)

    rpar(irp_i) = state_in % i
    rpar(irp_j) = state_in % j
    rpar(irp_k) = state_in % k

    ener_offset = eos_state_in % e

    y(net_ienuc) = ener_offset

    ! Pass through whether we are doing self-heating.

    if (burning_mode == 0 .or. burning_mode == 2) then
       rpar(irp_self_heat) = -ONE
    else if (burning_mode == 1 .or. burning_mode == 3) then
       rpar(irp_self_heat) = ONE
    else
       call amrex_error("Error: unknown burning_mode in vode_integrator.f90.")
    endif

    ! Copy in the zone size.

    rpar(irp_dx) = state_in % dx

    ! Set the sound crossing time.

    rpar(irp_t_sound) = state_in % dx / eos_state_in % cs

    ! Set the time offset -- this converts between the local integration 
    ! time and the simulation time

    rpar(irp_t0) = time

    ! If we are using the dT_crit functionality and therefore doing a linear
    ! interpolation of the specific heat in between EOS calls, do a second
    ! EOS call here to establish an initial slope.

    rpar(irp_Told) = eos_state_in % T

    if (dT_crit < 1.0d19) then

       eos_state_temp = eos_state_in
       eos_state_temp % T = eos_state_in % T * (ONE + sqrt(epsilon(ONE)))

       call eos(eos_input_rt, eos_state_temp)

       rpar(irp_dcvdt) = (eos_state_temp % cv - eos_state_in % cv) / (eos_state_temp % T - eos_state_in % T)
       rpar(irp_dcpdt) = (eos_state_temp % cp - eos_state_in % cp) / (eos_state_temp % T - eos_state_in % T)

    endif

    ! Save the initial state.

    rpar(irp_y_init:irp_y_init + neqs - 1) = y

    ! Call the integration routine.

    call dvode(f_rhs, neqs, y, local_time, local_time + dt, &
               ITOL, rtol, atol, ITASK, &
               istate, IOPT, rwork, LRW, iwork, LIW, jac, MF_JAC, rpar, ipar)

    ! If we are using hybrid burning and the energy release was negative (or we failed),
    ! re-run this in self-heating mode.

    if ( burning_mode == 2 .and. &
         (y(net_ienuc) - ener_offset < ZERO .or. &
          istate < 0) ) then

       rpar(irp_self_heat) = ONE

       istate = 1

       rwork(:) = ZERO
       iwork(:) = 0

       iwork(6) = 150000

       local_time = ZERO

       call eos_to_vode(eos_state_in, y, rpar)

       rpar(irp_Told) = eos_state_in % T

       if (dT_crit < 1.0d19) then

          rpar(irp_dcvdt) = (eos_state_temp % cv - eos_state_in % cv) / (eos_state_temp % T - eos_state_in % T)
          rpar(irp_dcpdt) = (eos_state_temp % cp - eos_state_in % cp) / (eos_state_temp % T - eos_state_in % T)

       endif

       y(net_ienuc) = ener_offset

       call dvode(f_rhs, neqs, y, local_time, local_time + dt, &
                  ITOL, rtol, atol, ITASK, &
                  istate, IOPT, rwork, LRW, iwork, LIW, jac, MF_JAC, rpar, ipar)

    endif

    ! If we still failed, print out the current state of the integration.
    ! VODE does not always fail even though it can lead to unphysical states,
    ! so add some sanity checks that trigger a retry even if VODE thinks
    ! the integration was successful.

    if (istate < 0) then
       integration_failed = .true.
    end if

    if (y(net_itemp) < ZERO) then
       integration_failed = .true.
    end if

    if (any(y(1:nspec_evolve) < -failure_tolerance)) then
       integration_failed = .true.
    end if

    if (any(y(1:nspec_evolve) > 1.d0 + failure_tolerance)) then
       integration_failed = .true.
    end if

    if (integration_failed) then
       print *, 'ERROR: integration failed in net'
       print *, 'istate = ', istate
       print *, 'time = ', local_time
       print *, 'dens = ', state_in % rho
       print *, 'temp start = ', state_in % T
       print *, 'xn start = ', state_in % xn
       print *, 'temp current = ', y(net_itemp)
       print *, 'xn current = ', y(1:nspec_evolve), rpar(irp_nspec:irp_nspec+n_not_evolved-1)
       print *, 'energy generated = ', y(net_ienuc) - ener_offset

       state_out % success = .false.
       return
    endif

    ! Subtract the energy offset
    y(net_ienuc) = y(net_ienuc) - ener_offset

    ! Store the final data, and then normalize abundances.
    call vode_to_burn(y, rpar, state_out)

    ! Get the number of RHS calls and jac evaluations from the VODE
    ! work arrays
    state_out % n_rhs = iwork(12)
    state_out % n_jac = iwork(13)

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

    ! set the integration time for any diagnostics
    state_out % time = time + dt

    if (burner_verbose) then

       ! Print out some integration statistics, if desired.

       print *, 'integration summary: '
       print *, 'dens: ', state_out % rho, ' temp: ', state_out % T, &
                ' energy released: ', state_out % e - state_in % e
       print *, 'number of steps taken: ', iwork(11)
       print *, 'number of f evaluations: ', iwork(12)

    endif

  end subroutine vode_integrator

end module vode_integrator_module
