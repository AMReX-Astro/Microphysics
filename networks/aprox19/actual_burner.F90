module actual_burner_module

  use amrex_constants_module
  use amrex_error_module
  use eos_module
  use eos_type_module
  use network
  use burn_type_module
#ifdef NSE
  use nse_module
  use nse_check_module
#endif
  use amrex_fort_module, only : rt => amrex_real
  implicit none

contains

  subroutine actual_burner(state_in, state_out, dt, time)

    use integrator_module, only: integrator
    use extern_probin_module, only : eta
    implicit none

    type (burn_t),       intent(in   ) :: state_in
    type (burn_t),       intent(inout) :: state_out
    real(rt)        ,    intent(in   ) :: dt, time

#ifdef NSE
    integer :: nse_check
    real(rt) :: deltaq, enuc
    real(rt) :: dq_out, abar_out, dyedt
    type(eos_t) :: eos_state
#endif

    !$gpu

#ifndef NSE
    call integrator(state_in, state_out, dt, time)

#else

    call in_nse(state_in, nse_check)

    if (nse_check == 0) then

       ! integrate the aprox19 network

       call integrator(state_in, state_out, dt, time)

       ! update the auxiliary state
       state_out % aux(iye) = sum(state_out % xn(:) * zion(:) * aion_inv(:))
       state_out % aux(iabar) = ONE / (sum(state_out % xn(:) * aion_inv(:)))
       state_out % aux(ibea) = (sum(state_out % xn(:) * bion(:) * aion_inv(:)))

    else

       ! use the NSE table
       call nse_interp(state_in % T, state_in % rho, state_in % aux(iye), &
                       abar_out, dq_out, dyedt, state_out % xn(:))

       ! update Ye
       state_out % aux(iye) = state_in % aux(iye) + dt * dyedt

       ! now get the composition from the table using the updated Ye
       call nse_interp(state_in % T, state_in % rho, state_out % aux(iye), &
                       abar_out, dq_out, dyedt, state_out % xn(:))

       ! this is MeV / nucleon
       deltaq = dq_out - state_in % aux(ibea)

       ! under-relaxation / inertia (see Ma et el. 2013)
       deltaq = eta * deltaq

       ! convert the energy to erg / g
       enuc = deltaq * mev2erg * avo

       ! now update the temperature based on the energy deposition
       eos_state % rho = state_in % rho
       eos_state % T = state_in % T
       eos_state % e = state_in % e + enuc
       eos_state % xn(:) = state_out % xn(:)
       eos_state % aux(iye) = state_out % aux(iye)
       eos_state % aux(iabar) = abar_out
       eos_state % aux(ibea) = state_in % aux(ibea) + deltaq

       call eos(eos_input_re, eos_state)

       ! now call the table one last time, with the updated T
       call nse_interp(eos_state % T, state_in % rho, eos_state % aux(iye), &
                       abar_out, dq_out, dyedt, state_out % xn(:))

       deltaq = eta * deltaq
       enuc = deltaq * mev2erg * avo

       state_out % aux(ibea) = state_in % aux(ibea) + deltaq
       state_out % aux(iabar) = abar_out

       state_out % e = enuc + state_in % e

       state_out % success = .true.
       state_out % n_rhs = 0
       state_out % n_jac = 0

       state_out % time = time + dt

    end if

#endif

  end subroutine actual_burner



  subroutine actual_burner_init()

    use integrator_module, only: integrator_init

    implicit none

    call integrator_init()

  end subroutine actual_burner_init

end module actual_burner_module
