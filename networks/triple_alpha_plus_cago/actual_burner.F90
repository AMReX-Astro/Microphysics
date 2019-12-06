module actual_burner_module

  use burn_type_module
  use network
  use microphysics_type_module

contains

  subroutine actual_burner_init()

    use integrator_module, only: integrator_init

    implicit none

    call integrator_init()

  end subroutine actual_burner_init



  subroutine actual_burner(state_in, state_out, dt, time)

    !$acc routine seq

    use integrator_module, only: integrator

    implicit none

    type (burn_t),    intent(in   ) :: state_in
    type (burn_t),    intent(inout) :: state_out
    real(rt), intent(in   ) :: dt, time

    !$gpu

    call integrator(state_in, state_out, dt, time)

  end subroutine actual_burner



  subroutine get_enuc_T_sensitivity(dens, temp, X, denucdT)

    ! Calculate the energy generation rate's temperature sensitivity
    ! Used for diagnostic purposes only

    use amrex_fort_module, only : rt => amrex_real
    use rates_module
    use screen_module
    use dydt_module

    implicit none

    real(kind=rt), intent(IN   ) :: dens, temp, X(nspec)
    real(kind=rt), intent(  OUT) :: denucdT

    real(kind=rt) :: ymol(nspec)
    real(kind=rt) :: rates(nrates), dratesdt(nrates)
    real(kind=rt) :: dXdotdT(nspec)
    integer :: k

    !$gpu

    ! calculate ymol
    ymol = X * aion_inv

    ! get the d/dT(dX/dt) info, dydt(dratesdT) gives us this
    call make_rates(temp, dens, rates, dratesdt)
    call screen(temp, dens, ymol, rates, dratesdt)
    call dydt(ymol, dratesdt, dXdotdT)

    ! calculate temperature sensitivity
    denucdT = - sum(dXdotdT*ebin)

  end subroutine get_enuc_T_sensitivity

end module actual_burner_module
