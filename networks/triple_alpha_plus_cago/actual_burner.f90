module actual_burner_module

  use burn_type_module
  use network

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
    double precision, intent(in   ) :: dt, time

    call integrator(state_in, state_out, dt, time)

  end subroutine actual_burner



  subroutine get_enuc_T_sensitivity(dens, temp, X, denucdT)
    
    ! Calculate the energy generation rate's temperature sensitivity
    ! Used for diagnostic purposes only

    use rates_module
    use screen_module
    use dydt_module

    implicit none

    real(kind=dp_t), intent(IN   ) :: dens, temp, X(nspec)
    real(kind=dp_t), intent(  OUT) :: denucdT

    real(kind=dp_t) :: ymol(nspec)
    real(kind=dp_t) :: rates(nrates), dratesdt(nrates)
    real(kind=dp_t) :: dXdotdT(nspec)
    integer :: k

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
