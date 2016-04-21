module actual_burner_module

  use burn_type_module
  use network
  use actual_burner_data

contains

  subroutine actual_burner_init()

    use integrator_module, only: integrator_init

    implicit none

    reac_names(ir3a_)   = "3agc"   !     3 He4 --> C12
    reac_names(ircago_) = "cago"   ! C12 + He4 --> O16

    call integrator_init()

  end subroutine actual_burner_init



  subroutine actual_burner(state_in, state_out, dt, time)

    use integrator_module, only: integrator

    implicit none

    type (burn_t),    intent(in   ) :: state_in
    type (burn_t),    intent(inout) :: state_out
    double precision, intent(in   ) :: dt, time

    call integrator(state_in, state_out, dt, time)

  end subroutine actual_burner



  subroutine ener_gener_rate(dydt, enuc)

    implicit none

    double precision :: dydt(nspec), enuc

    enuc = sum(dydt(:) * ebin(:))

  end subroutine ener_gener_rate




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
    ymol = X / aion

    ! get the d/dT(dX/dt) info, dydt(dratesdT) gives us this
    call make_rates(temp, dens, rates, dratesdt)
    call screen(temp, dens, ymol, rates, dratesdt)
    call dydt(ymol, dratesdt, dXdotdT)

    ! calculate temperature sensitivity
    denucdT = - sum(dXdotdT*ebin)    

  end subroutine get_enuc_T_sensitivity

end module actual_burner_module
