! This is the equation of state for a polytropic fluid:
! P = K rho^gamma
!
! The internal energy and pressure are related via a gamma law:
!
! P = (gamma - 1) rho e
!
! gamma and K are fixed quantities for the run, and must either be
! supplied by the user or selected from a list of available options.
! Currently, we have fully degenerate ionized gases (both relativistic
! and non-relativistic), where the pressure is supplied by electrons.
!
! Note that here we define the mean number of electrons per ion as:
!
!   1/mu_e = sum_k { X_k Z_k / A_k }
!
! This is assumed to be constant for the degenerate gases.

module specific_eos_module

  use bl_types
  use bl_space
  use bl_error_module
  use bl_constants_module
  use network, only: nspec, aion, zion
  use eos_type_module
  use eos_data_module

  implicit none

  double precision :: gamma_const, gm1, K_const
  double precision :: mu_e
  integer          :: polytrope

contains

  subroutine specific_eos_init

    use extern_probin_module, only: polytrope_gamma, polytrope_K, polytrope_type, polytrope_mu_e

    implicit none
 
    ! Available pre-defined polytrope options:

    ! 1: Non-relativistic, fully degenerate electron gas
    ! 2: Relativistic, fully degenerate electron gas 

    if (polytrope_type > 0) then
      mu_e = polytrope_mu_e

      polytrope = polytrope_type
      if (polytrope .eq. 1) then
        gamma_const = FIVE3RD
        K_const     = 9.9154d12 ! (3 / pi)^(2/3) * h^2 / (20 * m_e * m_p^(5/3))
        K_const     = K_const / mu_e**gamma_const
      elseif (polytrope .eq. 2) then
        gamma_const = FOUR3RD
        K_const     = 1.2316d15 ! (3 / pi)^(1/3) * h c / (8 * m_p^(4/3))
        K_const     = K_const / mu_e**gamma_const
      else
        call bl_error('EOS: Polytrope type currently not defined')
      endif
    elseif (polytrope_gamma .gt. ZERO .and. polytrope_K .gt. ZERO) then
      gamma_const = polytrope_gamma
      K_const     = polytrope_K
      mu_e        = TWO ! This will not be used
    else
      call bl_error('EOS: Neither polytrope type nor both gamma and K are defined')
    endif

    gm1 = gamma_const - ONE

    initialized = .true.
 
  end subroutine specific_eos_init



  !---------------------------------------------------------------------------
  ! Public interfaces 
  !---------------------------------------------------------------------------

  subroutine eos_get_polytrope_parameters(polytrope_out,gamma_out,K_out,mu_e_out)

    integer,          intent(out) :: polytrope_out
    double precision, intent(out) :: gamma_out, K_out, mu_e_out

    polytrope_out = polytrope
    gamma_out     = gamma_const
    K_out         = K_const
    mu_e_out      = mu_e

  end subroutine eos_get_polytrope_parameters

  subroutine eos_set_polytrope_parameters(polytrope_in,gamma_in,K_in,mu_e_in)

    integer,          intent(in) :: polytrope_in
    double precision, intent(in) :: gamma_in, K_in, mu_e_in

    polytrope   = polytrope_in
    gamma_const = gamma_in
    K_const     = K_in
    mu_e        = mu_e_in

  end subroutine eos_set_polytrope_parameters



  !---------------------------------------------------------------------------
  ! The main interface
  !---------------------------------------------------------------------------
  subroutine specific_eos(eosfail, state, input)

    implicit none

    logical,           intent(in   ) :: eosfail
    type (eos_t),      intent(inout) :: state(:)
    integer,           intent(in   ) :: input

    ! Local variables
    double precision :: dens, temp, enth, pres, eint, entr

    integer :: j, N

    N = size(state)

    if (.not. initialized) call bl_error('EOS: not initialized')

    do j = 1, N

       dens = state(j) % rho
       temp = state(j) % T
       pres = state(j) % p
       enth = state(j) % h
       eint = state(j) % e
       entr = state(j) % s

       select case (input)

       !-------------------------------------------------------------------------
       ! Now do the calculations. In every case,
       ! make sure we have pressure, density, energy, and enthalpy.
       ! Relevant equations:
       ! h   = e + p / rho = (p / rho) * gamma / (gamma - 1) = e * gamma
       ! p   = K * (rho ** gamma) = (gamma - 1) * rho * e
       ! rho = (p / K)**(1 / gamma)
       ! e   = h - p / rho = (p / rho) / (gamma - 1)         = h / gamma
       !-------------------------------------------------------------------------

       case (eos_input_rh)

          ! dens, enthalpy, and xmass are inputs

          ! Solve for the pressure and energy:

          pres = enth * dens * gm1 / gamma_const
          eint = enth / gamma_const


       case (eos_input_rt)

          ! dens, temp, and xmass are inputs

          ! Solve for the pressure, energy and enthalpy:

          pres = K_const * dens**gamma_const
          enth = pres / dens * gamma_const / gm1
          eint = enth / gamma_const


       case (eos_input_tp)

          ! temp, pres, and xmass are inputs

          ! Solve for the density, energy and enthalpy:

          dens = (pres / K_const)**(ONE / gamma_const)
          enth = pres / dens * gamma_const / gm1
          eint = enth / gamma_const


       case (eos_input_rp)

          ! dens, pres, and xmass are inputs

          ! Solve for the enthalpy and energy:

          enth = (pres / dens) * gamma_const / gm1
          eint = enth / gamma_const


       case (eos_input_re)

          ! dens, energy, and xmass are inputs

          ! Solve for the pressure and enthalpy:

          pres = gm1 * dens * eint


       case (eos_input_ps)

          ! pressure, entropy and xmass are inputs

          ! Solve for the density, energy and enthalpy:

          dens = (pres / K_const)**(ONE / gamma_const)
          enth = pres / dens * gamma_const / gm1
          eint = enth / gamma_const



       case (eos_input_ph)

          ! pressure, enthalpy and xmass are inputs

          ! Solve for the density and energy:

          eint = enth / gamma_const
          dens = (pres / K_const)**(ONE / gamma_const)



       case (eos_input_th)

          ! temperature, enthalpy and xmass are inputs

          ! Solve for the density, energy and pressure:

          eint = enth / gamma_const
          dens = (gm1 / gamma_const * enth / K_const)**(ONE / gm1)
          pres = gm1 * dens * eint



       case default

          call bl_error('EOS: invalid input.')

       end select

       !-------------------------------------------------------------------------
       ! Now we have all relevant quantities, regardless of the inputs.
       !-------------------------------------------------------------------------

       state(j) % T   = temp
       state(j) % rho = dens
       state(j) % h   = enth
       state(j) % s   = entr
       state(j) % e   = eint
       state(j) % p   = pres

       ! Compute the thermodynamic derivatives and specific heats 
       state(j) % dPdT = ZERO
       state(j) % dPdr = gamma_const * pres / dens
       state(j) % dedT = ZERO
       state(j) % dedr = pres / (dens * dens)
       state(j) % dsdT = ZERO
       state(j) % dsdr = ZERO
       state(j) % dhdT = ZERO
       state(j) % dhdr = state(j) % dedr + gm1 * pres / dens**2

       state(j) % cv = state(j) % dedT
       state(j) % cp = gamma_const * state(j) % cv

       state(j) % gam1 = gamma_const

       ! Compute dPdX, dedX, dhdX.

       state(j) % dpdA = - state(j) % p / state(j) % abar
       state(j) % dpdZ =   state(j) % p / (ONE + state(j) % zbar)

       state(j) % dedA = - state(j) % e / state(j) % abar
       state(j) % dedZ =   state(j) % e / (ONE + state(j) % zbar)

       ! sound speed
       state(j) % cs = sqrt(gamma_const*pres/dens)

    enddo

  end subroutine specific_eos

end module specific_eos_module
