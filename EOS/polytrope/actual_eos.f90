! This is the equation of state for a polytropic fluid:
! P = K rho^gamma
!
! The internal energy is given by a gamma law:
!
! e = (P / rho) * (1 / (gamma - 1))
!
! Unlike the gamma law EOS, e is always a dependent variable
! that is directly determined by the fluid density. This guarantees
! that the fluid always obeys the polytropic relationship.
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

module actual_eos_module

  use bl_types
  use bl_space
  use bl_error_module
  use bl_constants_module
  use network, only: nspec, aion, zion
  use eos_type_module

  implicit none

  character (len=64), public :: eos_name = "polytrope"
  
  double precision, save :: gamma_const, K_const
  double precision, save :: mu_e
  integer         , save :: polytrope

  double precision, save :: gm1, polytrope_index

contains

  subroutine actual_eos_init

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

    polytrope_index = ONE / (gamma_const - ONE)

  end subroutine actual_eos_init



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
  subroutine actual_eos(input, state)

    implicit none

    integer,      intent(in   ) :: input
    type (eos_t), intent(inout) :: state

    ! Local variables
    double precision :: dens, temp, enth, pres, eint, entr

    dens = state % rho
    temp = state % T
    pres = state % p
    enth = state % h
    eint = state % e
    entr = state % s

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

       pres = K_const * dens**gamma_const
       enth = eint * gamma_const

    case (eos_input_ps)

       ! pressure, entropy and xmass are inputs

       ! Solve for the density, energy and enthalpy:

       dens = (pres / K_const)**(ONE / gamma_const)
       enth = pres / dens * gamma_const / gm1
       eint = enth / gamma_const



    case (eos_input_ph)

       ! pressure, enthalpy and xmass are inputs

       ! Solve for the density and energy:

       dens = (pres / K_const)**(ONE / gamma_const)       
       eint = (pres / dens) * ONE / gm1



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

    state % T   = temp
    state % rho = dens
    state % h   = enth
    state % s   = entr
    state % e   = eint
    state % p   = pres

    ! Compute the thermodynamic derivatives and specific heats 
    state % dPdT = ZERO
    state % dPdr = gamma_const * pres / dens
    state % dedT = ZERO
    state % dedr = pres / (dens * dens)
    state % dsdT = ZERO
    state % dsdr = ZERO
    state % dhdT = ZERO
    state % dhdr = state % dedr + gm1 * pres / dens**2

    state % dpde = ZERO
    state % dpdr_e = gamma_const * pres / dens
    
    state % cv = state % dedT
    state % cp = gamma_const * state % cv

    state % gam1 = gamma_const

    ! Compute dPdX, dedX, dhdX.

    state % dpdA = - state % p / state % abar
    state % dpdZ =   state % p / (ONE + state % zbar)

    state % dedA = - state % e / state % abar
    state % dedZ =   state % e / (ONE + state % zbar)

    ! sound speed
    state % cs = sqrt(gamma_const*pres/dens)

  end subroutine actual_eos

end module actual_eos_module
