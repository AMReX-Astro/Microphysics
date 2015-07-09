! This is a constant gamma equation of state, using an ideal gas.
!
! The gas may either be completely ionized or completely neutral.
!
! The ratio of specific heats (gamma) is allowed to vary.  NOTE: the
! expression for entropy is only valid for an ideal MONATOMIC gas
! (gamma = 5/3).  

module specific_eos_module

  use bl_types
  use bl_space
  use bl_error_module
  use bl_constants_module
  use network, only: nspec, aion, zion
  use eos_type_module
  use eos_data_module

  implicit none

  double precision  :: gamma_const

contains

  subroutine specific_eos_init

    use extern_probin_module, only: eos_gamma

    implicit none
 
    ! constant ratio of specific heats
    if (eos_gamma .gt. 0.d0) then
       gamma_const = eos_gamma
    else
       gamma_const = FIVE3RD
    end if
 
    initialized = .true.
 
  end subroutine specific_eos_init



  subroutine specific_eos(eosfail, state, input)

    use fundamental_constants_module, only: k_B, n_A, hbar

    implicit none

    logical,           intent(inout) :: eosfail
    type (eos_t),      intent(inout) :: state(:)
    integer,           intent(in   ) :: input

    ! Local variables
    double precision :: dens, temp

    ! Get the mass of a nucleon from Avogadro's number.
    double precision, parameter :: m_nucleon = ONE / n_A

    integer :: j, N

    N = size(state)

    if (.not. initialized) call bl_error('EOS: not initialized')

    !-------------------------------------------------------------------------
    ! For all EOS input modes EXCEPT eos_input_rt, first compute dens
    ! and temp as needed from the inputs.
    !-------------------------------------------------------------------------

    do j = 1, N

       select case (input)

       case (eos_input_rt)

          ! dens, temp and xmass are inputs

          ! We don't need to do anything here
          temp = state(j) % T
          dens = state(j) % rho


       case (eos_input_rh)

          ! dens, enthalpy, and xmass are inputs

          ! Solve for the temperature:
          ! h = e + p/rho = (p/rho)*[1 + 1/(gamma-1)] = (p/rho)*gamma/(gamma-1)
          dens = state(j) % rho
          temp = (state(j) % h * state(j) % mu * m_nucleon / k_B)*(gamma_const - ONE)/gamma_const


       case (eos_input_tp)

          ! temp, pres, and xmass are inputs

          ! Solve for the density:
          ! p = rho k T / (mu m_nucleon)
          dens = state(j) % p * state(j) % mu * m_nucleon / (k_B * state(j) % T)
          temp = state(j) % T


       case (eos_input_rp)

          ! dens, pres, and xmass are inputs

          ! Solve for the temperature:
          ! p = rho k T / (mu m_nucleon)
          dens = state(j) % rho
          temp = state(j) % p * state(j) % mu * m_nucleon / (k_B * state(j) % rho)


       case (eos_input_re)

          ! dens, energy, and xmass are inputs

          ! Solve for the temperature
          ! e = k T / [(mu m_nucleon)*(gamma-1)]
          dens = state(j) % rho
          temp = state(j) % e * state(j) % mu * m_nucleon * (gamma_const-ONE) / k_B


       case (eos_input_ps)

          ! pressure entropy, and xmass are inputs

          ! Solve for the temperature
          ! Invert Sackur-Tetrode eqn (below) using 
          ! rho = p mu m_nucleon / (k T)
          temp = state(j) % p**(TWO/FIVE) * &
               ( TWO*M_PI*hbar*hbar/(state(j) % mu*m_nucleon) )**(THREE/FIVE) * &
               dexp(TWO*state(j) % mu*m_nucleon*state(j) % s/(FIVE*k_B) - ONE) / k_B

          ! Solve for the density
          ! rho = p mu m_nucleon / (k T)
          dens = state(j) % p * state(j) % mu * m_nucleon / (k_B * temp)



       case (eos_input_ph)

          ! pressure, enthalpy and xmass are inputs

          ! Solve for temperature and density
          dens = state(j) % p / state(j) % h * gamma_const / (gamma_const - ONE)
          temp = state(j) % p * state(j) % mu * m_nucleon / (k_B * dens)



       case (eos_input_th)

          ! temperature, enthalpy and xmass are inputs

          ! This system is underconstrained.

          call bl_error('EOS: eos_input_th is not a valid input for the gamma law EOS.')



       case default

          call bl_error('EOS: invalid input.')

       end select

       !-------------------------------------------------------------------------
       ! Now we have the density and temperature (and mass fractions /
       ! mu), regardless of the inputs.
       !-------------------------------------------------------------------------

       state(j) % T   = temp
       state(j) % rho = dens

       ! Compute the pressure simply from the ideal gas law, and the
       ! specific internal energy using the gamma-law EOS relation.
       state(j) % p = dens*k_B*temp/(state(j) % mu*m_nucleon)
       state(j) % e = state(j) % p/(gamma_const - ONE)/dens

       ! enthalpy is h = e + p/rho
       state(j) % h = state(j) % e + state(j) % p / dens

       ! entropy (per gram) of an ideal monoatomic gas (the Sackur-Tetrode equation)
       ! NOTE: this expression is only valid for gamma = 5/3.
       state(j) % s = (k_B/(state(j) % mu*m_nucleon))*(2.5_dp_t + &
            log( ( (state(j) % mu*m_nucleon)**2.5/dens )*(k_B*temp)**1.5_dp_t / (TWO*M_PI*hbar*hbar)**1.5_dp_t ) )

       ! Compute the thermodynamic derivatives and specific heats 
       state(j) % dpdT = state(j) % p / temp
       state(j) % dpdr = state(j) % p / dens
       state(j) % dedT = state(j) % e / temp
       state(j) % dedr = ZERO
       state(j) % dsdT = THREE / TWO * k_B / (state(j) % mu * m_nucleon) / temp
       state(j) % dsdr = - k_B / (state(j) % mu * m_nucleon) / dens
       state(j) % dhdT = state(j) % dedT + state(j) % dpdT / dens
       state(j) % dhdr = ZERO

       state(j) % cv = state(j) % dedT
       state(j) % cp = gamma_const * state(j) % cv

       state(j) % gam1 = gamma_const

       state(j) % dpdr_e = state(j) % dpdr - state(j) % dpdT * state(j) % dedr / state(j) % dedT
       state(j) % dpde   = state(j) % dpdT / state(j) % dedT

       ! sound speed
       state(j) % cs = sqrt(gamma_const * state(j) % p / dens)

       state(j) % dpdA = - state(j) % p / state(j) % abar
       state(j) % dedA = - state(j) % e / state(j) % abar

       if (assume_neutral) then
         state(j) % dpdZ = ZERO
         state(j) % dedZ = ZERO
       else
         state(j) % dpdZ = state(j) % p / (ONE + state(j) % zbar)
         state(j) % dedZ = state(j) % e / (ONE + state(j) % zbar)
       endif

    enddo

  end subroutine specific_eos

end module specific_eos_module
