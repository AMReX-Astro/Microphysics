! This is a multi-gamma EOS.  Each species can have its own gamma, but
! otherwise, they all act as an ideal gas.

! Note: at the moment, it is not clear what the proper expression for
! a multi-gamma entropy should be, so do not rely on the entropy.

module specific_eos_module

  use bl_types
  use bl_space
  use bl_error_module
  use bl_constants_module
  use network, only: nspec, aion, zion
  use eos_type_module
  use eos_data_module

  implicit none

  double precision :: gammas(nspec)

contains

  subroutine specific_eos_init

    use extern_probin_module, only: eos_gamma_default, &
                                    species_a_name, species_a_gamma, &
                                    species_b_name, species_b_gamma, &
                                    species_c_name, species_c_gamma

    implicit none
 
    integer :: idx

    ! Set the gammas for the species -- we have some runtime parameters
    ! that can override the default gammas for a few named species.
    gammas(:) = eos_gamma_default

    print *, ""
    print *, "In specific_eos_init, species, gamma: "
    print *, "species a: ", trim(species_a_name), species_a_gamma
    print *, "species b: ", trim(species_b_name), species_b_gamma
    print *, "species c: ", trim(species_c_name), species_c_gamma
    print *, ""

    idx = network_species_index(species_a_name)
    if (idx > 0) gammas(idx) = species_a_gamma

    idx = network_species_index(species_b_name)
    if (idx > 0) gammas(idx) = species_b_gamma

    idx = network_species_index(species_c_name)
    if (idx > 0) gammas(idx) = species_c_gamma

    initialized = .true.
 
  end subroutine specific_eos_init



  subroutine specific_eos(input, state)

    use fundamental_constants_module, only: k_B, n_A, hbar

    implicit none

    integer,           intent(in   ) :: input
    type (eos_t),      intent(inout) :: state(:)

    ! Local variables
    double precision :: sumY_gm1, sumYg_gm1
    double precision :: dens, temp

    ! Get the mass of a nucleon from Avogadro's number.
    double precision, parameter :: m_nucleon = ONE / n_A

    integer :: j, N

    N = size(state)

    if (.not. initialized) call bl_error('EOS: not initialized')

    do j = 1, N

       ! Special gamma factors

       sumY_gm1  = sum(state(j) % xn(:)/(aion(:)*(gammas(:)-ONE)))
       sumYg_gm1 = sum(state(j) % xn(:)*gammas(:)/(aion(:)*(gammas(:)-ONE)))

       !-------------------------------------------------------------------------
       ! For all EOS input modes EXCEPT eos_input_rt, first compute dens
       ! and temp as needed from the inputs.
       !-------------------------------------------------------------------------

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
          temp = (state(j) % h * m_nucleon / k_B)/sumYg_gm1


       case (eos_input_tp)

          ! temp, pres, and xmass are inputs

          ! Solve for the density:
          ! p = rho k T / (abar m_nucleon)
          dens = state(j) % p * state(j) % abar * m_nucleon / (k_B * state(j) % T)
          temp = state(j) % T


       case (eos_input_rp)

          ! dens, pres, and xmass are inputs

          ! Solve for the temperature:
          ! p = rho k T / (mu m_nucleon)
          dens = state(j) % rho
          temp = state(j) % p * state(j) % abar * m_nucleon / (k_B * state(j) % rho)


       case (eos_input_re)

          ! dens, energy, and xmass are inputs

          ! Solve for the temperature
          ! e = k T / [(mu m_nucleon)*(gamma-1)]
          dens = state(j) % rho
          temp = state(j) % e * m_nucleon / (k_B * sumY_gm1)


       case (eos_input_ps)

          ! pressure entropy, and xmass are inputs
          call bl_error("eos_input_ps is not supported")


       case (eos_input_ph)

          ! pressure, enthalpy and xmass are inputs

          ! Solve for temperature and density
          dens = state(j) % p * state(j) % abar / state(j) % h * sumYg_gm1
          temp = state(j) % p * state(j) % abar * m_nucleon / (k_B * dens)


       case (eos_input_th)

          ! temperature, enthalpy and xmass are inputs
          call bl_error("eos_input_th is not supported")


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
       state(j) % p = dens*k_B*temp/(state(j) % abar*m_nucleon)
       state(j) % e = k_B*temp*sumY_gm1/m_nucleon

       ! Enthalpy is h = e + p/rho
       state(j) % h = state(j) % e + state(j) % p / dens

       ! Entropy (per gram) -- not implemented.
       state(j) % s = ZERO

       ! Compute the thermodynamic derivatives and specific heats 
       state(j) % dpdT = state(j) % p / temp
       state(j) % dpdr = state(j) % p / dens
       state(j) % dedT = state(j) % e / temp
       state(j) % dedr = ZERO
       state(j) % dsdT = ZERO
       state(j) % dsdr = ZERO
       state(j) % dhdT = state(j) % dedT + state(j) % dpdT / dens
       state(j) % dhdr = 0.0d0

       state(j) % cv = state(j) % dedT
       state(j) % cp = state(j) % h / state(j) % T

       state(j) % gam1 = state(j) % cp / state(j) % cv

       state(j) % dpdr_e = state(j) % dpdr - state(j) % dpdT * state(j) % dedr / state(j) % dedT
       state(j) % dpde   = state(j) % dpdT / state(j) % dedT


       ! These need to be worked out.
       state(j) % dpdA = ZERO
       state(j) % dedA = ZERO

       state(j) % dpdZ = ZERO
       state(j) % dedZ = ZERO

       ! These expressions will not work at present because we overwrite them
       ! at the end of the EOS call. We need to move the multi-fluid property
       ! up to the main EOS modules.

       state(j) % dpdX(:) = state(j) % rho * k_B * state(j) % T / (m_nucleon * aion(:))
       state(j) % dedX(:) = k_B * state(j) % T / (m_nucleon * aion(:) * (gammas(:) - ONE))

       state(j) % dhdX(:) = state(j) % dedX(:) + (state(j) % p / state(j) % rho**2 - state(j) % dedr) &
                                            *  state(j) % dpdX(:) / state(j) % dpdr

       ! Sound speed
       state(j) % cs = sqrt(state(j) % gam1 * state(j) % p / dens)

    enddo

  end subroutine specific_eos

end module specific_eos_module
