! This is a multi-gamma EOS.  Each species can have its own gamma, but
! otherwise, they all act as an ideal gas.

! Note: at the moment, it is not clear what the proper expression for
! a multi-gamma entropy should be, so do not rely on the entropy.

module actual_eos_module

  use bl_types
  use bl_space
  use bl_error_module
  use bl_constants_module
  use network, only: nspec, aion, zion
  use eos_type_module
  use eos_data_module

  implicit none

  character (len=64) :: eos_name = "multigamma"
  
  double precision :: gammas(nspec)

contains

  subroutine actual_eos_init

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
    print *, "In actual_eos_init, species, gamma: "
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
 
  end subroutine actual_eos_init



  subroutine actual_eos(input, state)

    use fundamental_constants_module, only: k_B, n_A, hbar

    implicit none

    integer,             intent(in   ) :: input
    type (eos_t_vector), intent(inout) :: state

    ! Local variables
    double precision :: sumY_gm1, sumYg_gm1
    double precision :: dens, temp

    ! Get the mass of a nucleon from Avogadro's number.
    double precision, parameter :: m_nucleon = ONE / n_A

    integer :: j

    do j = 1, state % N

       ! Special gamma factors

       sumY_gm1  = sum(state % xn(j,:)/(aion(:)*(gammas(:)-ONE)))
       sumYg_gm1 = sum(state % xn(j,:)*gammas(:)/(aion(:)*(gammas(:)-ONE)))

       !-------------------------------------------------------------------------
       ! For all EOS input modes EXCEPT eos_input_rt, first compute dens
       ! and temp as needed from the inputs.
       !-------------------------------------------------------------------------

       select case (input)

       case (eos_input_rt)

          ! dens, temp and xmass are inputs

          ! We don't need to do anything here
          temp = state % T(j)
          dens = state % rho(j)


       case (eos_input_rh)

          ! dens, enthalpy, and xmass are inputs

          ! Solve for the temperature:
          ! h = e + p/rho = (p/rho)*[1 + 1/(gamma-1)] = (p/rho)*gamma/(gamma-1)
          dens = state % rho(j)
          temp = (state % h(j) * m_nucleon / k_B)/sumYg_gm1


       case (eos_input_tp)

          ! temp, pres, and xmass are inputs

          ! Solve for the density:
          ! p = rho k T / (abar m_nucleon)
          dens = state % p(j) * state % abar(j) * m_nucleon / (k_B * state % T(j))
          temp = state % T(j)


       case (eos_input_rp)

          ! dens, pres, and xmass are inputs

          ! Solve for the temperature:
          ! p = rho k T / (mu m_nucleon)
          dens = state % rho(j)
          temp = state % p(j) * state % abar(j) * m_nucleon / (k_B * state % rho(j))


       case (eos_input_re)

          ! dens, energy, and xmass are inputs

          ! Solve for the temperature
          ! e = k T / [(mu m_nucleon)*(gamma-1)]
          dens = state % rho(j)
          temp = state % e(j) * m_nucleon / (k_B * sumY_gm1)


       case (eos_input_ps)

          ! pressure entropy, and xmass are inputs
          call bl_error("eos_input_ps is not supported")


       case (eos_input_ph)

          ! pressure, enthalpy and xmass are inputs

          ! Solve for temperature and density
          dens = state % p(j) * state % abar(j) / state % h(j) * sumYg_gm1
          temp = state % p(j) * state % abar(j) * m_nucleon / (k_B * dens)


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

       state % T(j)   = temp
       state % rho(j) = dens

       ! Compute the pressure simply from the ideal gas law, and the
       ! specific internal energy using the gamma-law EOS relation.
       state % p(j) = dens*k_B*temp/(state % abar(j)*m_nucleon)
       state % e(j) = k_B*temp*sumY_gm1/m_nucleon

       ! Enthalpy is h = e + p/rho
       state % h(j) = state % e(j) + state % p(j) / dens

       ! Entropy (per gram) -- not implemented.
       state % s(j) = ZERO

       ! Compute the thermodynamic derivatives and specific heats 
       state % dpdT = state % p(j) / temp
       state % dpdr = state % p(j) / dens
       state % dedT(j) = state % e(j) / temp
       state % dedr(j) = ZERO
       state % dsdT(j) = ZERO
       state % dsdr(j) = ZERO
       state % dhdT(j) = state % dedT(j) + state % dpdT(j) / dens
       state % dhdr(j) = ZERO

       state % cv(j) = state % dedT(j)
       state % cp(j) = state % h(j) / state % T(j)

       state % gam1(j) = state % cp(j) / state % cv(j)

       state % dpdr_e(j) = state % dpdr(j) - state % dpdT(j) * state % dedr(j) / state % dedT(j)
       state % dpde(j)   = state % dpdT(j) / state % dedT(j)


       ! These need to be worked out.
       state % dpdA(j) = ZERO
       state % dedA(j) = ZERO

       state % dpdZ(j) = ZERO
       state % dedZ(j) = ZERO

       ! These expressions will not work at present because we overwrite them
       ! at the end of the EOS call. We need to move the multi-fluid property
       ! up to the main EOS modules.

       state % dpdX(j,:) = state % rho(j) * k_B * state % T(j) / (m_nucleon * aion(:))
       state % dedX(j,:) = k_B * state % T(j) / (m_nucleon * aion(:) * (gammas(:) - ONE))

       state % dhdX(j,:) = state % dedX(j,:) + (state % p(j) / state % rho(j)**2 - state % dedr(j)) &
                                            *  state % dpdX(j,:) / state % dpdr(j)

       ! Sound speed
       state % cs(j) = sqrt(state % gam1(j) * state % p(j) / dens)

    enddo

  end subroutine actual_eos

end module actual_eos_module
