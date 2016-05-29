module actual_network

  use actual_network_data

  implicit none

  ! M12_chamulak is the effective number of C12 nuclei destroyed per
  ! reaction
  double precision, parameter :: M12_chamulak = 2.93d0

contains

  subroutine actual_network_init

    implicit none

    spec_names(ic12_)  = "carbon-12"
    spec_names(io16_)  = "oxygen-16"
    spec_names(iash_)  = "ash"

    short_spec_names(ic12_) = "C12"
    short_spec_names(io16_) = "O16"
    short_spec_names(iash_) = "ash"

    ! the ash from C12 burning according to Chamulak et al. is a mixture
    ! of C13, O16, Ne20, and Na23.   Ne20 + alpha results 60% of the time,
    ! while Na23 + p is the other 40%.   Fusing 6 C12 will result in
    ! 1.2 Na23, 1.2 O16 (from the alpha), 1.8 Ne20, and 1.8 C13.
    ! The ash state will have an A and Z corresponding to this mixture.
    aion(ic12_)  = 12.0d0
    aion(io16_)  = 16.0d0
    aion(iash_)  = 18.0d0

    zion(ic12_)  = 6.0d0
    zion(io16_)  = 8.0d0
    zion(iash_)  = 8.8d0

  end subroutine actual_network_init



  subroutine get_ebin(density, ebin)

    use bl_constants_module, only: ZERO
    use fundamental_constants_module

    implicit none

    double precision :: density, ebin(nspec)

    double precision :: rho9, q_eff

    ebin = ZERO

    ! Chamulak et al. provide the q-value resulting from C12 burning,
    ! given as 3 different values (corresponding to 3 different densities).
    ! Here we do a simple quadratic fit to the 3 values provided (see
    ! Chamulak et al., p. 164, column 2).

    ! our convention is that the binding energies are negative.  We convert
    ! from the MeV values that are traditionally written in astrophysics
    ! papers by multiplying by 1.e6 eV/MeV * 1.60217646e-12 erg/eV.  The
    ! MeV values are per nucleus, so we divide by aion to make it per
    ! nucleon and we multiple by Avogardo's # (6.0221415e23) to get the
    ! value in erg/g
    rho9 = density/1.0e9_dp_t

    ! q_eff is effective heat evolved per reaction (given in MeV)
    q_eff = 0.06_dp_t*rho9**2 + 0.02_dp_t*rho9 + 8.83_dp_t

    ! convert from MeV to ergs / gram.  Here M12_chamulak is the
    ! number of C12 nuclei destroyed in a single reaction and 12.0 is
    ! the mass of a single C12 nuclei.  Also note that our convention
    ! is that binding energies are negative.
    ebin(iC12_) = -q_eff*MeV2eV*eV2erg*n_A/(M12_chamulak*12.0_dp_t)

  end subroutine get_ebin



  subroutine actual_network_finalize

    implicit none

    ! Nothing to do here.

  end subroutine actual_network_finalize

end module actual_network
