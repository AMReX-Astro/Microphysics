module weaklib_type_module

  use amrex_fort_module, only: rt => amrex_real
  use eos_type_module
  use fundamental_constants_module, only: k_B, ev2erg, MeV2eV  

  implicit none

  public

  real(rt), parameter :: temp_conv = k_B / ev2erg / MeV2eV
  real(rt), parameter :: energy_shift = ZERO

  type :: weaklib_eos_t
     real(rt) :: density
     real(rt) :: temperature
     real(rt) :: electron_fraction
     
     real(rt) :: pressure
     real(rt) :: entropy_per_baryon
     real(rt) :: specific_internal_energy
     
     real(rt) :: electron_chemical_potential
     real(rt) :: proton_chemical_potential
     real(rt) :: neutron_chemical_potential
     
     real(rt) :: proton_mass_fraction
     real(rt) :: neutron_mass_fraction
     real(rt) :: alpha_mass_fraction
     real(rt) :: heavy_mass_fraction
     
     real(rt) :: gamma_one
     real(rt) :: speed_of_sound
  end type weaklib_eos_t

contains

  subroutine eos_to_weaklib(eos_state, weaklib_state)

    use amrex_constants_module, only: ZERO

    implicit none

    type (eos_t), intent(in) :: eos_state
    type (weaklib_eos_t), intent(out) :: weaklib_state

    call convert_to_table_format(eos_state)

    weaklib_state % density = eos_state % rho
    weaklib_state % temperature = eos_state % T
    weaklib_state % electron_fraction = eos_state % ye
    
    weaklib_state % pressure = eos_state % p
    weaklib_state % entropy_per_baryon = ZERO
    weaklib_state % specific_internal_energy = eos_state % e

    weaklib_state % electron_chemical_potential = ZERO
    weaklib_state % proton_chemical_potential = ZERO
    weaklib_state % neutron_chemical_potential = ZERO

    weaklib_state % proton_mass_fraction = ZERO
    weaklib_state % neutron_mass_fraction = ZERO
    weaklib_state % alpha_mass_fraction = ZERO
    weaklib_state % heavy_mass_fraction = ZERO

    weaklib_state % gamma_one = ZERO

  end subroutine eos_to_weaklib


  subroutine weaklib_to_eos(weaklib_state, eos_state)

    implicit none

    type (weaklib_eos_t), intent(in) :: weaklib_state

    ! eos_state is intent(inout) so we don't reset
    ! quantities that are calculated outside this EOS
    type (eos_t), intent(inout) :: eos_state

    eos_state % rho = weaklib_state % density
    eos_state % T   = weaklib_state % temperature
    eos_state % ye  = weaklib_state % electron_fraction
    
    eos_state % p = weaklib_state % pressure
    eos_state % s = weaklib_state % entropy_per_baryon    
    eos_state % e = weaklib_state % specific_internal_energy

    eos_state % gam1 = weaklib_state % gamma_one

    call convert_from_table_format(eos_state)

    ! construct enthalpy
    eos_state % h = eos_state % e + eos_state % p/eos_state % rho

    ! construct sound speed using the first adiabatic index
    eos_state % cs = sqrt(eos_state % gam1 * eos_state % p / eos_state % rho)

  end subroutine weaklib_to_eos


  ! Convert from the units used in Castro to the units of the table.
  subroutine convert_to_table_format(state)

    use fundamental_constants_module, only: k_B, ev2erg, MeV2eV, n_A
    use amrex_constants_module, only: ZERO, ONE

    implicit none

    type(eos_t), intent(inout) :: state

    ! the weaklib tables use some log10 variables, as well as 
    ! units of MeV for temperature and chemical potential, and k_B / baryon 
    ! for entropy
    if (state % rho > ZERO) then
       state % rho = dlog10(state % rho)
    endif
    if (state % p > ZERO) then
       state % p = dlog10(state % p)
    endif
    if (state % e > ZERO) then
       state % e = dlog10(state % e - energy_shift)
    endif
    if (state % T > ZERO) then
       state % T = dlog10(state % T * temp_conv)
    endif
    ! assuming baryon mass to be ~ 1 amu = 1/N_A
    state % s = state % s * k_B / n_A

  end subroutine convert_to_table_format

  
  ! this converts from the units used in the table to the units of Castro
  subroutine convert_from_table_format(state)

    use fundamental_constants_module
    use amrex_constants_module, only: TEN

    implicit none

    type(eos_t), intent(inout) :: state

    state % rho = TEN**state % rho
    state % p   = TEN**state % p
    state % e   = TEN**state % e + energy_shift
    state % T = (TEN**state % T) / temp_conv
    state % s = state % s * n_A / k_B

  end subroutine convert_from_table_format
  
end module weaklib_type_module
