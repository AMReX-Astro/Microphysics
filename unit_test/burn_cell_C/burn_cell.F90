! Burn a single cell and output time-series integration data.

subroutine burn_cell(name, namlen) bind(C, name="burn_cell")

  use amrex_error_module
  use amrex_constants_module
  use amrex_fort_module, only : rt => amrex_real

  use extern_probin_module, only: run_prefix, small_temp, small_dens
  use extern_probin_module
  use burn_type_module
  use actual_burner_module
  use microphysics_module
  use eos_type_module, only : eos_get_small_temp, eos_get_small_dens, eos_t, &
                              eos_input_rt
  use eos_module
  use network

  implicit none

  integer, intent(in) :: namlen
  integer, intent(in) :: name(namlen)

  type (burn_t)       :: burn_state_in, burn_state_out
  type (eos_t)        :: eos_state_in, eos_state_out

  real (rt)    :: tmax, energy, time, dt
  integer             :: numsteps, i, istate

  character (len=256) :: params_file
  integer             :: params_file_unit

  character (len=256) :: out_directory_name
  character (len=256) :: out_name
  character (len=6)   :: out_num

  ! Starting conditions for integration
  real (rt)    :: density, temperature, massfractions(nspec)

  ! Useful for evaluating final values
  real (rt)    :: eos_energy_generated, eos_energy_rate

  namelist /cellparams/ tmax, numsteps, density, temperature, massfractions

  ! runtime
  call runtime_init(name, namlen)

  ! microphysics
  call microphysics_init(small_temp=small_temp, small_dens=small_dens)
  call eos_get_small_temp(small_temp)
  call eos_get_small_dens(small_dens)

  print *, "small_temp = ", small_temp
  print *, "small_dens = ", small_dens

  ! Set mass fractions to sanitize inputs for them
  massfractions = -1.0d0

  ! Get initial conditions for the burn
  call get_command_argument(1, value = params_file)

  open(newunit=params_file_unit, file=params_file, status="old", action="read")
  read(unit=params_file_unit, nml=cellparams)
  close(unit=params_file_unit)

  ! Make sure user set all the mass fractions to values in the interval [0, 1]
  do i = 1, nspec
     if (massfractions(i) .lt. 0 .or. massfractions(i) .gt. 1) then
        call amrex_error('mass fraction for ' // short_spec_names(i) // ' not initialized in the interval [0,1]!')
     end if
  end do

  ! Echo initial conditions at burn and fill burn state input
  write(*,*) 'Maximum Time (s): ', tmax
  write(*,*) 'Number of time subdivisions: ', numsteps
  write(*,*) 'State Density (g/cm^3): ', density
  write(*,*) 'State Temperature (K): ', temperature
  do i = 1, nspec
     write(*,*) 'Mass Fraction (', short_spec_names(i), '): ', massfractions(i)
  end do

  burn_state_in%T   = temperature
  burn_state_in%rho = density
  burn_state_in%xn(:) = massfractions(:)

  ! normalize -- just in case
  !call normalize_abundances_burn(burn_state_in)

  ! Initialize initial energy to zero
  burn_state_in % e = ZERO
  burn_state_out % e = ZERO
  energy = ZERO

  ! output initial burn type data
  time = ZERO

  ! call the EOS to set initial e
  call burn_to_eos(burn_state_in, eos_state_in)
  call eos(eos_input_rt, eos_state_in)
  call eos_to_burn(eos_state_in, burn_state_in)

  call actual_burner(burn_state_in, burn_state_out, dt, time)
  energy = energy + burn_state_out % e

  ! call the EOS to check consistency of integrated e
  call burn_to_eos(burn_state_out, eos_state_out)
  call eos(eos_input_rt, eos_state_out)

  write(*,*) "------------------------------------"
  write(*,*) "Completed burn to: ", burn_state_out % time, " seconds:"
  write(*,*) " - Hnuc = ", burn_state_out % e / dt
  write(*,*) " - integrated e = ", eos_state_in % e + energy
  write(*,*) " - EOS e(rho, T) = ", eos_state_out % e
  write(*,*) " - integrated/EOS percent diff. = ", 100.0d0 * (eos_state_in % e + energy - eos_state_out % e)/eos_state_out % e

  ! output burn type data
  call write_burn_t(burn_state_out)

  write(*,*) "------------------------------------"
  write(*,*) "EOS e(rho, T) initial = ", eos_state_in % e
  write(*,*) "EOS e(rho, T) final = ", eos_state_out % e
  eos_energy_generated = eos_state_out % e - eos_state_in % e
  write(*,*) "EOS e(rho, T) generated = ", eos_energy_generated
  eos_energy_rate = (eos_state_out % e - eos_state_in % e)/tmax
  write(*,*) "EOS e(rho, T) generation rate = ", eos_energy_rate
  write(*,*) "Integrator total generated energy: ", energy
  write(*,*) "Integrator average energy generation rate: ", energy/tmax
  write(*,*) "(integrator - EOS)/EOS percent diff for generated energy: ", 100.0d0 * (energy - eos_energy_generated)/eos_energy_generated
  write(*,*) "(integrator - EOS)/EOS percent diff for energy gen. rate: ", 100.0d0 * (energy/tmax - eos_energy_rate)/eos_energy_rate

  call microphysics_finalize()

end subroutine burn_cell

subroutine write_burn_t(burnt)
  use network
  use burn_type_module

  implicit none

  ! Writes contents of burn_t type burnt to file named fname
  type(burn_t), intent(in) :: burnt
  character(len=20), parameter :: DPFMT = '(E30.16E5)'
  character(len=20) :: VDPFMT = ''

  integer :: i, j

  write(*, fmt=*) '! Burn Type Data'

  write(*, fmt=*) 'nspec:'
  write(*, fmt=*) nspec

  write(*, fmt=*) 'neqs:'
  write(*, fmt=*) neqs

  write(*, fmt=*) 'short_spec_names:'
  do i = 1, nspec
     write(*, fmt=*) short_spec_names(i)
  end do

  write(*, fmt=*) 'rho:'
  write(*, fmt=DPFMT) burnt % rho

  write(*, fmt=*) 'T:'
  write(*, fmt=DPFMT) burnt % T

  write(*, fmt=*) 'e:'
  write(*, fmt=DPFMT) burnt % e

  write(VDPFMT, '("(", I0, "E30.16E5", ")")') nspec
  write(*, fmt=*) 'xn:'
  write(*, fmt=VDPFMT) (burnt % xn(i), i = 1, nspec)

  write(*, fmt=*) 'cv:'
  write(*, fmt=DPFMT) burnt % cv

  write(*, fmt=*) 'cp:'
  write(*, fmt=DPFMT) burnt % cp

  write(*, fmt=*) 'y_e:'
  write(*, fmt=DPFMT) burnt % y_e

  write(*, fmt=*) 'eta:'
  write(*, fmt=DPFMT) burnt % eta

  write(*, fmt=*) 'cs:'
  write(*, fmt=DPFMT) burnt % cs

  write(*, fmt=*) 'dx:'
  write(*, fmt=DPFMT) burnt % dx

  write(*, fmt=*) 'abar:'
  write(*, fmt=DPFMT) burnt % abar

  write(*, fmt=*) 'zbar:'
  write(*, fmt=DPFMT) burnt % zbar

  write(*, fmt=*) 'T_old:'
  write(*, fmt=DPFMT) burnt % T_old

  write(*, fmt=*) 'dcvdT:'
  write(*, fmt=DPFMT) burnt % dcvdT

  write(*, fmt=*) 'dcpdT:'
  write(*, fmt=DPFMT) burnt % dcpdT

  write(VDPFMT, '("(", I0, "E30.16E5", ")")') neqs
  write(*, fmt=*) 'ydot:'
  write(*, fmt=VDPFMT) (burnt % ydot(i), i = 1, neqs)

  write(*, fmt=*) 'jac:'
  do i = 1, neqs
     write(*, fmt=VDPFMT) (burnt % jac(i,j), j = 1, neqs)
  end do

  write(*, fmt=*) 'self_heat:'
  write(*, fmt=*) burnt % self_heat

  write(*, fmt=*) 'i:'
  write(*, fmt=*) burnt % i

  write(*, fmt=*) 'j:'
  write(*, fmt=*) burnt % j

  write(*, fmt=*) 'k:'
  write(*, fmt=*) burnt % k

  write(*, fmt=*) 'n_rhs:'
  write(*, fmt=*) burnt % n_rhs

  write(*, fmt=*) 'n_jac:'
  write(*, fmt=*) burnt % n_jac

  write(*, fmt=*) 'time:'
  write(*, fmt=DPFMT) burnt % time

end subroutine write_burn_t
