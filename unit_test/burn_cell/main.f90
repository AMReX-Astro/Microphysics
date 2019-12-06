! Burn a single cell and output time-series integration data.

program burn_cell

  use amrex_error_module
  use amrex_constants_module
  use microphysics_type_module

  use fabio_module, only: fabio_mkdir
  use probin_module, only: run_prefix, small_temp, small_dens
  use runtime_init_module
  use extern_probin_module
  use burn_type_module
  use actual_burner_module
  use microphysics_module
  use eos_type_module, only : eos_get_small_temp, eos_get_small_dens, eos_t, &
                              eos_input_rt
  use eos_module
  use network
  use build_info_module

  implicit none

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
  call runtime_init(.true.)

  ! microphysics
  call microphysics_init(small_temp=small_temp, small_dens=small_dens)
  call eos_get_small_temp(small_temp)
  print *, "small_temp = ", small_temp
  call eos_get_small_dens(small_dens)
  print *, "small_dens = ", small_dens

  ! Set mass fractions to sanitize inputs for them
  massfractions = -1.0e0_rt

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

  ! Create directory to hold output
  out_directory_name = trim(run_prefix) // 'output'
  call fabio_mkdir(trim(out_directory_name), istate)
  if (istate /= 0) then
     call amrex_error('output directory "' // trim(out_directory_name) // '" could not be created.')
  end if

  ! output initial burn type data
  time = ZERO

  ! call the EOS to set initial e
  call burn_to_eos(burn_state_in, eos_state_in)
  call eos(eos_input_rt, eos_state_in)
  call eos_to_burn(eos_state_in, burn_state_in)

  write(out_num,'(I6.6)') 0
  out_name = trim(out_directory_name) // '/' // trim(run_prefix) // out_num
  burn_state_in % time = time
  call write_burn_t(out_name, burn_state_in)
  
  dt = tmax/numsteps

  do i = 1, numsteps
     ! Do burn
     call actual_burner(burn_state_in, burn_state_out, dt, time)
     energy = energy + burn_state_out % e
     burn_state_in = burn_state_out
     burn_state_in % e = ZERO

     ! call the EOS to check consistency of integrated e
     call burn_to_eos(burn_state_out, eos_state_out)
     call eos(eos_input_rt, eos_state_out)

     write(*,*) "------------------------------------"
     write(*,*) "Completed burn to: ", burn_state_out % time, " seconds:"
     write(*,*) " - Hnuc = ", burn_state_out % e / dt
     write(*,*) " - integrated e = ", eos_state_in % e + energy
     write(*,*) " - EOS e(rho, T) = ", eos_state_out % e
     write(*,*) " - integrated/EOS percent diff. = ", 100.0e0_rt * (eos_state_in % e + energy - eos_state_out % e)/eos_state_out % e
     
     ! output burn type data
     write(out_num,'(I6.6)') i
     out_name = trim(out_directory_name) // '/' // trim(run_prefix) // out_num
     call write_burn_t(out_name, burn_state_out)

     time = burn_state_out % time
  end do

  write(*,*) "------------------------------------"
  write(*,*) "EOS e(rho, T) initial = ", eos_state_in % e
  write(*,*) "EOS e(rho, T) final = ", eos_state_out % e
  eos_energy_generated = eos_state_out % e - eos_state_in % e
  write(*,*) "EOS e(rho, T) generated = ", eos_energy_generated
  eos_energy_rate = (eos_state_out % e - eos_state_in % e)/tmax
  write(*,*) "EOS e(rho, T) generation rate = ", eos_energy_rate
  write(*,*) "Integrator total generated energy: ", energy
  write(*,*) "Integrator average energy generation rate: ", energy/tmax
  write(*,*) "(integrator - EOS)/EOS percent diff for generated energy: ", 100.0e0_rt * (energy - eos_energy_generated)/eos_energy_generated
  write(*,*) "(integrator - EOS)/EOS percent diff for energy gen. rate: ", 100.0e0_rt * (energy/tmax - eos_energy_rate)/eos_energy_rate

  call microphysics_finalize()

end program burn_cell

subroutine write_burn_t(fname, burnt)
  use network
  use burn_type_module

  implicit none
  
  ! Writes contents of burn_t type burnt to file named fname
  character(len=256), intent(in) :: fname
  type(burn_t), intent(in) :: burnt
  integer, parameter :: burn_t_unit = 10
  character(len=20), parameter :: DPFMT = '(E30.16E5)'
  character(len=20) :: VDPFMT = ''
  
  integer :: i, j

  open(unit=burn_t_unit, file=fname, action='WRITE')
  write(unit=burn_t_unit, fmt=*) '! Burn Type Data'

  write(unit=burn_t_unit, fmt=*) 'nspec:'
  write(unit=burn_t_unit, fmt=*) nspec
  
  write(unit=burn_t_unit, fmt=*) 'neqs:'
  write(unit=burn_t_unit, fmt=*) neqs

  write(unit=burn_t_unit, fmt=*) 'short_spec_names:'
  do i = 1, nspec
     write(unit=burn_t_unit, fmt=*) short_spec_names(i)
  end do
  
  write(unit=burn_t_unit, fmt=*) 'rho:'
  write(unit=burn_t_unit, fmt=DPFMT) burnt % rho

  write(unit=burn_t_unit, fmt=*) 'T:'
  write(unit=burn_t_unit, fmt=DPFMT) burnt % T

  write(unit=burn_t_unit, fmt=*) 'e:'
  write(unit=burn_t_unit, fmt=DPFMT) burnt % e

  write(VDPFMT, '("(", I0, "E30.16E5", ")")') nspec
  write(unit=burn_t_unit, fmt=*) 'xn:'
  write(unit=burn_t_unit, fmt=VDPFMT) (burnt % xn(i), i = 1, nspec)

  write(unit=burn_t_unit, fmt=*) 'cv:'
  write(unit=burn_t_unit, fmt=DPFMT) burnt % cv

  write(unit=burn_t_unit, fmt=*) 'cp:'
  write(unit=burn_t_unit, fmt=DPFMT) burnt % cp

  write(unit=burn_t_unit, fmt=*) 'y_e:'
  write(unit=burn_t_unit, fmt=DPFMT) burnt % y_e

  write(unit=burn_t_unit, fmt=*) 'eta:'
  write(unit=burn_t_unit, fmt=DPFMT) burnt % eta

  write(unit=burn_t_unit, fmt=*) 'cs:'
  write(unit=burn_t_unit, fmt=DPFMT) burnt % cs

  write(unit=burn_t_unit, fmt=*) 'dx:'
  write(unit=burn_t_unit, fmt=DPFMT) burnt % dx

  write(unit=burn_t_unit, fmt=*) 'abar:'
  write(unit=burn_t_unit, fmt=DPFMT) burnt % abar

  write(unit=burn_t_unit, fmt=*) 'zbar:'
  write(unit=burn_t_unit, fmt=DPFMT) burnt % zbar

  write(unit=burn_t_unit, fmt=*) 'T_old:'
  write(unit=burn_t_unit, fmt=DPFMT) burnt % T_old

  write(unit=burn_t_unit, fmt=*) 'dcvdT:'
  write(unit=burn_t_unit, fmt=DPFMT) burnt % dcvdT

  write(unit=burn_t_unit, fmt=*) 'dcpdT:'
  write(unit=burn_t_unit, fmt=DPFMT) burnt % dcpdT

  write(VDPFMT, '("(", I0, "E30.16E5", ")")') neqs
  write(unit=burn_t_unit, fmt=*) 'ydot:'
  write(unit=burn_t_unit, fmt=VDPFMT) (burnt % ydot(i), i = 1, neqs)

  write(unit=burn_t_unit, fmt=*) 'jac:'
  do i = 1, neqs
     write(unit=burn_t_unit, fmt=VDPFMT) (burnt % jac(i,j), j = 1, neqs)
  end do

  write(unit=burn_t_unit, fmt=*) 'self_heat:'
  write(unit=burn_t_unit, fmt=*) burnt % self_heat

  write(unit=burn_t_unit, fmt=*) 'i:'
  write(unit=burn_t_unit, fmt=*) burnt % i
  
  write(unit=burn_t_unit, fmt=*) 'j:'
  write(unit=burn_t_unit, fmt=*) burnt % j

  write(unit=burn_t_unit, fmt=*) 'k:'
  write(unit=burn_t_unit, fmt=*) burnt % k

  write(unit=burn_t_unit, fmt=*) 'n_rhs:'
  write(unit=burn_t_unit, fmt=*) burnt % n_rhs

  write(unit=burn_t_unit, fmt=*) 'n_jac:'
  write(unit=burn_t_unit, fmt=*) burnt % n_jac

  write(unit=burn_t_unit, fmt=*) 'time:'
  write(unit=burn_t_unit, fmt=DPFMT) burnt % time
  
  close(unit=burn_t_unit)
end subroutine write_burn_t

