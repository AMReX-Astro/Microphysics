! Burn a single cell and output time-series integration data.

program burn_cell

  use bl_constants_module
  use bl_types
  use probin_module, only: run_prefix, small_temp, small_dens
  use burn_type_module
  use actual_burner_module
  use microphysics_module
  use eos_type_module, only : eos_get_small_temp, eos_get_small_dens
  use eos_module
  use network
  use build_info_module

  implicit none

  type (burn_t)       :: burn_state_in, burn_state_out
  type (burn_t), allocatable :: burn_state_history(:)

  real (kind=dp_t)    :: tmax, energy, dt
  integer             :: ntimes, i

  character (len=256) :: out_name

  ! microphysics
  call microphysics_init(small_temp=small_temp, small_dens=small_dens)
  call eos_get_small_temp(small_temp)
  print *, "small_temp = ", small_temp
  call eos_get_small_dens(small_dens)
  print *, "small_dens = ", small_dens

  ! Fill burn state input
  write(*,*) 'Maximum Time (s): '
  read(*,*) tmax
  write(*,*) 'Number of time subdivisions: '
  read(*,*) ntimes
  write(*,*) 'State Density (g/cm^3): '
  read(*,*) burn_state_in%rho
  write(*,*) 'State Temperature (K): '
  read(*,*) burn_state_in%T
  do i = 1, nspec
     write(*,*) 'Mass Fraction (', spec_names(i), '): '
     read(*,*) burn_state_in%xn(i)
  end do

  ! normalize -- just in case
  !call normalize_abundances_burn(burn_state_in)

  ! Allocate burn history
  allocate( burn_state_history(ntimes+1) )
  
  ! Initialize initial energy to zero
  burn_state_in % e = ZERO
  burn_state_out % e = ZERO
  energy = ZERO

  dt = tmax/ntimes
  burn_state_history(1) = burn_state_in
  
  do i = 1, ntimes
     call actual_burner(burn_state_in, burn_state_out, dt, ZERO)
     energy = energy + burn_state_out % e
     write(*,'(E25.16)') energy
     burn_state_history(i+1) = burn_state_out
     burn_state_in = burn_state_out
     burn_state_in % e = ZERO
  end do

  ! output burn history
  out_name = trim(run_prefix) // "test_react." // trim(integrator_dir)

  ! TODO

  ! Free burn history and finalize
  deallocate( burn_state_history )
  call microphysics_finalize()

end program burn_cell
