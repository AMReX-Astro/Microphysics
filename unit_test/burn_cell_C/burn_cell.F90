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

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer, intent(in) :: namlen
  integer, intent(in) :: name(namlen)

  type (burn_t)       :: burn_state_in, burn_state_out
  type (eos_t)        :: eos_state_in, eos_state_out

  real(rt)     :: energy, time, dt
  integer             :: i, istate

  character (len=256) :: params_file
  integer             :: params_file_unit

  character (len=256) :: out_directory_name
  character (len=256) :: out_name
  character (len=6)   :: out_num

  ! Starting conditions for integration
  real(rt)     :: massfractions(nspec)

  ! Useful for evaluating final values
  real(rt)     :: eos_energy_generated, eos_energy_rate

  ! runtime
  call runtime_init(name, namlen)

  ! microphysics
  call microphysics_init(small_temp=small_temp, small_dens=small_dens)
  call eos_get_small_temp(small_temp)
  call eos_get_small_dens(small_dens)

  print *, "small_temp = ", small_temp
  print *, "small_dens = ", small_dens

  ! Set mass fractions to sanitize inputs for them
  massfractions = -1.0e0_rt

  ! Make sure user set all the mass fractions to values in the interval [0, 1]
  do i = 1, nspec
     select case (i)
     case (1)
        massfractions(i) = X1
     case (2)
        massfractions(i) = X2
     case (3)
        massfractions(i) = X3
     case (4)
        massfractions(i) = X4
     case (5)
        massfractions(i) = X5
     case (6)
        massfractions(i) = X6
     case (7)
        massfractions(i) = X7
     case (8)
        massfractions(i) = X8
     case (9)
        massfractions(i) = X9
     case (10)
        massfractions(i) = X10
     case (11)
        massfractions(i) = X11
     case (12)
        massfractions(i) = X12
     case (13)
        massfractions(i) = X13
     case (14)
        massfractions(i) = X14
     case (15)
        massfractions(i) = X15
     case (16)
        massfractions(i) = X16
     case (17)
        massfractions(i) = X17
     case (18)
        massfractions(i) = X18
     case (19)
        massfractions(i) = X19
     case (20)
        massfractions(i) = X20
     case (21)
        massfractions(i) = X21
     end select

     if (massfractions(i) .lt. 0 .or. massfractions(i) .gt. 1) then
        call amrex_error('mass fraction for ' // short_spec_names(i) // ' not initialized in the interval [0,1]!')
     end if
  end do

  ! Echo initial conditions at burn and fill burn state input
  write(*,*) 'Maximum Time (s): ', tmax
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

  dt = tmax
  call actual_burner(burn_state_in, burn_state_out, dt, time)
  energy = energy + burn_state_out % e

  ! call the EOS to check consistency of integrated e
  call burn_to_eos(burn_state_out, eos_state_out)
  call eos(eos_input_rt, eos_state_out)

  write(*,*) "------------------------------------"
  write(*,*) "successful? ", burn_state_out % success
  write(*,*) "Completed burn to: ", burn_state_out % time, " seconds:"
  write(*,*) " - Hnuc = ", burn_state_out % e / dt
  write(*,*) " - integrated e = ", eos_state_in % e + energy
  write(*,*) " - EOS e(rho, T) = ", eos_state_out % e
  write(*,*) " - integrated/EOS percent diff. = ", 100.0e0_rt * (eos_state_in % e + energy - eos_state_out % e)/eos_state_out % e

  ! output burn type data
  !call write_burn_t(burn_state_out)

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

end subroutine burn_cell

subroutine write_burn_t(burnt)
  use network
  use burn_type_module

  use amrex_fort_module, only : rt => amrex_real
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
