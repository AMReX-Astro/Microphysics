! Evaluate the reaction rates for given thermodynamic conditions

program test_rates

  use amrex_error_module
  use amrex_constants_module
  use amrex_fort_module, only : rt => amrex_real

  use fabio_module, only: fabio_mkdir
  use probin_module, only: run_prefix, small_temp, small_dens
  use runtime_init_module
  use extern_probin_module
  use burn_type_module
  use actual_burner_module
  use actual_rhs_module, only: rate_eval_t, evaluate_rates
  use microphysics_module
  use network
  use build_info_module

  implicit none

  type (burn_t)       :: burn_state
  type (rate_eval_t), allocatable :: rate_state(:)

  integer             :: i, irho, itemp, istate, density_npts, temperature_npts, numpoints

  character (len=256) :: params_file
  integer             :: params_file_unit

  real(rt)     :: density_lo, density_hi, delta_density, &
                  temperature_lo, temperature_hi, delta_temperature, &
                  massfractions(nspec)

  namelist /cellparams/ density_lo, density_hi, density_npts, &
                        temperature_lo, temperature_hi, temperature_npts, &
                        massfractions

  ! runtime
  call runtime_init(.true.)

  ! microphysics
  call microphysics_init(small_temp=small_temp, small_dens=small_dens)

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

  ! Calculate grid deltas
  delta_density = (density_hi - density_lo)/(density_npts - 1)
  delta_temperature = (temperature_hi - temperature_lo)/(temperature_npts - 1)

  ! Echo initial conditions at burn and fill burn state input
  write(*,*) 'State Density Range (g/cm^3): ', density_lo, ' - ', density_hi
  write(*,*) 'Number of density points: ', density_npts
  write(*,*) 'Delta Density (g/cm^3): ', delta_density
  write(*,*) ''
  write(*,*) 'State Temperature Range (K): ', temperature_lo, ' - ', temperature_hi
  write(*,*) 'Number of temperature points: ', temperature_npts
  write(*,*) 'Delta Temperature (K): ', delta_temperature  
  write(*,*) ''
  do i = 1, nspec
     write(*,*) 'Mass Fraction (', short_spec_names(i), '): ', massfractions(i)
  end do

  ! Allocate memory for rate evaluations
  numpoints = density_npts * temperature_npts
  write(*,*) 'Number of total points in rho, T: ', numpoints
  allocate(rate_state(numpoints))

  i = 1
  do irho = 1, density_npts
     do itemp = 1, temperature_npts

        burn_state % T   = temperature_lo + delta_temperature * (itemp - 1)
        burn_state % rho = density_lo + delta_density * (irho - 1)
        burn_state % xn(:) = massfractions(:)

        ! normalize -- just in case
        !call normalize_abundances_burn(burn_state)

        ! Initialize initial energy to zero
        burn_state % e = ZERO

        ! Evaluate rates
        call evaluate_rates(burn_state, rate_state(i))

        i = i + 1
     end do
  end do

  call microphysics_finalize()

  call write_rate_grid()

contains

  subroutine write_rate_grid
    use network, only: nrates, nrat_tabular, i_rate, i_drate_dt, i_scor, i_dscor_dt, nspec

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    character(len=256) :: fname
    integer :: file_unit, j, jstart, jend
    integer, parameter :: jentries = 2 + 5*nrates + 2*nrat_tabular
    character(len=50) :: VDPFMT = ''
    real(rt)  :: output_vector(jentries)

    fname = trim(run_prefix) // "_rates"
    open(newunit=file_unit, file=fname, action='WRITE')

    write(unit=file_unit, fmt=*) '! Evaluated Rates'
    write(unit=file_unit, fmt=*) '! The following fields are of size Number of Rates: '
    write(unit=file_unit, fmt=*) '!  rates, drates_dt, screening, dscreening_dt, screened_rates'
    write(unit=file_unit, fmt=*) '! The following fields are of size Number of Tabular Rates: '
    write(unit=file_unit, fmt=*) '!  tabular_dqweak, tabular_epart'

    write(unit=file_unit, fmt=*) 'Number of Rates: ', nrates
    write(unit=file_unit, fmt=*) 'Number of Tabular Rates: ', nrat_tabular

    write(unit=file_unit, fmt=*) 'Mass Fractions: '
    write(VDPFMT, '("(", I0, "E30.16E5", ")")') nspec
    write(unit=file_unit, fmt=VDPFMT) (massfractions(i), i = 1, nspec)

    write(VDPFMT, '("(", I0, "A18", ")")') 9
    write(unit=file_unit, fmt=VDPFMT) 'rho', 'T', 'rates', 'drates_dt', &
                                      'screening', 'dscreening_dt', &
                                      'screened_rates', 'tabular_dqweak', 'tabular_epart'

    write(VDPFMT, '("(", I0, "E30.16E5", ")")') jentries
    i = 1
    do irho = 1, density_npts
       do itemp = 1, temperature_npts
          output_vector(1) = density_lo + delta_density * (irho - 1)
          output_vector(2) = temperature_lo + delta_temperature * (itemp - 1)

          jstart = 3
          jend = jstart + nrates - 1
          output_vector(jstart:jend) = rate_state(i) % unscreened_rates(i_rate, :)

          jstart = jend + 1
          jend = jstart + nrates - 1
          output_vector(jstart:jend) = rate_state(i) % unscreened_rates(i_drate_dt, :)

          jstart = jend + 1
          jend = jstart + nrates - 1
          output_vector(jstart:jend) = rate_state(i) % unscreened_rates(i_scor, :)

          jstart = jend + 1
          jend = jstart + nrates - 1
          output_vector(jstart:jend) = rate_state(i) % unscreened_rates(i_dscor_dt, :)

          jstart = jend + 1
          jend = jstart + nrates - 1
          output_vector(jstart:jend) = rate_state(i) % screened_rates(:)

          jstart = jend + 1
          jend = jstart + nrat_tabular - 1
          output_vector(jstart:jend) = rate_state(i) % dqweak(:)

          jstart = jend + 1
          jend = jstart + nrat_tabular - 1
          output_vector(jstart:jend) = rate_state(i) % epart(:)

          write(unit=file_unit, fmt=VDPFMT) (output_vector(j), j = 1, jentries)
          
          i = i + 1
       end do
    end do
    close(unit=file_unit)
  end subroutine write_rate_grid
  
end program test_rates

subroutine write_rate_t(fname, rate_state)
  use network, only: nrates, nrat_tabular, i_rate, i_drate_dt, i_scor, i_dscor_dt
  use actual_rhs_module, only: rate_eval_t

  use amrex_fort_module, only : rt => amrex_real
  implicit none
  
  ! Writes contents of rate_eval_t type rate_state to file named fname
  character(len=256), intent(in) :: fname
  type(rate_eval_t), intent(in)  :: rate_state
  integer, parameter :: file_unit = 10
  character(len=50) :: VDPFMT = ''
  
  integer :: i

  open(unit=file_unit, file=fname, action='WRITE')
  write(unit=file_unit, fmt=*) '! Rate Eval Type Data'

  write(unit=file_unit, fmt=*) 'unscreened rates:'
  write(VDPFMT, '("(", I0, "E30.16E5", ")")') nrates
  write(unit=file_unit, fmt=VDPFMT) (rate_state % unscreened_rates(i_rate, i), i = 1, nrates)
  
  write(unit=file_unit, fmt=*) 'unscreened drates_dt:'
  write(VDPFMT, '("(", I0, "E30.16E5", ")")') nrates
  write(unit=file_unit, fmt=VDPFMT) (rate_state % unscreened_rates(i_drate_dt, i), i = 1, nrates)

  write(unit=file_unit, fmt=*) 'screening factors:'
  write(VDPFMT, '("(", I0, "E30.16E5", ")")') nrates
  write(unit=file_unit, fmt=VDPFMT) (rate_state % unscreened_rates(i_scor, i), i = 1, nrates)
  
  write(unit=file_unit, fmt=*) 'screening dfactors_dt:'
  write(VDPFMT, '("(", I0, "E30.16E5", ")")') nrates
  write(unit=file_unit, fmt=VDPFMT) (rate_state % unscreened_rates(i_dscor_dt, i), i = 1, nrates)

  write(unit=file_unit, fmt=*) 'screened rates:'
  write(VDPFMT, '("(", I0, "E30.16E5", ")")') nrates
  write(unit=file_unit, fmt=VDPFMT) (rate_state % screened_rates(i), i = 1, nrates)

  write(unit=file_unit, fmt=*) '(tabular) weak rate dq:'
  write(VDPFMT, '("(", I0, "E30.16E5", ")")') nrat_tabular
  write(unit=file_unit, fmt=VDPFMT) (rate_state % dqweak(i), i = 1, nrat_tabular)

  write(unit=file_unit, fmt=*) '(tabular) weak rate eparticle:'
  write(VDPFMT, '("(", I0, "E30.16E5", ")")') nrat_tabular
  write(unit=file_unit, fmt=VDPFMT) (rate_state % epart(i), i = 1, nrat_tabular)

  close(unit=file_unit)
end subroutine write_rate_t
