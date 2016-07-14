program testburn

  use bl_types
  use bl_constants_module
  use network
  use eos_module
  use actual_burner_module
  use microphysics_module

  implicit none

  real(kind=dp_t) :: dens, temp, dt
  real(kind=dp_t), dimension(nspec) :: Xin
  type(burn_t) :: state_in, state_out 

  call microphysics_init()

  dens = 2.6e9_dp_t
  temp = 6.e8_dp_t

  Xin(ic12_)  = 0.5_dp_t
  Xin(io16_)  = 0.5_dp_t
  Xin(iash_)  = 0.0_dp_t

  dt = 0.06_dp_t

  
  print *, 'calling the burner...'

  state_in % rho = dens
  state_in % T = temp
  state_in % e = ZERO
  state_in % xn(:) = Xin(:)

  call actual_burner(state_in, state_out, dt, ZERO)

  print *, 'done!'

  print *, 'Xin:  ', state_in % xn(:)
  print *, 'Xout: ', state_out % xn(:)
  print *, 'rho_Hnuc: ', dens * (state_out % e - state_in % e) / dt
  print *, 'Hnuc (erg/g/s): ', (state_out % e - state_in % e) / dt

  call microphysics_finalize()

end program testburn
