program testburn

  use microphysics_type_module
  use network
  use microphysics_module
  use actual_burner_module

  implicit none

  real(rt) :: dens, temp, t, dt
  real(rt), dimension(nspec) :: Xin
  type(burn_t) :: state_in, state_out

  call microphysics_init()

  dens = 2.6e9_rt
  temp = 6.e8_rt

  Xin(ic12) = 0.5_rt
  Xin(io16) = 0.5_rt
  Xin(img24) = 0.0_rt

  t = 0.0_rt
  dt = 0.06_rt

  state_in % rho = dens
  state_in % T = temp
  state_in % e = ZERO
  state_in % xn(:) = Xin(:)

  call actual_burner(state_in, state_out, dt, t)

  print *, 'Xin:  ', state_in % xn(:)
  print *, 'Xout: ', state_out % xn(:)
  print *, 'rho Hnuc (erg/s/cm**3): ', dens * state_out % e / dt

  call microphysics_finalize()

end program testburn
