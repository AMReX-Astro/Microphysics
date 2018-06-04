program testburn

  use bl_types
  use bl_constants_module
  use network
  use microphysics_module
  use actual_burner_module

  implicit none

  real(kind=dp_t) :: dens, temp, t, dt
  real(kind=dp_t), dimension(nspec) :: Xin
  type(burn_t) :: state_in, state_out

  call microphysics_init()

  dens = 2.6e9_dp_t
  temp = 6.e8_dp_t

  Xin(ic12) = 0.5_dp_t
  Xin(io16) = 0.5_dp_t
  Xin(img24) = 0.0_dp_t

  t = 0.0_dp_t
  dt = 0.06_dp_t

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
