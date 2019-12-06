program testburn

  use microphysics_type_module, only: rt, ZERO
  use network
  use eos_module
  use actual_burner_module
  use microphysics_module
  use extern_probin_module

  implicit none

  real(rt) :: dens, temp, dt
  real(rt), dimension(nspec) :: Xin
  type(burn_t) :: state_in, state_out

  call microphysics_init()

  dens = 2.6e9_rt
  temp = 6.e8_rt

  Xin(ic12)  = 0.5_rt
  Xin(io16)  = 0.5_rt
  Xin(iash)  = 0.0_rt

  dt = 0.06_rt


  print *, 'calling the burner...', nspec, nspec_evolve, neqs

  state_in % rho = dens
  state_in % T = temp
  state_in % e = ZERO
  state_in % xn(:) = Xin(:)

  jacobian = 2

  call actual_burner(state_in, state_out, dt, ZERO)

  print *, 'done!'

  print *, 'Xin:  ', state_in % xn(:)
  print *, 'Xout: ', state_out % xn(:)
  print *, 'rho_Hnuc: ', dens * (state_out % e - state_in % e) / dt
  print *, 'Hnuc (erg/g/s): ', (state_out % e - state_in % e) / dt

  call microphysics_finalize()

end program testburn
