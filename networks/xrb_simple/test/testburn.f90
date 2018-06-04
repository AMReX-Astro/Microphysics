program testburn

  use bl_types
  use bl_constants_module
  use network
  use eos_module
  use actual_burner_module
  use burn_type_module
  use microphysics_module
  use extern_probin_module

  implicit none

  real(kind=dp_t) :: dens, temp, dt
  real(kind=dp_t), dimension(nspec) :: Xin
  type(burn_t) :: state_in, state_out
  integer :: n

  call microphysics_init()

  dens = 6.558e5_dp_t
  temp = 7.108e8_dp_t

  Xin(io14) = 1.027e-3_dp_t
  Xin(io15) = 3.558e-3_dp_t
  Xin(ine18) = 2.788e-2_dp_t
  Xin(isi25) = 1.5735e-2_dp_t
  Xin(ihe4) = 0.2624e0_dp_t
  Xin(ih1) = 0.6894e0_dp_t

  dt = 0.01_dp_t

  print *, 'calling the burner...'

  integrate_temperature = .false.

  state_in % rho = dens
  state_in % T = temp
  state_in % e = ZERO
  state_in % xn(:) = Xin(:)

  call actual_burner(state_in, state_out, dt, ZERO)

  print *, 'done!'

1000 format(g20.10, 1x, g20.10)

  print *, 'Xin / Xout:'
  do n = 1, nspec
     print 1000, state_in % xn(n), state_out % xn(n)
  enddo

  print *, 'Hnuc: ', (state_out % e - state_in % e) / dt

  call microphysics_finalize()

end program testburn
