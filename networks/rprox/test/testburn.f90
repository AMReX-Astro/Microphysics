program testburn

  use microphysics_type_module

  use network
  use eos_module
  use actual_burner_module
  use burn_type_module
  use microphysics_module
  use extern_probin_module

  implicit none

  real(rt) :: dens, temp, dt
  real(rt), dimension(nspec) :: Xin
  type(burn_t) :: state_in, state_out
  integer :: n

  call microphysics_init()

  dens = 1737938.5689184519e0_rt
  temp = 158489319.24611109e0_rt     

  Xin = [5.0000000000000003e-002_rt, &
         0.0000000000000000_rt, &
         0.0000000000000000_rt, &
         0.0000000000000000_rt, &
         0.0000000000000000_rt, &
         0.0000000000000000_rt, &
         0.0000000000000000_rt, &
         0.0000000000000000_rt, &
         0.25000000000000000_rt, &
         0.69999999999999996_rt]


  dt = 1.e-4_rt
  !use_timestep_estimator = .true.
  !scaling_method = 2
  !ode_scale_floor = 1.d-12
  burning_mode = 0

  print *, 'calling the burner...'

  jacobian = 2
  centered_diff_jac = .true.

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
