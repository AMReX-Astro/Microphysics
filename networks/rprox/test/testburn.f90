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

  dens = 1737938.5689184519
  temp = 158489319.24611109     

  Xin = [5.0000000000000003E-002, &
         0.0000000000000000, &
         0.0000000000000000, &
         0.0000000000000000, &
         0.0000000000000000, &
         0.0000000000000000, &
         0.0000000000000000, &
         0.0000000000000000, &
         0.25000000000000000, &
         0.69999999999999996]


  dt = 1.e-4_dp_t
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
