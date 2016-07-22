program testburn

  use bl_types
  use bl_constants_module
  use microphysics_module
  use network
  use eos_module
  use actual_burner_module
  use extern_probin_module, only : jacobian
  use runtime_init_module, only : runtime_pretty_print

  implicit none

  real(kind=dp_t) :: dens, temp, dt
  real(kind=dp_t), dimension(nspec) :: Xin
  type(burn_t) :: state_in, state_out

  integer :: n
  integer :: of

  call microphysics_init()

  dens = 1.5e7_dp_t
  temp = 3.0e8_dp_t

  Xin(:) = ZERO
  Xin(ihe4)  = ONE
  Xin(ic12)  = ZERO

  dens = 100000000.d0
  temp = 2705847633.d0

  Xin(:) = ZERO
  Xin(ihe4) = 0.45d0
  Xin(ic12) = 0.2d0
  Xin(io16) = 0.2d0
  Xin(ine20) = 0.15d0

  dt = 0.0001_dp_t

  jacobian = 2

  open(newunit=of, file="testburn-params.out", status="replace", action="write")
  call runtime_pretty_print(of)
  close(unit=of)

  print *, 'calling the burner...'

  state_in % rho = dens
  state_in % T = temp
  state_in % e = ZERO
  state_in % xn(:) = Xin(:)

  call actual_burner(state_in, state_out, dt, ZERO)

  print *, 'done!'

1000 format(g25.15, 1x, g25.15)

  print *, 'Xin / Xout:'
  do n = 1, nspec
     print 1000, state_in % xn(n), state_out % xn(n)
  enddo

  print *, 'rho_Hnuc: ', dens * (state_out % e - state_in % e) /dt

  call microphysics_finalize()

end program testburn
