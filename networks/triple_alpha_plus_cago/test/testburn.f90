program testburn

  use microphysics_type_module, only: rt, ZERO
  use microphysics_module
  use network
  use eos_module
  use actual_burner_module
  use runtime_init_module, only : runtime_pretty_print

  implicit none

  real(rt) :: dens, temp, dt
  real(rt), dimension(nspec) :: Xin
  type(burn_t) :: state_in, state_out 

  integer :: n
  integer :: of

  call microphysics_init()

  dens = 1.e7_rt
  temp = 367821127.2_rt

  Xin(ihe4)  = 0.8_rt
  Xin(ic12)  = 0.1_rt
  Xin(io16)  = 0.1_rt
  Xin(ife56) = 0.0_rt

  dt = 0.001_rt

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

  print *, 'rho_Hnuc: ', dens * (state_out % e - state_in % e) / dt

  print *, "number of RHS evaluations", state_out % n_rhs

  call microphysics_finalize()

end program testburn
