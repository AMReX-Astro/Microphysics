program testburn

  use amrex_constants_module
  use amrex_fort_module, only : rt => amrex_real
  use microphysics_module
  use network
  use eos_module
  use actual_burner_module
  use extern_probin_module
  use runtime_init_module, only : runtime_pretty_print

  implicit none

  real(rt) :: dens, temp, dt
  real(rt), dimension(nspec) :: Xin
  type(burn_t) :: state_in, state_out

  integer :: n
  integer :: of

  call microphysics_init()

  dens =    10000.0e0_rt
  temp =    4999999999.9999990e0_rt
  Xin(:) = 1.e-10_rt
  Xin(ihe4) = 1.0_rt - (nspec-1)*1.e-10_dp_t

  dt = 0.001_rt

  jacobian = 2
  centered_diff_jac = .true.

  !p_age = 1
  !jac_age = 1

  !dt_min = 0.0

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
