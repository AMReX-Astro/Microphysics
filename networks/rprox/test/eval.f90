program eval

  use amrex_constants_module
  use microphysics_type_module

  use network
  use eos_module
  use actual_burner_module
  use burn_type_module
  use microphysics_module
  use actual_network, only: short_spec_names
  use actual_rhs_module, only: actual_rhs
  
  implicit none

  real(rt) :: dens, temp
  real(rt), dimension(nspec) :: Xin
  type(burn_t) :: state_in
  integer :: n
  
  call microphysics_init()

  dens = 6.558e5_rt
  temp = 7.108e8_rt

  Xin(:) = ZERO
  Xin(io14) = 1.027e-3_rt
  Xin(io15) = 3.558e-3_rt
  Xin(if17) = 2.788e-2_rt
  Xin(is30) = 1.5735e-2_rt
  Xin(ihe4) = 0.2624e0_rt
  Xin(ih1) = 0.6894e0_rt

  print *, 'calling the burner RHS ...'

  state_in % rho = dens
  state_in % T = temp
  state_in % e = ZERO
  state_in % xn(:) = Xin(:)

  call actual_rhs(state_in)

  do n = 1, nspec
     print *, short_spec_names(n), state_in % xn(n), state_in % ydot(n)
  enddo

  call microphysics_finalize()

end program eval
