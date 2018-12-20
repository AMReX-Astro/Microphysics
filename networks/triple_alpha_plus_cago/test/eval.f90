program eval

  use amrex_constants_module
  use amrex_fort_module, only : rt => amrex_real
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

  dens = 1.e7_rt
  temp = 367821127.2_rt

  Xin(ihe4)  = 0.8_rt
  Xin(ic12)  = 0.1_rt
  Xin(io16)  = 0.1_rt
  Xin(ife56) = 0.0_rt

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
