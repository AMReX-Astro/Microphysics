subroutine do_burn() bind (C)

  use network
  use eos_module
  use burner_module
  use burn_type_module
  use actual_burner_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  type (burn_t) :: state_in, state_out

  real(rt)         :: time = 0.0, dt = 1.0e-6_rt

  type (eos_t) :: eos_state

  character (len=32) :: probin_file
  integer :: probin_pass(32)
  integer :: i

  probin_file = "probin"
  do i = 1, len(trim(probin_file))
     probin_pass(i) = ichar(probin_file(i:i))
  enddo

  call runtime_init(probin_pass(1:len(trim(probin_file))), len(trim(probin_file)))

  call network_init()
  call burner_init()
  call eos_init()

  state_in % rho       = 1.4311401611205835e7_rt
  state_in % T         = 4.6993994016410122e9_rt

  state_in % xn(ihe4)  = 4.2717633762309063e-3_rt
  state_in % xn(ic12)  = 2.4502021307478711e-5_rt
  state_in % xn(io16)  = 1.2059146851610723e-4_rt
  state_in % xn(ine20) = 5.4419551339421394e-7_rt
  state_in % xn(img24) = 2.5178594678377961e-4_rt
  state_in % xn(isi28) = 3.5998829467937532e-1_rt
  state_in % xn(ini56) = 1.6962109761781674e-1_rt

  ! Add the remainder to iron
  state_in % xn(ini56) = state_in % xn(ini56) + ONE - sum(state_in % xn)

  call burn_to_eos(state_in, eos_state)
  call normalize_abundances(eos_state)
  call eos(eos_input_rt, eos_state)
  call eos_to_burn(eos_state, state_in)

  print *, "rho_in: ", state_in % rho
  print *, "T_in: ", state_in % T
  print *, "X_in: ", state_in % xn

  state_out = state_in
  
  call actual_burner(state_in, state_out, dt, time)
  
  print *, 'done!'

  print *, "rho_out: ", state_out % rho
  print *, "T_out: ", state_out % T
  print *, "X_out: ", state_out % xn

  print *, "Energy change: ", state_out % e - state_in % e

end subroutine do_burn

