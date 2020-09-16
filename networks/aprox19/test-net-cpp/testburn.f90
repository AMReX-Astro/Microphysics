subroutine do_burn() bind (C)

  use network
  use eos_module
  use eos_type_module
  use burner_module
  use actual_burner_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  type (burn_t) :: state_in, state_out

  real(rt)         :: time = 0.0_rt, dt = 1.0e-6_rt

  type (eos_t) :: eos_state

  state_in % rho    = 2.0e7_rt
  state_in % T      = 8.0e9_rt

  state_in % xn(:)  = ONE / nspec

  call burn_to_eos(state_in, eos_state)
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
