subroutine do_burn() bind (C)

  use network
  use actual_rhs_module, only: actual_rhs_init
  use eos_module
  use eos_type_module
  use burner_module
  use actual_burner_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  type (burn_t) :: state_in, state_out

  real(rt)         :: time = 0.0_rt, dt = 2.5e-4_rt

  type (eos_t) :: eos_state

  integer :: i
  real(rt)         :: start, finish


  state_in % rho       = 19134047.319619529e0_rt
  state_in % T         = 2091535024.9521415e0_rt

  state_in % xn(ihe4)  = 5.7883378521362706e-005_rt
  state_in % xn(ic12)  = 0.39011705806048469e0_rt
  state_in % xn(io16)  = 0.45214212005419885e0_rt 
  state_in % xn(ine20) = 0.12165391998665546e0_rt
  state_in % xn(img24) = 3.4000133879978744e-002_rt
  state_in % xn(isi28) = 2.0242732826176763e-003_rt   
  state_in % xn(is32)  = 4.6061465566806610e-006_rt   
  state_in % xn(iar36) = 5.2044648939936637e-009_rt
  state_in % xn(ica40) = 2.5216462475548033e-012_rt   
  state_in % xn(iti44) = 1.0000646156113865e-012_rt
  state_in % xn(icr48) = 1.0000116722154411e-012_rt   
  state_in % xn(ife52) = 1.0000093693211233e-012_rt   
  state_in % xn(ini56) = 1.0000083281605682e-012_rt

  print *, "rho_in: ", state_in % rho
  print *, "T_in: ", state_in % T
  print *, "X_in: ", state_in % xn

  ! We need an initial energy consistent with this
  ! temperature.

  call burn_to_eos(state_in, eos_state)
  call normalize_abundances(eos_state)
  call eos(eos_input_rt, eos_state)
  call eos_to_burn(eos_state, state_in)

  state_in % e = 0.0e0_rt

  state_out = state_in

  call actual_burner(state_in, state_out, dt, time)

  print *, 'done!'

  print *, "rho_out: ", state_out % rho
  print *, "T_out: ", state_out % T
  print *, "X_out: ", state_out % xn

  print *, "Energy change: ", state_out % e - state_in % e

end subroutine do_burn
