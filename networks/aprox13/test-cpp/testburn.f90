subroutine do_burn() bind (C)

  use network
  use eos_module
  use burner_module
  use actual_burner_module

  implicit none

  type (burn_t) :: state_in, state_out

  double precision :: time = 0.0, dt = 2.5d-4

  type (eos_t) :: eos_state

  character (len=32) :: probin_file
  integer :: probin_pass(32)
  integer :: i
  double precision :: start, finish

  probin_file = "probin"
  do i = 1, len(trim(probin_file))
     probin_pass(i) = ichar(probin_file(i:i))
  enddo

  call runtime_init(probin_pass(1:len(trim(probin_file))), len(trim(probin_file)))

  call network_init()
  call burner_init()
  call eos_init()

  state_in % rho       = 19134047.319619529d0
  state_in % T         = 2091535024.9521415d0

  state_in % xn(ihe4)  = 5.7883378521362706d-005
  state_in % xn(ic12)  = 0.39011705806048469d0
  state_in % xn(io16)  = 0.45214212005419885d0 
  state_in % xn(ine20) = 0.12165391998665546d0
  state_in % xn(img24) = 3.4000133879978744d-002
  state_in % xn(isi28) = 2.0242732826176763d-003   
  state_in % xn(is32)  = 4.6061465566806610d-006   
  state_in % xn(iar36) = 5.2044648939936637d-009
  state_in % xn(ica40) = 2.5216462475548033d-012   
  state_in % xn(iti44) = 1.0000646156113865d-012
  state_in % xn(icr48) = 1.0000116722154411d-012   
  state_in % xn(ife52) = 1.0000093693211233d-012   
  state_in % xn(ini56) = 1.0000083281605682d-012

  print *, "rho_in: ", state_in % rho
  print *, "T_in: ", state_in % T
  print *, "X_in: ", state_in % xn

  ! We need an initial energy consistent with this
  ! temperature.

  call burn_to_eos(state_in, eos_state)
  call normalize_abundances(eos_state)
  call eos(eos_input_rt, eos_state)
  call eos_to_burn(eos_state, state_in)

  state_in % e = 0.0d0

  state_out = state_in

  call actual_burner(state_in, state_out, dt, time)

  print *, 'done!'

  print *, "rho_out: ", state_out % rho
  print *, "T_out: ", state_out % T
  print *, "X_out: ", state_out % xn

  print *, "Energy change: ", state_out % e - state_in % e

end subroutine do_burn
