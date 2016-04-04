subroutine do_burn() bind (C)

  use network
  use eos_module
  use burner_module
  use actual_burner_module

  implicit none

  type (burn_t) :: state_in, state_out

  double precision :: time = 0.0, dt = 1.25d-3

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

  state_in % rho       = 8.3476247460418558d6
  state_in % T         = 1.1253503737951503d9

  state_in % xn(ihe4)  = 4.4384367964996223d-9
  state_in % xn(ic12)  = 4.9936731456374661d-1
  state_in % xn(io16)  = 4.9957381026327208d-1
  state_in % xn(ine20) = 1.0581542663100497d-3
  state_in % xn(img24) = 7.1617670399329966d-7
  state_in % xn(isi28) = 2.8453065131927444d-10
  state_in % xn(is32)  = 1.0000113624896530d-12
  state_in % xn(iar36) = 1.0000000000854599d-12
  state_in % xn(ica40) = 1.0000000000052922d-12
  state_in % xn(iti44) = 1.0000000000002189d-12
  state_in % xn(icr48) = 9.9999999999999998d-13
  state_in % xn(ife52) = 1.0000000000000081d-12
  state_in % xn(ini56) = 1.0000000000000010d-12


  print *, "rho_in: ", state_in % rho
  print *, "T_in: ", state_in % T
  print *, "X_in: ", state_in % xn

  ! We need an initial energy consistent with this
  ! temperature.

  call burn_to_eos(state_in, eos_state)
  call normalize_abundances(eos_state)
  call eos(eos_input_rt, eos_state)
  call eos_to_burn(eos_state, state_in)

  state_out = state_in

  call actual_burner(state_in, state_out, dt, time)

  print *, 'done!'

  print *, "rho_out: ", state_out % rho
  print *, "T_out: ", state_out % T
  print *, "X_out: ", state_out % xn

  print *, "Energy change: ", state_out % e - state_in % e

end subroutine do_burn

