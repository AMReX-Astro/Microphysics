subroutine do_burn() bind (C)

  use network
  use eos_module
  use burner_module
  use actual_burner_module

  implicit none

  type (burn_t) :: state_in, state_out

  double precision :: time = 0.0, dt = 3.125d-5

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

  state_in % rho       = 2.2933907037204478d7
  state_in % T         = 5.3757726209940176d9

  state_in % xn(ihe4)  = 5.8278210351460392d-2
  state_in % xn(ic12)  = 6.8334315299861196d-7
  state_in % xn(io16)  = 2.4687369464161301d-6
  state_in % xn(ine20) = 4.0170490921084967d-8
  state_in % xn(img24) = 1.1540561661029377d-5
  state_in % xn(isi28) = 1.3242394917027779d-2
  state_in % xn(is32)  = 2.0543812858533062d-2
  state_in % xn(iar36) = 1.5796293927450564d-2
  state_in % xn(ica40) = 2.7812814474560826d-2
  state_in % xn(iti44) = 7.6968156633668591d-4
  state_in % xn(icr48) = 5.2978509402914617d-3
  state_in % xn(ife52) = 6.1483728449529436d-2
  state_in % xn(ini56) = 7.9676047970255859d-1

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
