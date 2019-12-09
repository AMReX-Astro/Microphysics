subroutine do_burn() bind (C)

  use network
  use eos_module
  use burner_module
  use actual_burner_module
  use microphysics_type_module, only: rt

  implicit none

  type (burn_t) :: state_in, state_out

  real(rt) :: time = 0.0_rt, dt = 2.5e-4_rt

  type (eos_t) :: eos_state

  character (len=32) :: probin_file
  integer :: probin_pass(32)
  integer :: i
  real(rt) :: start, finish

  probin_file = "probin"
  do i = 1, len(trim(probin_file))
     probin_pass(i) = ichar(probin_file(i:i))
  enddo

  call runtime_init(probin_pass(1:len(trim(probin_file))), len(trim(probin_file)))

  call network_init()
  call burner_init()
  call eos_init()

  state_in % rho       = 1.0e7_rt
  state_in % T         = 1.0e9_rt

  state_in % xn(1)     = 0.8e0_rt
  state_in % xn(2)     = 0.1e0_rt
  state_in % xn(3)     = 0.1e0_rt

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

  !$acc parallel copyin(state_in, dt, time) copy(state_out)

  call actual_burner(state_in, state_out, dt, time)

  !$acc end parallel

  print *, 'done!'

  print *, "rho_out: ", state_out % rho
  print *, "T_out: ", state_out % T
  print *, "X_out: ", state_out % xn

  print *, "Energy change: ", state_out % e - state_in % e

end subroutine do_burn
