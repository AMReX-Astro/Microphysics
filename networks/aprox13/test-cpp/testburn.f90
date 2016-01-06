subroutine do_burn() bind (C)

  use network
  use network_indices
  use eos_module
  use burner_module
  use actual_burner_module

  implicit none

  type (eos_t) :: state_in, state_out

  call network_init()
  call eos_init()

  state_in % rho = 1.4210128725795556d7
  state_in % T = 4.7053584415339775d9

  state_in % xn(ihe4) = 3.4573088050620371d-3
  state_in % xn(ic12) = 2.8224859252252550d-5
  state_in % xn(io16) = 1.3939441952623305d-4
  state_in % xn(ine20) = 5.7447866213696762d-7
  state_in % xn(img24) = 2.7057281109966992d-4
  state_in % xn(isi28) = 0.35451909718671237d0
  state_in % xn(is32) = 0.25344593450184338d0
  state_in % xn(iar36) = 8.3977223489996505d-2
  state_in % xn(ica40) = 7.4009021542773307d-2
  state_in % xn(iti44) = 5.8498004479900112d-4
  state_in % xn(icr48) = 2.5591436360638212d-3
  state_in % xn(ife52) = 2.1193774212954306d-2
  state_in % xn(ini56) = 0.20581475041307201d0

  call normalize_abundances(state_in)

  call eos(eos_input_rt, state_in)

  call burner(state_in, state_out, 1.0d-6, 0.0d0)

  print *, 'done!'

  print *, 'Xout: ', state_out%xn(:)

end subroutine do_burn

