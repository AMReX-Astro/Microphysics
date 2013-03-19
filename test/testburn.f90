program testburn

  use bl_types
  use bl_constants_module
  use network
  use eos_module
  use burner_module

  implicit none
  integer :: k
  real(kind=dp_t) :: dens, temp, dt, rho_Hnuc
  real(kind=dp_t), dimension(nspec) :: Xin, Xout, rho_omegadot

  integer :: ihe4, ic12

  call network_init()
  call eos_init()

  ihe4  = network_species_index("he4")
  ic12  = network_species_index("c12")

  dens = 1.5e7_dp_t
  temp = 3.0e8_dp_t

  Xin(:) = ZERO
  Xin(ihe4)  = ONE
  Xin(ic12)  = ZERO

  dt = 0.01_dp_t

  
  print *, 'calling the burner...'


  call burner(dens, temp, Xin, dt, Xout, rho_omegadot, rho_Hnuc)

  print *, 'done!'

  print *, 'Xin:  ', Xin
  print *, 'Xout: ', Xout
  print *, 'rho_omegadot: ', rho_omegadot
  print *, 'rho_Hnuc: ', rho_Hnuc

end program testburn
