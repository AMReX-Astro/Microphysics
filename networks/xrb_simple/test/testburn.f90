program testburn

  use eos_module
!  use eos_type_module
  use network
  use burner_module
  use network_indices
  use bl_types
  use bl_constants_module, only: ZERO

  implicit none

  type(eos_t) :: state_in, state_out

  real(kind=dp_t) :: dens, temp, dt, rho_Hnuc
  real(kind=dp_t), dimension(nspec) :: Xin, Xout, rho_omegadot

  call network_init()
  call eos_init()

  dens = 6.558e5_dp_t
  temp = 7.108e8_dp_t

  Xin(io14) = 1.027e-3_dp_t
  Xin(io15) = 3.558e-3_dp_t
  Xin(ine18) = 2.788e-2_dp_t
  Xin(isi25) = 1.5735e-2_dp_t
  Xin(ihe4) = 0.2624e0_dp_t
  Xin(ih1) = 0.6894e0_dp_t

  dt = 0.01_dp_t

  Xout = ZERO

  print *, 'calling the burner...'

  call burner(dens, temp, Xin, dt, Xout, rho_omegadot, rho_Hnuc)

  print *, 'done!'

  print *, 'Xin:      ', Xin
  print *, 'sum:      ', sum(Xin)
  print *, 'Xout:     ', Xout
  print *, 'sum:      ', sum(Xout)
  print *, 'rho_Hnuc: ', rho_Hnuc
  print *, 'enuc:     ', rho_Hnuc*dt/dens

end program testburn
