program testburn

  use bl_types
  use bl_constants_module
  use network
  use eos_module
  use eos_type_module
  use burner_module
  use rpar_indices
  
  implicit none

  real(kind=dp_t) :: dens, temp
  real(kind=dp_t), dimension(nspec) :: Xin
  real(kind=dp_t), dimension(nspec+1) :: y, ydot
  real(kind=dp_t) :: enucdot

  real(kind=dp_t), allocatable :: rpar(:)
  integer :: ipar

  real(kind=dp_t) :: c_p, dhdx(nspec), T_eos, dT_crit

  integer :: n

  integer :: ihe4, ic12, io16, ife56

  type(eos_t) :: eos_state
  
  call network_init()
  call eos_init()

  dens = 2.e9_dp_t
  temp = 7.e8_dp_t

  ihe4  = network_species_index("helium-4")
  ic12  = network_species_index("carbon-12")
  io16  = network_species_index("oxygen-16")
  ife56 = network_species_index("iron-56")

  Xin(ihe4)  = 0.5_dp_t
  Xin(ic12)  = 0.25_dp_t
  Xin(io16)  = 0.25_dp_t
  Xin(ife56) = 0.0_dp_t


  eos_state%rho = dens
  eos_state%T = temp
  eos_state%xn(:) = Xin(:)

  call eos(eos_input_rt, eos_state)

  print *, 'evaluating the RHS...'

  ! load the state
  y(1) = Xin(ihe4)
  y(2) = Xin(ic12)
  y(3) = Xin(io16)
  y(4) = Xin(ife56)
  y(5) = temp

  ! load the rpar array
  allocate(rpar(n_rpar_comps))
  
  rpar(irp_dens) = dens
  rpar(irp_Teos) = temp
  rpar(irp_Tcrit) = 0.01d0
  rpar(irp_dhdX:irp_dhdX-1+nspec) = eos_state%dhdX(:)
  rpar(irp_cp) = eos_state%cp
  rpar(irp_Y56) = Xin(ife56)/aion(ife56)
  
  call f_rhs(nspec+1, ZERO, y, ydot, rpar, ipar)


  print *, 'done!'

  print *, 'Xin:  ', Xin
  print *, 'rhs:  ', ydot

  ! compute the energy release/s (erg/g/s)
  enucdot = ZERO
  do n = 1, nspec
     enucdot = enucdot - ebin(n)*ydot(n)
  enddo
  print *, 'enucdot = ', enucdot

end program testburn
