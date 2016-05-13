program testburn

  use bl_types
  use bl_constants_module
  use bl_error_module
  use network
  use eos_module
  use burner_module
  use burner_aux_module

  implicit none

  real(kind=dp_t) :: dens, temp
  real(kind=dp_t), dimension(nspec) :: Xin
  real(kind=dp_t), dimension(nspec+1) :: y, ydot
  type(burn_t) :: state
  
  integer :: ic12, ine20, ina23
  integer :: n

  call network_init()
  call eos_init()

  ic12 = network_species_index("carbon-12")
  ine20 = network_species_index("neon-20")
  ina23 = network_species_index("sodium-23")

  if (ic12 < 0 .or. ine20 < 0 .or. ina23 < 0) then
     call bl_error("ERROR: species index not defined")
  endif

  state%rho = 1.0e7_dp_t
  state%T = 1.0e9_dp_t
  state%y_e = 0.5_dp_t
  state%xn(:) = 0.0_dp_t
  state%xn(ic12) = 1.0_dp_t

  print *, 'evaluating the RHS...'

  call actual_rhs(state)
  
  print *, 'done!'
  print *, 'xn:  ', xn
  print *, 'rhs:  ', ydot

end program testburn
