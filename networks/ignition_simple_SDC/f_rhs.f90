subroutine f_rhs(n, t, y, ydot, rpar, ipar)

  use bl_types
  use bl_constants_module
  use network
  use bl_error_module
  
  use burner_aux_module, only : sdc_rhoX_pass, sdc_rhoh_pass, p0_pass

  implicit none

  ! our convention is that y(1:nspec) are the density-weighted species
  ! (in the same order as defined in network.f90, and y(nspec+1) is
  ! (rho h).
  integer :: n
  real(kind=dp_t) :: y(n), ydot(n)
  real(kind=dp_t) :: rpar
  integer :: ipar

  real(kind=dp_t) :: t

  real(kind=dp_t) :: rho_Hnuc

  integer, save :: ic12, io16, img24

  logical, save :: firstCall = .true.

  if (firstCall) then
     ic12 = network_species_index("carbon-12")
     io16 = network_species_index("oxygen-16")
     img24 = network_species_index("magnesium-24")

     firstCall = .false.
  end if

  ! get the RHS
  call f_rhs_instantaneous_reaction_rates(n, t, y, ydot, rho_Hnuc, rpar, ipar)

  ! now make ydots refer to rhoX and include the sdc sources
  ydot(ic12)  = ydot(ic12)  + sdc_rhoX_pass(ic12)
  ydot(io16)  = ydot(io16)  + sdc_rhoX_pass(io16)
  ydot(img24) = ydot(img24) + sdc_rhoX_pass(img24)
  ydot(nspec+1) = rho_Hnuc + sdc_rhoh_pass

  return

end subroutine f_rhs

subroutine jac(neq, t, y, ml, mu, pd, nrpd, rpar, ipar)

  use bl_types
  use bl_constants_module
  use network

  ! we get the thermodynamic state through the burner_aux module -- we freeze
  ! these to the values are the top of the timestep to avoid costly
  ! EOS calls

  implicit none

  integer        , intent(IN   ) :: neq, ml, mu, nrpd, ipar
  real(kind=dp_t), intent(IN   ) :: y(neq), rpar, t
  real(kind=dp_t), intent(  OUT) :: pd(neq,neq)

  ! initialize
  pd(:,:)  = ZERO

  return
end subroutine jac

