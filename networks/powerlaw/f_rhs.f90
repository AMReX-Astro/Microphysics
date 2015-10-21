subroutine f_rhs(n, t, y, ydot, rpar, ipar)

  use bl_types
  use bl_constants_module
  use network
  use rpar_indices
  use extern_probin_module
  
  implicit none

  ! our convention is that y(1:nspec) are the species (in the same
  ! order as defined in network.f90, and y(nspec+1) is the temperature
  integer :: n
  double precision :: y(n), ydot(n)

  double precision :: rpar(n_rpar_comps)
  integer :: ipar

  double precision :: t
  double precision :: xfueltmp

  double precision :: dens, temp, rate

  ! note that we don't really need to integrate this unless we 
  ! include a temperature-evolution equation in the mix.  Otherwise,
  ! it is analytic

  xfueltmp = max(y(ifuel_), ZERO)
  dens = rpar(irp_dens)
  temp = rpar(irp_temp)

  rate = rtilde*(dens/rho_burn_ref)*xfueltmp**2*(temp/T_burn_ref)**nu
  ydot(ifuel_) = -rate
  ydot(iash_)  =  rate
  ydot(iinert_)  =  ZERO

end subroutine f_rhs


subroutine jac(neq, t, y, ml, mu, pd, nrpd, rpar, ipar)

  use bl_types
  use bl_constants_module
  use network

  ! EOS calls

  implicit none

  integer        , intent(IN   ) :: neq, ml, mu, nrpd, ipar
  double precision, intent(IN   ) :: y(neq), rpar, t
  double precision, intent(  OUT) :: pd(neq,neq)

end subroutine jac

