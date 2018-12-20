subroutine f_rhs(n, t, y, ydot, rpar, ipar)

  use amrex_constants_module
  use network  

  implicit none

  ! our convention is that y(1:nspec) are the species
  integer :: n
  double precision :: y(n), ydot(n)

  double precision :: rpar
  integer :: ipar

  double precision :: t
  double precision :: xfueltmp

  xfueltmp = max(y(ifuel_),0.d0)

  ydot(ifuel_) = -xfueltmp*(ONE-xfueltmp)
  ydot(iash_)  =  xfueltmp*(ONE-xfueltmp)

end subroutine f_rhs


subroutine jac(neq, t, y, ml, mu, pd, nrpd, rpar, ipar)

  use amrex_constants_module
  use network

  implicit none

  integer         , intent(IN   ) :: neq, ml, mu, nrpd, ipar
  double precision, intent(IN   ) :: y(neq), rpar, t
  double precision, intent(  OUT) :: pd(neq,neq)

  ! we get the thermodynamic state through the burner_aux module -- we freeze
  ! these to the values are the top of the timestep to avoid costly

end subroutine jac

