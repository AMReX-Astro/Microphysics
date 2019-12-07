subroutine f_rhs(n, t, y, ydot, rpar, ipar)

  use amrex_constants_module
  use network  

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  ! our convention is that y(1:nspec) are the species
  integer :: n
  real(rt)         :: y(n), ydot(n)

  real(rt)         :: rpar
  integer :: ipar

  real(rt)         :: t
  real(rt)         :: xfueltmp

  xfueltmp = max(y(ifuel_),0.e0_rt)

  ydot(ifuel_) = -xfueltmp*(ONE-xfueltmp)
  ydot(iash_)  =  xfueltmp*(ONE-xfueltmp)

end subroutine f_rhs


subroutine jac(neq, t, y, ml, mu, pd, nrpd, rpar, ipar)

  use amrex_constants_module
  use network

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer         , intent(IN   ) :: neq, ml, mu, nrpd, ipar
  real(rt)        , intent(IN   ) :: y(neq), rpar, t
  real(rt)        , intent(  OUT) :: pd(neq,neq)

  ! we get the thermodynamic state through the burner_aux module -- we freeze
  ! these to the values are the top of the timestep to avoid costly

end subroutine jac

