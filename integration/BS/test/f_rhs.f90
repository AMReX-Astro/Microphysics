! The example RHS from VODE

subroutine f_rhs(neq, t, y, ydot)
  integer, intent(in) :: neq
  double precision t, y(neq), ydot(neq)

  ydot(1) = -.04d0*y(1) + 1.d4*y(2)*y(3)
  ydot(3) = 3.d7*y(2)*y(2)
  ydot(2) = -ydot(1) - ydot(3)
  
end subroutine f_rhs
 
subroutine jac(neq, t, y, dfdy)
  integer, intent(in) :: neq
  double precision rpar, t, y(neq), dfdy(neq,neq)

  dfdy(1,1) = -.04d0
  dfdy(1,2) = 1.d4*y(3)
  dfdy(1,3) = 1.d4*y(2)
  dfdy(2,1) = .04d0
  dfdy(2,3) = -dfdy(1,3)
  dfdy(3,2) = 6.d7*y(2)
  dfdy(2,2) = -dfdy(1,2) - dfdy(3,2)

end subroutine jac
