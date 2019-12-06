! The example RHS from VODE

subroutine f_rhs(neq, t, y, ydot)

  use microphysics_type_module, only: rt

  integer, intent(in) :: neq
  real(rt) t, y(neq), ydot(neq)

  ydot(1) = -.04e0_rt*y(1) + 1.e4_rt*y(2)*y(3)
  ydot(3) = 3.e7_rt*y(2)*y(2)
  ydot(2) = -ydot(1) - ydot(3)
  
end subroutine f_rhs
 
subroutine jac(neq, t, y, dfdy)

  use microphysics_type_module, only: rt

  integer, intent(in) :: neq
  real(rt) rpar, t, y(neq), dfdy(neq,neq)
  dfdy(:,:) = 0.0e0_rt

  dfdy(1,1) = -.04e0_rt
  dfdy(1,2) = 1.e4_rt*y(3)
  dfdy(1,3) = 1.e4_rt*y(2)
  dfdy(2,1) = .04e0_rt
  dfdy(2,3) = -dfdy(1,3)
  dfdy(3,2) = 6.e7_rt*y(2)
  dfdy(2,2) = -dfdy(1,2) - dfdy(3,2)

end subroutine jac
