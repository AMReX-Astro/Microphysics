module tvode_rhs_module

  implicit none

contains

  SUBROUTINE FEX (NEQ, T, Y, YDOT, RPAR, IPAR)
    use bl_types, only: dp_t

    integer    :: NEQ, IPAR(:)
    real(dp_t) :: RPAR(:), T, Y(NEQ)
    real(dp_t), pointer :: YDOT(:)

    YDOT(1) = -.04D0*Y(1) + 1.D4*Y(2)*Y(3)
    YDOT(3) = 3.D7*Y(2)*Y(2)
    YDOT(2) = -YDOT(1) - YDOT(3)
    RETURN
  END SUBROUTINE FEX

  SUBROUTINE JEX (NEQ, T, Y, ML, MU, PD, NRPD, RPAR, IPAR)
    use bl_types, only: dp_t

    integer    :: NRPD, NEQ, ML, MU, IPAR(:)
    real(dp_t) :: PD(NRPD,NEQ), RPAR(:), T, Y(NEQ)

    PD(1,1) = -.04D0
    PD(1,2) = 1.D4*Y(3)
    PD(1,3) = 1.D4*Y(2)
    PD(2,1) = .04D0
    PD(2,3) = -PD(1,3)
    PD(3,2) = 6.D7*Y(2)
    PD(2,2) = -PD(1,2) - PD(3,2)
    RETURN
  END SUBROUTINE JEX

end module tvode_rhs_module
