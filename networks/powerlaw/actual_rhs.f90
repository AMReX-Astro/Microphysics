module actual_rhs_module

  implicit none

contains

  subroutine actual_rhs(neq, time, y, ydot, rpar)

    use bl_types
    use bl_constants_module
    use network
    use vode_data
    use rpar_indices
    use extern_probin_module, only: f_act, T_burn_ref, rho_burn_ref, rtilde, nu

    implicit none

    integer          :: neq
    double precision :: time
    double precision :: y(neq), ydot(neq)
    double precision :: rpar(n_rpar_comps)

    double precision :: xfueltmp
    double precision :: dens, temp, rate

    ydot = ZERO

    xfueltmp = max(y(ifuel_) * aion(ifuel_), ZERO)
    dens     = rpar(irp_dens)
    temp     = y(net_itemp)

    ! Rate is expressed in mass fraction form
    
    if (temp < f_act * T_burn_ref) then
       rate = ZERO
    else
       rate = rtilde * (dens/rho_burn_ref) * xfueltmp**2 * (temp/T_burn_ref)**nu
    endif

    ydot(ifuel_)  = -rate
    ydot(iash_)   =  rate
    ydot(iinert_) =  ZERO

    ! Convert back to molar form

    ydot(1:nspec) = ydot(1:nspec) / aion
    
    ydot(net_ienuc) = -sum(ebin(:) * ydot(1:nspec))

    call temperature_rhs(neq, y, ydot, rpar)

  end subroutine actual_rhs



  ! At present the analytical Jacobian is not implemented.

  subroutine actual_jac(neq, time, y, pd, rpar)

    use network
    use rpar_indices
    use bl_constants_module, only: ZERO

    implicit none

    integer          :: neq
    double precision :: time
    double precision :: y(neq), pd(neq, neq)
    double precision :: rpar(n_rpar_comps)

    pd(:,:) = ZERO

  end subroutine actual_jac

end module actual_rhs_module
