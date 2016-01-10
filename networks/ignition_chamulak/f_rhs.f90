subroutine f_rhs(n, time, y, ydot, rpar, ipar)

  use bl_types
  use bl_constants_module
  use network
  use eos_module
  use network_indices
  use rpar_indices
  use screening_module, only: screenz
  
  implicit none

  ! our convention is that y(1:nspec) are the species (in the same
  ! order as defined in network.f90, and y(nspec+1) is the temperature
  integer,          intent(in   ) :: n, ipar
  double precision, intent(in   ) :: y(n)
  double precision, intent(in   ) :: time
  double precision, intent(  out) :: ydot(n)
  double precision, intent(inout) :: rpar(*)

  double precision :: ymass(nspec)

  double precision :: dens, c_p, dhdX(nspec), X_O16
  double precision :: temp, T9, T9a, dT9dt, dT9adt

  double precision :: rate, dratedt
  double precision :: sc1212, dsc1212dt
  double precision :: xc12tmp

  double precision :: scratch, dscratchdt

  double precision :: a, b, dadt, dbdt

  double precision, parameter :: FIVE6TH = FIVE / SIX


  ! we freeze these to the values are the top of the timestep to avoid costly
  ! EOS calls
  dens = rpar(irp_dens)
  temp = y(nspec_advance+1)

  c_p     = rpar(irp_cp)
  dhdX(:) = rpar(irp_dhdX:irp_dhdX-1+nspec)
  X_O16   = rpar(irp_o16)

  ! compute the molar fractions -- needed for the screening
  ymass(ic12_) = y(1)/aion(ic12_)
  ymass(io16_) = X_O16/aion(io16_)
  ymass(iash_) = (ONE - y(1) - X_O16)/aion(iash_)


  ! call the screening routine
  call screenz(temp,dens,SIX,SIX,TWELVE,TWELVE,ymass,aion,zion,nspec,     &
               sc1212, dsc1212dt)


  ! compute some often used temperature constants
  T9     = temp/1.d9
  dT9dt  = ONE/1.d9
  T9a    = T9/(ONE + 0.0396d0*T9)
  dT9adt = (T9a / T9 - (T9a / (ONE + 0.0396d0*T9)) * 0.0396d0) * dT9dt

  ! compute the CF88 rate
  scratch    = T9a**THIRD
  dscratchdt = THIRD * T9a**(-TWO3RD) * dT9adt

  a       = 4.27d26*T9a**FIVE6TH*T9**(-1.5d0)
  dadt    = FIVE6TH * (a/T9a) * dT9adt - 1.5d0 * (a/T9) * dT9dt

  b       = dexp(-84.165d0/scratch - 2.12d-3*T9*T9*T9)
  dbdt    = (84.165d0 * dscratchdt/ scratch**TWO       &
             - THREE * 2.12d-3 * T9 * T9 * dT9dt) * b

  rate    = a *  b
  dratedt = dadt * b + a * dbdt

  ! The change in number density of C12 is
  ! d(n12)/dt = - M12_chamulak * 1/2 (n12)**2 <sigma v>
  !
  ! where <sigma v> is the average of the relative velocity times the
  ! cross section for the reaction, and the factor accounting for the
  ! total number of particle pairs has a 1/2 because we are
  ! considering a reaction involving identical particles (see Clayton
  ! p. 293).  Finally, the -M12_chamulak means that for each reaction,
  ! we lose M12_chamulak C12 nuclei (for a single rate, C12+C12,
  ! M12_chamulak would be 2.  In Chamulak et al. (2008), they say a
  ! value of 2.93 captures the energetics from a larger network
  !
  ! Switching over to mass fractions, using n = rho X N_A/A, where N_A is
  ! Avogadro's number, and A is the mass number of the nucleon, we get
  !
  ! d(X12)/dt = -M12_chamulak * 1/2 (X12)**2 rho N_A <sigma v> / A12
  !
  ! The quantity [N_A <sigma v>] is what is tabulated in Caughlin and Fowler.
  !
  ! we will always refer to the species by integer indices that come from
  ! the network module -- this makes things robust to a shuffling of the
  ! species ordering

  xc12tmp = max(y(ic12_),0.d0)
  ydot(ic12_) = -TWELFTH*HALF*M12_chamulak*dens*sc1212*rate*xc12tmp**2

  ! now compute the change in temperature, using the evolution equation
  ! dT/dt = -(1/c_p) sum_k (xi_k + q_k) omega_k
  !
  ! we make use of the fact that omega(Mg24) = - omega(C12), and that
  ! omega(O16) = 0 in our simplified burner
  ydot(nspec_advance+1) =  ( (dhdx(iash_) - dhdx(ic12_)) + &
                              get_ebin_value(dens) )*ydot(ic12_)/c_p

  ! for Mr. Jacobian
  rpar(irp_rate)      = rate
  rpar(irp_dratedt)   = dratedt
  rpar(irp_sc1212)    = sc1212
  rpar(irp_dsc1212dt) = dsc1212dt
  rpar(irp_xc12tmp)   = xc12tmp

end subroutine f_rhs


subroutine jac(neq, t, y, ml, mu, pd, nrpd, rpar, ipar)

  use bl_types
  use bl_constants_module
  use network
  use network_indices
  use rpar_indices

  implicit none

  integer        , intent(IN   ) :: neq, ml, mu, nrpd, ipar
  double precision, intent(IN   ) :: y(neq), rpar(*)
  double precision, intent(IN   ) :: t
  double precision, intent(  OUT) :: pd(neq,neq)

  double precision :: dens, c_p, dhdX(nspec), X_O16
  double precision :: rate, dratedt, scorr, dscorrdt, xc12tmp

  integer :: itemp

  dens    = rpar(irp_dens)
  c_p     = rpar(irp_cp)
  dhdX(:) = rpar(irp_dhdX:irp_dhdX-1+nspec)
  X_O16   = rpar(irp_o16)

  rate     = rpar(irp_rate)
  dratedt  = rpar(irp_dratedt)
  scorr    = rpar(irp_sc1212)
  dscorrdt = rpar(irp_dsc1212dt)
  xc12tmp  = rpar(irp_xc12tmp)


  ! initialize
  pd(:,:)  = ZERO

  itemp = nspec_advance + 1


  ! carbon jacobian elements
  pd(ic12_, ic12_) = -TWO*TWELFTH*M12_chamulak*HALF*dens*scorr*rate*xc12tmp


  ! add the temperature derivatives: df(y_i) / dT
  pd(ic12_,itemp) = -TWELFTH*M12_chamulak*HALF* &
                                     (dens*rate*xc12tmp**2*dscorrdt  &
                                    + dens*scorr*xc12tmp**2*dratedt)


  ! add the temperature jacobian elements df(T) / df(y)
  pd(itemp,ic12_) =  ( (dhdX(iash_) - dhdX(ic12_)) + &
                      get_ebin_value(dens) )*pd(ic12_,ic12_)/c_p


  ! add df(T) / dT
  pd(itemp,itemp) = ( (dhdX(iash_) - dhdX(ic12_)) + &
                       get_ebin_value(dens) )*pd(ic12_,itemp)/c_p


end subroutine jac
