! see top level interfaces file for description of subroutine function
!
! Dean Townsley, Alan Calder 2006-8
!
! original interpolation kernel by Alan Calder 2006


subroutine NSE_finalAtPres(qbar_nse,sumyi_nse,approxtemp,edot,Yedot, Ye, pres, hmq)
  
  use NSE_data, ONLY: nhmq, hmq_grid, lpres_grid, dhmq, dlpres, nlpres, &
       & p_dYe, p_Ye_grid, p_nYe, p_ltemp_tab, p_lqbar_tab, p_ledot_tab, &
       & p_lmYedot_tab, p_lAbar_tab

  implicit none
  
  real, intent(IN)    :: Ye, pres, hmq
  real, intent(OUT)   :: qbar_nse,sumyi_nse,approxtemp,edot,Yedot
  
  integer :: hmq_a, lpres_a, Ye_a
  
  real :: lpres, hmq_v
  real :: te111, te211, te121, te112, te221, te212, te122, te222
  real :: qb111, qb211, qb121, qb112, qb221, qb212, qb122, qb222
  real :: ed111, ed211, ed121, ed112, ed221, ed212, ed122, ed222
  real :: yd111, yd211, yd121, yd112, yd221, yd212, yd122, yd222
  real :: ab111, ab211, ab121, ab112, ab221, ab212, ab122, ab222
  
  real :: c1, c2, c3, abar
  
  
  !write (6,*) 'working at ye, pres, hmq', Ye, pres, hmq

  !! find the location in the table grid
  !! remember that we use the log of pressure
  
  hmq_v = hmq
  lpres = log10(pres) 
  
  ! we want to truncate at the "bottom" edges (in pressure and hmq)
  !instead of extrapolate.  Should be pure nickel down here.

  if ( hmq_v < hmq_grid(nhmq)   )   hmq_v = hmq_grid(nhmq)
  if ( lpres < lpres_grid(1) ) lpres = lpres_grid(1)

  ! find v_a  such that v_grid(v_a) <= v < v_grid(v_a+1)
  hmq_a   = floor((hmq_v - hmq_grid(1))/dhmq)+1
  lpres_a = floor((lpres - lpres_grid(1))/dlpres)+1
  Ye_a   = floor((Ye - p_Ye_grid(1))/p_dYe)+1

  ! treat edges (allows extrapolation for cases not handled above)
  if (hmq_a < 1) hmq_a = 1
  if (hmq_a >= nhmq) hmq_a = nhmq - 1
  ! lower limit accounted for above
  if (lpres_a >= nlpres) lpres_a = nlpres - 1

  if (Ye_a >= p_nYe) Ye_a = p_nYe - 1
  if (Ye_a < 1) Ye_a = 1

  !! calculate the coefficients
  
  c1 = ( hmq_v - hmq_grid(hmq_a) )/                &
       ( hmq_grid(hmq_a+1) - hmq_grid(hmq_a) )
  c2 = ( lpres - lpres_grid(lpres_a)) /              &
       ( lpres_grid(lpres_a+1) - lpres_grid(lpres_a) )
  c3 = ( Ye - p_Ye_grid(Ye_a) )/                  &
       ( p_Ye_grid(Ye_a+1) - p_Ye_grid(Ye_a) )

  !! build the local cubes
  te111 = p_ltemp_tab(hmq_a,lpres_a,Ye_a)
  te211 = p_ltemp_tab(hmq_a+1,lpres_a,Ye_a)
  te121 = p_ltemp_tab(hmq_a,lpres_a+1,Ye_a)
  te112 = p_ltemp_tab(hmq_a,lpres_a,Ye_a+1)
  te221 = p_ltemp_tab(hmq_a+1,lpres_a+1,Ye_a)
  te212 = p_ltemp_tab(hmq_a+1,lpres_a,Ye_a+1)
  te122 = p_ltemp_tab(hmq_a,lpres_a+1,Ye_a+1)
  te222 = p_ltemp_tab(hmq_a+1,lpres_a+1,Ye_a+1)
  
  qb111 = p_lqbar_tab(hmq_a,lpres_a,Ye_a)
  qb211 = p_lqbar_tab(hmq_a+1,lpres_a,Ye_a)
  qb121 = p_lqbar_tab(hmq_a,lpres_a+1,Ye_a)
  qb112 = p_lqbar_tab(hmq_a,lpres_a,Ye_a+1)
  qb221 = p_lqbar_tab(hmq_a+1,lpres_a+1,Ye_a)
  qb212 = p_lqbar_tab(hmq_a+1,lpres_a,Ye_a+1)
  qb122 = p_lqbar_tab(hmq_a,lpres_a+1,Ye_a+1)
  qb222 = p_lqbar_tab(hmq_a+1,lpres_a+1,Ye_a+1)
  
  ed111 = p_ledot_tab(hmq_a,lpres_a,Ye_a)
  ed211 = p_ledot_tab(hmq_a+1,lpres_a,Ye_a)
  ed121 = p_ledot_tab(hmq_a,lpres_a+1,Ye_a)
  ed112 = p_ledot_tab(hmq_a,lpres_a,Ye_a+1)
  ed221 = p_ledot_tab(hmq_a+1,lpres_a+1,Ye_a)
  ed212 = p_ledot_tab(hmq_a+1,lpres_a,Ye_a+1)
  ed122 = p_ledot_tab(hmq_a,lpres_a+1,Ye_a+1)
  ed222 = p_ledot_tab(hmq_a+1,lpres_a+1,Ye_a+1)
  
  yd111 = p_lmYedot_tab(hmq_a,lpres_a,Ye_a)
  yd211 = p_lmYedot_tab(hmq_a+1,lpres_a,Ye_a)
  yd121 = p_lmYedot_tab(hmq_a,lpres_a+1,Ye_a)
  yd112 = p_lmYedot_tab(hmq_a,lpres_a,Ye_a+1)
  yd221 = p_lmYedot_tab(hmq_a+1,lpres_a+1,Ye_a)
  yd212 = p_lmYedot_tab(hmq_a+1,lpres_a,Ye_a+1)
  yd122 = p_lmYedot_tab(hmq_a,lpres_a+1,Ye_a+1)
  yd222 = p_lmYedot_tab(hmq_a+1,lpres_a+1,Ye_a+1)
  
  ab111 = p_lAbar_tab(hmq_a,lpres_a,Ye_a)
  ab211 = p_lAbar_tab(hmq_a+1,lpres_a,Ye_a)
  ab121 = p_lAbar_tab(hmq_a,lpres_a+1,Ye_a)
  ab112 = p_lAbar_tab(hmq_a,lpres_a,Ye_a+1)
  ab221 = p_lAbar_tab(hmq_a+1,lpres_a+1,Ye_a)
  ab212 = p_lAbar_tab(hmq_a+1,lpres_a,Ye_a+1)
  ab122 = p_lAbar_tab(hmq_a,lpres_a+1,Ye_a+1)
  ab222 = p_lAbar_tab(hmq_a+1,lpres_a+1,Ye_a+1)
  
  !! now interpolate

   approxtemp  =    c3                          &
       *( (1.0-c1)*(1.0-c2)    &
       *te112    &
       +      c1*(1.0-c2)    &
       *te212    &
       +      c2*(1.0-c1)    &
       *te122    &
       +      c1*c2          &
       *te222    &
       )                      &
       +(1.0-c3)                    &
       *( (1.0-c1)*(1.0-c2)    &
       *te111    &
       +      c1*(1.0-c2)    &
       *te211    &
       +      c2*(1.0-c1)    &
       *te121    &
       +      c1*c2          &
       *te221    &
       )
  
  approxtemp = 10.e0**approxtemp
  
  qbar_nse  =    c3                          &
       *( (1.0-c1)*(1.0-c2)    &
       *qb112    &
       +      c1*(1.0-c2)    &
       *qb212    &
       +      c2*(1.0-c1)    &
       *qb122    &
       +      c1*c2          &
       *qb222    &
       )                      &
       +(1.0-c3)                    &
       *( (1.0-c1)*(1.0-c2)    &
       *qb111    &
       +      c1*(1.0-c2)    &
       *qb211    &
       +      c2*(1.0-c1)    &
       *qb121    &
       +      c1*c2          &
       *qb221    &
       )
  
  qbar_nse = 10.e0**qbar_nse
  
  edot  =    c3                          &
       *( (1.0-c1)*(1.0-c2)    &
       *ed112    &
       +      c1*(1.0-c2)    &
       *ed212    &
       +      c2*(1.0-c1)    &
       *ed122    &
       +      c1*c2          &
       *ed222    &
       )                      &
       +(1.0-c3)                    &
       *( (1.0-c1)*(1.0-c2)    &
       *ed111    &
       +      c1*(1.0-c2)    &
       *ed211    &
       +      c2*(1.0-c1)    &
       *ed121    &
       +      c1*c2          &
       *ed221    &
       )
  
  edot = 10.e0**edot
  
  Yedot  =    c3                         &
       *( (1.0-c1)*(1.0-c2)    &
       *yd112    &
       +      c1*(1.0-c2)    &
       *yd212    &
       +      c2*(1.0-c1)    &
       *yd122    &
       +      c1*c2          &
       *yd222    &
       )                      &
       +(1.0-c3)                    &
       *( (1.0-c1)*(1.0-c2)    &
       *yd111    &
       +      c1*(1.0-c2)    &
       *yd211    &
       +      c2*(1.0-c1)    &
       *yd121    &
       +      c1*c2          &
       *yd221    &
       )
  
  ! note that Yedot is a negative quantity
  Yedot = - 10.e0**Yedot
  
  abar  =    c3                          &
       *( (1.0-c1)*(1.0-c2)    &
       *ab112    &
       +      c1*(1.0-c2)    &
       *ab212    &
       +      c2*(1.0-c1)    &
       *ab122    &
       +      c1*c2          &
       *ab222    &
       )                      &
       +(1.0-c3)                    &
       *( (1.0-c1)*(1.0-c2)    &
       *ab111    &
       +      c1*(1.0-c2)    &
       *ab211    &
       +      c2*(1.0-c1)    &
       *ab121    &
       +      c1*c2          &
       *ab221    &
       )
  
  abar = 10.e0**abar

  sumyi_nse = 1.0/abar
  
  !write (6,*) 'returning qbar_nse,sumyi_nse,approxtemp,edot,Yedot', &
  !              qbar_nse,sumyi_nse,approxtemp,edot,Yedot
  return

end subroutine
