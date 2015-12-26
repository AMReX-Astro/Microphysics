
module NSE_data

  implicit none

  ! pressure table  (independent variables Ye, Pressure,
  !   H-q enthalpy minus nuclear binding energy, both per gram)
  ! we are gridding interpolation in log pressure, linear in the others
  integer,save  :: p_nYe, nlpres, nhmq
  double precision,save :: p_dYe, dlpres, dhmq
  double precision, save, dimension(:), allocatable :: p_Ye_grid, lpres_grid, hmq_grid
  ! all these are stored as log10s and mYedot is log10(-yedot)
  ! the order of indices is hmq, lpres, ye
  ! this matches the order in the datafile
  double precision, save, dimension(:,:,:), allocatable :: p_ltemp_tab, p_lqbar_tab, &
       p_ledot_tab, p_lmYedot_tab, p_lAbar_tab

  ! density table  (independent variables Ye, density,
  !   eint-q internal energy minus nuclear binding energy, both per gram)
  ! we are gridding interpolaton in log density, linear in the others
  integer,save  :: d_nYe, nldens, nemq
  double precision,save :: d_dYe, dldens, demq
  double precision, save, dimension(:), allocatable :: d_Ye_grid, ldens_grid, emq_grid
  ! the order of indices is emq, ldens, ye
  double precision, save, dimension(:,:,:), allocatable :: d_ltemp_tab, d_lqbar_tab, &
       d_ledot_tab, d_lmYedot_tab, d_lAbar_tab

contains

  subroutine NSE_init()

    use bl_error_module

    implicit none

    character (len=100),save :: prestablename, denstablename
    integer :: i, j, k, istat
    double precision    :: rtemp1, rtemp2, rtemp3, rtemp4, rtemp5
    double precision    :: rtemp6, rtemp7, rtemp8, rtemp9


    ! get table names

    prestablename = "nse_pres_hmq_table.txt"
    denstablename = "nse_dens_emq_table.txt"


    !-------------------------------------------------
    ! read table for pressure-based final state lookup
    !-------------------------------------------------

    open (unit=21,file=prestablename,status='OLD',iostat=istat)
    if (istat /= 0) call bl_error("Unable to open nse pressure table")

    read(21,*) p_nYe
    read(21,*) nlpres
    read(21,*) nhmq
    read(21,*) 

    ! space for coordinate grid
    allocate(p_Ye_grid(p_nYe),STAT=istat)
    if (istat /= 0) call bl_error("Cannot allocate p_Ye_grid in NSE_init")
    allocate(lpres_grid(nlpres),STAT=istat)
    if (istat /= 0) call bl_error("Cannot allocate lpres_grid in NSE_init")
    allocate(hmq_grid(nhmq),STAT=istat)
    if (istat /= 0) call bl_error("Cannot allocate hmq_grid in NSE_init")

    ! space for tables
    allocate(p_ltemp_tab(nhmq,nlpres,p_nYe),STAT=istat)
    if (istat /= 0) call bl_error("Cannot allocate p_qbartab in NSE_init")
    allocate(p_lqbar_tab(nhmq,nlpres,p_nYe),STAT=istat)
    if (istat /= 0) call bl_error("Cannot allocate p_qbartab in NSE_init")
    allocate(p_ledot_tab(nhmq,nlpres,p_nYe),STAT=istat)
    if (istat /= 0) call bl_error("Cannot allocate p_edottab in NSE_init")
    allocate(p_lmYedot_tab(nhmq,nlpres,p_nYe),STAT=istat)
    if (istat /= 0) call bl_error("Cannot allocate p_Yedottab in NSE_init")
    allocate(p_lAbar_tab(nhmq,nlpres,p_nYe),STAT=istat)
    if (istat /= 0) call bl_error("Cannot allocate p_Abartab in NSE_init")

    !! read the table, taking logs will be done separately
    do k = 1, p_nYe
       do j = 1, nlpres
          do i = 1, nhmq
             read(21,*) rtemp1, rtemp2, hmq_grid(i), rtemp4, p_ltemp_tab(i,j,k), &
                  rtemp6, p_lqbar_tab(i,j,k), p_lAbar_tab(i,j,k), rtemp9, p_ledot_tab(i,j,k), &
                  p_lmYedot_tab(i,j,k)

             ! avoid posttive Yedot
             if (p_lmYedot_tab(i,j,k).gt.0.0) p_lmYedot_tab(i,j,k) = -1.e-20

             ! empty portions of table are filled with values at the last good
             ! point up the column in hmq (actually done by propagating the value)
             if (p_lqbar_tab(i,j,k) == 0.0) then
                p_ltemp_tab(i,j,k) = p_ltemp_tab(i-1,j,k)
                p_lqbar_tab(i,j,k) = p_lqbar_tab(i-1,j,k)
                p_lAbar_tab(i,j,k) = p_lAbar_tab(i-1,j,k)
                p_ledot_tab(i,j,k) = p_ledot_tab(i-1,j,k)
                p_lmYedot_tab(i,j,k) = p_lmYedot_tab(i-1,j,k)
             endif

          enddo
          lpres_grid(j) = log10(rtemp2)
       enddo
       p_Ye_grid(k) = rtemp1
    enddo
    close(21)

    !! work with logs for temp, qbar, edot, -Yedot, abartables
    ! working with logs of Yedot. Change the sign to make it a positive
    ! quantity, then change it back after interpolating
    p_ltemp_tab = log10(p_ltemp_tab)
    p_lqbar_tab = log10(p_lqbar_tab)
    p_ledot_tab = log10(p_ledot_tab)
    p_lmYedot_tab = log10(-p_lmYedot_tab)
    p_lAbar_tab = log10(p_lAbar_tab)

    !! store the deltas of grid for fast index calculation
    p_dYe = (p_Ye_grid(p_nYe)-p_Ye_grid(1))/(p_nYe-1)
    dlpres = (lpres_grid(nlpres)-lpres_grid(1))/(nlpres-1)
    dhmq = (hmq_grid(nhmq)-hmq_grid(1))/(nhmq-1)



    !----------------------------------------------------------
    ! read table for density-based final state lookup
    !----------------------------------------------------------
    open (unit=21,file=denstablename,status='OLD',iostat=istat)
    if (istat /= 0) call bl_error("Unable to open nse density table")

    read(21,*) d_nYe
    read(21,*) nldens
    read(21,*) nemq
    read(21,*) 

    ! space for coordinate grid
    allocate(d_Ye_grid(d_nYe),STAT=istat)
    if (istat /= 0) call bl_error("Cannot allocate d_Ye_grid in NSE_init")
    allocate(ldens_grid(nldens),STAT=istat)
    if (istat /= 0) call bl_error("Cannot allocate ldens_grid in NSE_init")
    allocate(emq_grid(nemq),STAT=istat)
    if (istat /= 0) call bl_error("Cannot allocate emq_grid in NSE_init")

    ! space for tables
    allocate(d_ltemp_tab(nemq,nldens,d_nYe),STAT=istat)
    if (istat /= 0) call bl_error("Cannot allocate d_qbartab in NSE_init")
    allocate(d_lqbar_tab(nemq,nldens,d_nYe),STAT=istat)
    if (istat /= 0) call bl_error("Cannot allocate d_qbartab in NSE_init")
    allocate(d_ledot_tab(nemq,nldens,d_nYe),STAT=istat)
    if (istat /= 0) call bl_error("Cannot allocate d_edottab in NSE_init")
    allocate(d_lmYedot_tab(nemq,nldens,d_nYe),STAT=istat)
    if (istat /= 0) call bl_error("Cannot allocate d_Yedottab in NSE_init")
    allocate(d_lAbar_tab(nemq,nldens,d_nYe),STAT=istat)
    if (istat /= 0) call bl_error("Cannot allocate d_Abartab in NSE_init")

    !! read the table data
    do k = 1, d_nYe
       do j = 1, nldens
          do i = 1, nemq
             read(21,*) rtemp1, rtemp2, emq_grid(i), rtemp4, d_ltemp_tab(i,j,k), &
                  rtemp6, d_lqbar_tab(i,j,k), d_lAbar_tab(i,j,k), rtemp9, d_ledot_tab(i,j,k), &
                  d_lmYedot_tab(i,j,k)

             ! avoid posttive Yedot
             if (d_lmYedot_tab(i,j,k).gt.0.0) d_lmYedot_tab(i,j,k) = -1.e-20

             ! empty portions of table are filled with values at the last good
             ! point up the column in hmq (actually done by propagating the value)
             if (d_lqbar_tab(i,j,k) == 0.0) then
                d_ltemp_tab(i,j,k) = d_ltemp_tab(i-1,j,k)
                d_lqbar_tab(i,j,k) = d_lqbar_tab(i-1,j,k)
                d_lAbar_tab(i,j,k) = d_lAbar_tab(i-1,j,k)
                d_ledot_tab(i,j,k) = d_ledot_tab(i-1,j,k)
                d_lmYedot_tab(i,j,k) = d_lmYedot_tab(i-1,j,k)
             endif

          enddo
          ldens_grid(j) = log10(rtemp2)
       enddo
       d_Ye_grid(k) = rtemp1
    enddo
    close(21)


    ! work with logs for temp, qbar, edot, -Yedot, abartables
    ! working with logs of Yedot. Change the sign to make it a positive
    ! quantity, then change it back after interpolating

    d_ltemp_tab   = log10(d_ltemp_tab)
    d_lqbar_tab   = log10(d_lqbar_tab)
    d_ledot_tab   = log10(d_ledot_tab)
    d_lmYedot_tab = log10(-d_lmYedot_tab)
    d_lAbar_tab   = log10(d_lAbar_tab)

    !! store the deltas of grid for fast index calculation
    d_dYe  = (d_Ye_grid(d_nYe)-d_Ye_grid(1))/(d_nYe-1)
    dldens = (ldens_grid(nldens)-ldens_grid(1))/(nldens-1)
    demq   = (emq_grid(nemq)-emq_grid(1))/(nemq-1)


  end subroutine NSE_init


  
  subroutine NSE_finalAtDens(qbar_nse,sumyi_nse,approxtemp,edot,Yedot, Ye, dens, emq)

    implicit none

    double precision, intent(IN) :: Ye, dens, emq
    double precision, intent(OUT) :: qbar_nse,sumyi_nse,approxtemp,edot,Yedot

    integer :: emq_a, ldens_a, Ye_a

    double precision :: ldens, emq_v
    double precision :: te111, te211, te121, te112, te221, te212, te122, te222
    double precision :: qb111, qb211, qb121, qb112, qb221, qb212, qb122, qb222
    double precision :: ed111, ed211, ed121, ed112, ed221, ed212, ed122, ed222
    double precision :: yd111, yd211, yd121, yd112, yd221, yd212, yd122, yd222
    double precision :: ab111, ab211, ab121, ab112, ab221, ab212, ab122, ab222

    double precision :: c1, c2, c3, abar


    !! find the location in the table grid
    !! remember that we use the log of density

    emq_v = emq
    ldens = log10(dens) 


    ! we want to truncate at the "bottom" edges (in density and emq)
    !instead of extrapolate.  Should be pure nickel down here.
    ! emq is stored in drecrasing order

    if ( emq_v < emq_grid(nemq)   )   emq_v = emq_grid(nemq)
    if ( ldens < ldens_grid(1) ) ldens = ldens_grid(1)

    ! find v_a  such that v_grid(v_a) <= v < v_grid(v_a+1)
    emq_a   = floor((emq_v - emq_grid(1))/demq)+1
    ldens_a = floor((ldens - ldens_grid(1))/dldens)+1
    Ye_a   = floor((Ye - d_Ye_grid(1))/d_dYe)+1

    ! at upper edges we will extrapolate
    if (emq_a < 1) emq_a = 1
    if (emq_a >= nemq) emq_a = nemq - 1
    ! lower limit handled above
    if (ldens_a >= nldens) ldens_a = nldens - 1

    ! extrapolate both ways in Ye
    if (Ye_a >= d_nYe) Ye_a = d_nYe - 1
    if (Ye_a < 1) Ye_a = 1

    !! calculate the coefficients

    c1 = ( emq_v - emq_grid(emq_a) )/                &
         ( emq_grid(emq_a+1) - emq_grid(emq_a) )
    c2 = ( ldens - ldens_grid(ldens_a)) /              &
         ( ldens_grid(ldens_a+1) - ldens_grid(ldens_a) )
    c3 = ( Ye - d_Ye_grid(Ye_a) )/                  &
         ( d_Ye_grid(Ye_a+1) - d_Ye_grid(Ye_a) )

    !! build the local cubes
    te111 = d_ltemp_tab(emq_a,ldens_a,Ye_a)
    te211 = d_ltemp_tab(emq_a+1,ldens_a,Ye_a)
    te121 = d_ltemp_tab(emq_a,ldens_a+1,Ye_a)
    te112 = d_ltemp_tab(emq_a,ldens_a,Ye_a+1)
    te221 = d_ltemp_tab(emq_a+1,ldens_a+1,Ye_a)
    te212 = d_ltemp_tab(emq_a+1,ldens_a,Ye_a+1)
    te122 = d_ltemp_tab(emq_a,ldens_a+1,Ye_a+1)
    te222 = d_ltemp_tab(emq_a+1,ldens_a+1,Ye_a+1)

    qb111 = d_lqbar_tab(emq_a,ldens_a,Ye_a)
    qb211 = d_lqbar_tab(emq_a+1,ldens_a,Ye_a)
    qb121 = d_lqbar_tab(emq_a,ldens_a+1,Ye_a)
    qb112 = d_lqbar_tab(emq_a,ldens_a,Ye_a+1)
    qb221 = d_lqbar_tab(emq_a+1,ldens_a+1,Ye_a)
    qb212 = d_lqbar_tab(emq_a+1,ldens_a,Ye_a+1)
    qb122 = d_lqbar_tab(emq_a,ldens_a+1,Ye_a+1)
    qb222 = d_lqbar_tab(emq_a+1,ldens_a+1,Ye_a+1)

    ed111 = d_ledot_tab(emq_a,ldens_a,Ye_a)
    ed211 = d_ledot_tab(emq_a+1,ldens_a,Ye_a)
    ed121 = d_ledot_tab(emq_a,ldens_a+1,Ye_a)
    ed112 = d_ledot_tab(emq_a,ldens_a,Ye_a+1)
    ed221 = d_ledot_tab(emq_a+1,ldens_a+1,Ye_a)
    ed212 = d_ledot_tab(emq_a+1,ldens_a,Ye_a+1)
    ed122 = d_ledot_tab(emq_a,ldens_a+1,Ye_a+1)
    ed222 = d_ledot_tab(emq_a+1,ldens_a+1,Ye_a+1)

    yd111 = d_lmYedot_tab(emq_a,ldens_a,Ye_a)
    yd211 = d_lmYedot_tab(emq_a+1,ldens_a,Ye_a)
    yd121 = d_lmYedot_tab(emq_a,ldens_a+1,Ye_a)
    yd112 = d_lmYedot_tab(emq_a,ldens_a,Ye_a+1)
    yd221 = d_lmYedot_tab(emq_a+1,ldens_a+1,Ye_a)
    yd212 = d_lmYedot_tab(emq_a+1,ldens_a,Ye_a+1)
    yd122 = d_lmYedot_tab(emq_a,ldens_a+1,Ye_a+1)
    yd222 = d_lmYedot_tab(emq_a+1,ldens_a+1,Ye_a+1)

    ab111 = d_lAbar_tab(emq_a,ldens_a,Ye_a)
    ab211 = d_lAbar_tab(emq_a+1,ldens_a,Ye_a)
    ab121 = d_lAbar_tab(emq_a,ldens_a+1,Ye_a)
    ab112 = d_lAbar_tab(emq_a,ldens_a,Ye_a+1)
    ab221 = d_lAbar_tab(emq_a+1,ldens_a+1,Ye_a)
    ab212 = d_lAbar_tab(emq_a+1,ldens_a,Ye_a+1)
    ab122 = d_lAbar_tab(emq_a,ldens_a+1,Ye_a+1)
    ab222 = d_lAbar_tab(emq_a+1,ldens_a+1,Ye_a+1)


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
  end subroutine NSE_finalAtDens



  subroutine NSE_finalAtPres(qbar_nse,sumyi_nse,approxtemp,edot,Yedot, Ye, pres, hmq)

    implicit none

    double precision, intent(IN)    :: Ye, pres, hmq
    double precision, intent(OUT)   :: qbar_nse,sumyi_nse,approxtemp,edot,Yedot

    integer :: hmq_a, lpres_a, Ye_a

    double precision :: lpres, hmq_v
    double precision :: te111, te211, te121, te112, te221, te212, te122, te222
    double precision :: qb111, qb211, qb121, qb112, qb221, qb212, qb122, qb222
    double precision :: ed111, ed211, ed121, ed112, ed221, ed212, ed122, ed222
    double precision :: yd111, yd211, yd121, yd112, yd221, yd212, yd122, yd222
    double precision :: ab111, ab211, ab121, ab112, ab221, ab212, ab122, ab222

    double precision :: c1, c2, c3, abar


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

  end subroutine NSE_finalAtPres
  


  subroutine NSE_finalize()

    implicit none

    ! deallocate pressure-based table
    deallocate(p_Ye_grid)
    deallocate(lpres_grid)
    deallocate(hmq_grid)
    deallocate(p_ltemp_tab)
    deallocate(p_lqbar_tab)
    deallocate(p_ledot_tab)
    deallocate(p_lmYedot_tab)
    deallocate(p_lAbar_tab)

    ! deallocate density-based table
    deallocate(d_Ye_grid)
    deallocate(ldens_grid)
    deallocate(emq_grid)
    deallocate(d_ltemp_tab)
    deallocate(d_lqbar_tab)
    deallocate(d_ledot_tab)
    deallocate(d_lmYedot_tab)
    deallocate(d_lAbar_tab)

  end subroutine NSE_finalize
  
end module NSE_data
