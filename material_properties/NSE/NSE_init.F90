! This function reads in the tables used to interpolate properties of the NSE
! final state.  Two tables sets are read from separate files.
!
! Dean Townsley 2007,8
!

subroutine NSE_init()
  
  use NSE_data, ONLY: nemq, emq_grid, ldens_grid, demq, dldens, nldens, &
       & d_dYe, d_Ye_grid, d_nYe, &
       & d_ltemp_tab, d_lqbar_tab, d_ledot_tab, d_lmYedot_tab, d_lAbar_tab, &
       &  nhmq, hmq_grid, lpres_grid, dhmq, dlpres, nlpres, &
       & p_dYe, p_Ye_grid, p_nYe, p_ltemp_tab, p_lqbar_tab, p_ledot_tab, &
       & p_lmYedot_tab, p_lAbar_tab

  use bl_error_module
  
  implicit none

  character (len=100),save :: prestablename, denstablename
  integer :: i, j, k, istat
  real    :: rtemp1, rtemp2, rtemp3, rtemp4, rtemp5
  real    :: rtemp6, rtemp7, rtemp8, rtemp9


  ! get table names

  prestablename = "nse_pres_hmq_table.txt"
  denstablename = "nse_dens_hmq_table.txt"


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
