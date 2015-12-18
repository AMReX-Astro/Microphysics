! data structures for storing tables of NSE final states and metadata about them
!
! Dean Townsley 2007,8
!

module NSE_data

  implicit none

  ! pressure table  (independent variables Ye, Pressure,
  !   H-q enthalpy minus nuclear binding energy, both per gram)
  ! we are gridding interpolation in log pressure, linear in the others
  integer,save  :: p_nYe, nlpres, nhmq
  real,save :: p_dYe, dlpres, dhmq
  real, save, dimension(:), allocatable :: p_Ye_grid, lpres_grid, hmq_grid
  ! all these are stored as log10s and mYedot is log10(-yedot)
  ! the order of indices is hmq, lpres, ye
  ! this matches the order in the datafile
  real, save, dimension(:,:,:), allocatable :: p_ltemp_tab, p_lqbar_tab, &
       p_ledot_tab, p_lmYedot_tab, p_lAbar_tab

  ! density table  (independent variables Ye, density,
  !   eint-q internal energy minus nuclear binding energy, both per gram)
  ! we are gridding interpolaton in log density, linear in the others
  integer,save  :: d_nYe, nldens, nemq
  real,save :: d_dYe, dldens, demq
  real, save, dimension(:), allocatable :: d_Ye_grid, ldens_grid, emq_grid
  ! the order of indices is emq, ldens, ye
  real, save, dimension(:,:,:), allocatable :: d_ltemp_tab, d_lqbar_tab, &
       d_ledot_tab, d_lmYedot_tab, d_lAbar_tab

end module NSE_data
