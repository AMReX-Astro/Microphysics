! free stuff
!
! Dean Townsley 2008
!

subroutine NSE_finalize()

  use NSE_data

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

end subroutine
