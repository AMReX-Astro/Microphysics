module react_zones_module
#ifdef CUDA
  use cudafor
#endif
  use bl_types
  use bl_space
  use bl_constants_module, only: ZERO
  use network, only: nspec
  use burn_type_module, only: burn_t, normalize_abundances_burn
  use actual_burner_module, only : actual_burner
  use managed_probin_module, only: cu_tmax
  
  implicit none

  type :: pfidx_t
     integer :: itemp
     integer :: irho
     integer :: ispec
     integer :: ispec_old
     integer :: irodot
     integer :: irho_hnuc
     integer :: endrho
     integer :: endT
     integer :: endX
     integer :: ncomps
  end type pfidx_t
  
contains

#ifdef CUDA  
  attributes(global) &
#endif
   subroutine react_zones(state, pfidx, lo, hi)
    integer :: lo(MAX_SPACEDIM), hi(MAX_SPACEDIM)
    type(pfidx_t)   :: pfidx
    real(kind=dp_t) :: state(0:, 0:, 0:, :)
    type (burn_t)   :: burn_state_in, burn_state_out
    integer         :: ii, jj, kk, j    

#ifdef CUDA    
    ii = (blockIdx%x - 1) * blockDim % x + threadIdx % x - 1
    jj = (blockIdx%y - 1) * blockDim % y + threadIdx % y - 1
    kk = (blockIdx%z - 1) * blockDim % z + threadIdx % z - 1

    ! if (&
    !      ii >= lo(1) .and. ii <= hi(1) .and. &
    !      jj >= lo(2) .and. jj <= hi(2) .and. &
    !      kk >= lo(3) .and. kk <= hi(3)) then
    if (ii .eq. 1 .and. jj .eq. 1 .and. kk .eq. 1) then
#else
    ! do ii = lo(1), hi(1)
    ! do jj = lo(2), hi(2)
    ! do kk = lo(3), hi(3)
    ii = 1
    jj = 1
    kk = 1
#endif
       burn_state_in % rho = state(ii, jj, kk, pfidx % irho)
       burn_state_in % T = state(ii, jj, kk, pfidx % itemp)
       do j = 1, nspec
          burn_state_in % xn(j) = state(ii, jj, kk, pfidx % ispec_old + j - 1)
       enddo

       call normalize_abundances_burn(burn_state_in)

       ! the integrator doesn't actually care about the initial internal
       ! energy.
       burn_state_in % e = ZERO

       call actual_burner(burn_state_in, burn_state_out, cu_tmax, ZERO)

       do j = 1, nspec
          state(ii, jj, kk, pfidx % ispec + j - 1) = burn_state_out % xn(j)
       enddo

       do j=1, nspec
          ! an explicit loop is needed here to keep the GPU happy if running on GPU
          state(ii, jj, kk, pfidx % irodot + j - 1) = &
               (burn_state_out % xn(j) - burn_state_in % xn(j)) / cu_tmax
       enddo

       state(ii, jj, kk, pfidx % irho_hnuc) = &
            state(ii, jj, kk, pfidx % irho) * (burn_state_out % e - burn_state_in % e) / cu_tmax
#ifdef CUDA       
    end if
#else
    ! enddo
    ! enddo
    ! enddo
#endif
  end subroutine react_zones
  
end module react_zones_module
