module react_zones_module
#ifdef CUDA
  use cudafor
#endif
  use bl_types
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
  end type pfidx_t
  
contains

#ifdef CUDA  
  attributes(global) &
#endif
  subroutine react_zones(state, pfidx)
    real(kind=dp_t) :: state(:,:,:,:)  
    type(pfidx_t), value :: pfidx
    type (burn_t)   :: burn_state_in, burn_state_out
    integer         :: ii, jj, kk, j    
    integer         :: dim_state(4)
    
    dim_state = shape(state)

#ifdef CUDA    
    ii = (blockIdx%x - 1) * blockDim % x + threadIdx % x
    jj = (blockIdx%y - 1) * blockDim % y + threadIdx % y
    kk = (blockIdx%z - 1) * blockDim % z + threadIdx % z    

    if (&
         ii <= dim_state(1) .and. &
         jj <= dim_state(2) .and. &
         kk <= dim_state(3)) then
#else
    do ii = 1, dim_state(1)
    do jj = 1, dim_state(2)
    do kk = 1, dim_state(3)
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
    enddo
    enddo
    enddo
#endif
  end subroutine react_zones
  
end module react_zones_module
