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
  subroutine react_zones(state, pfidx, coffset, cend, sLength)
    implicit none
  
    integer, value, intent(in)  :: coffset, cend, sLength
    type(pfidx_t),   intent(in) :: pfidx
    real(kind=dp_t) &
#ifdef CUDA
         , device &
#endif     
         , intent(inout) :: state(:,:)
    type (burn_t)   :: burn_state_in, burn_state_out
    integer         :: ii, j

#ifdef CUDA    
    ii = (blockIdx % x - 1) * blockDim % x + threadIdx % x

    if (ii > cend) return
#else
    !$OMP PARALLEL DO PRIVATE(ii,j) &
    !$OMP PRIVATE(burn_state_in, burn_state_out) &
    !$OMP SCHEDULE(DYNAMIC,1)
    do ii = lo, hi
#endif
       ii = ii + coffset
       burn_state_in % rho = state(pfidx % irho, ii)
       burn_state_in % T = state(pfidx % itemp, ii)
       do j = 1, nspec
          burn_state_in % xn(j) = state(pfidx % ispec_old + j - 1, ii)
       enddo

       call normalize_abundances_burn(burn_state_in)

       ! the integrator doesn't actually care about the initial internal
       ! energy.
       burn_state_in % e = ZERO

       call actual_burner(burn_state_in, burn_state_out, cu_tmax, ZERO)

       do j = 1, nspec
          state(pfidx % ispec + j - 1, ii) = burn_state_out % xn(j)
       enddo

       do j=1, nspec
          ! an explicit loop is needed here to keep the GPU happy if running on GPU
          state(pfidx % irodot + j - 1, ii) = &
               (burn_state_out % xn(j) - burn_state_in % xn(j)) / cu_tmax
       enddo

       state(pfidx % irho_hnuc, ii) = &
            state(pfidx % irho, ii) * (burn_state_out % e - burn_state_in % e) / cu_tmax
#ifndef CUDA       
    enddo
    !$OMP END PARALLEL DO
#endif
  end subroutine react_zones
  
end module react_zones_module
