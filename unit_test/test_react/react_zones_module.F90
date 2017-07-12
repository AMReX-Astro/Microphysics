module react_zones_module
#ifdef CUDA
  use cudafor
#endif
  use bl_types
  use bl_space
  use bl_constants_module, only: ZERO
  use probin_module, only: tmax
  use network, only: nspec
  use burn_type_module, only: burn_t, normalize_abundances_burn
  use actual_burner_module, only : actual_burner
  
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

  type(pfidx_t), &
#ifdef CUDA
       managed, &
#endif
       allocatable :: pfidx
  
contains

#ifdef CUDA  
  attributes(global) subroutine react_zones(state, coffset, cend, sPitch, sLength)
#else
  subroutine react_zones(state, lo, hi)
#endif
      
    implicit none

#ifdef CUDA    
    integer, value, intent(in) :: coffset, cend, sPitch, sLength
#else
    integer, intent(in) :: lo, hi
#endif
    
    real(kind=dp_t), &
#ifdef CUDA
         device, intent(inout) :: state(1:sPitch, 1:sLength)
#else
                 intent(inout) :: state(1:pfidx % ncomps, lo:hi)
#endif
         
    type (burn_t)   :: burn_state_in, burn_state_out
    integer         :: ii, j

#ifdef CUDA 
    ii = (blockIdx % x - 1) * blockDim % x + threadIdx % x

    if (ii > cend) then
       return
    else
       ii = ii + coffset
    end if
#else
    !$OMP PARALLEL DO PRIVATE(ii,j) &
    !$OMP PRIVATE(burn_state_in, burn_state_out) &
    !$OMP SCHEDULE(DYNAMIC,1)
    do ii = lo, hi
#endif
       burn_state_in % rho = state(pfidx % irho, ii)
       burn_state_in % T = state(pfidx % itemp, ii)
       do j = 1, nspec
          burn_state_in % xn(j) = state(pfidx % ispec_old + j - 1, ii)
       enddo

       call normalize_abundances_burn(burn_state_in)

       ! the integrator doesn't actually care about the initial internal
       ! energy.
       burn_state_in % e = ZERO

       call actual_burner(burn_state_in, burn_state_out, tmax, ZERO)

       do j = 1, nspec
          state(pfidx % ispec + j - 1, ii) = burn_state_out % xn(j)
       enddo

       do j=1, nspec
          ! an explicit loop is needed here to keep the GPU happy if running on GPU
          state(pfidx % irodot + j - 1, ii) = &
               (burn_state_out % xn(j) - burn_state_in % xn(j)) / tmax
       enddo

       state(pfidx % irho_hnuc, ii) = &
            state(pfidx % irho, ii) * (burn_state_out % e - burn_state_in % e) / tmax
#ifndef CUDA       
    enddo
    !$OMP END PARALLEL DO
#endif
  end subroutine react_zones
  
end module react_zones_module
