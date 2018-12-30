module react_zones_module

  use variables
  use network
  use eos_type_module
  use eos_module
  use burn_type_module
  use amrex_fort_module, only : rt => amrex_real
  use amrex_constants_module
  use extern_probin_module
  use util_module
  use actual_burner_module

  implicit none

contains

  subroutine init_state(lo, hi, &
                        state, s_lo, s_hi, npts) bind(C, name="init_state")

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: s_lo(3), s_hi(3)
    real(rt), intent(inout) :: state(s_lo(1):s_hi(1), s_lo(2):s_hi(2), s_lo(3):s_hi(3), p % n_plot_comps)
    integer, intent(in) :: npts

    real(rt) :: dlogrho, dlogT, dmetal
    real(rt) :: temp_zone, dens_zone
    real(rt), allocatable :: xn_zone(:,:)

    integer :: ii, jj, kk
    real(rt) :: sum_X

    if (npts > 1) then
       dlogrho = (log10(dens_max) - log10(dens_min))/(npts - 1)
       dlogT   = (log10(temp_max) - log10(temp_min))/(npts - 1)
    else
       dlogrho = ZERO
       dlogT   = ZERO
    endif

    allocate(xn_zone(nspec, 0:npts-1))   ! this assumes that lo(3) = 0

    call get_xn(xn_zone)

    ! normalize -- just in case
    do kk = lo(3), hi(3)
       sum_X = sum(xn_zone(:, kk))
       xn_zone(:, kk) = xn_zone(:, kk)/sum_X
    enddo

    do kk = lo(3), hi(3)
       do jj = lo(2), hi(2)
          do ii = lo(1), hi(1)

             state(ii, jj, kk, p % itemp) = 10.0_rt**(log10(temp_min) + dble(jj)*dlogT)
             state(ii, jj, kk, p % irho)  = 10.0_rt**(log10(dens_min) + dble(ii)*dlogrho)
             state(ii, jj, kk, p % ispec_old:p % ispec_old+nspec-1) = max(xn_zone(:, kk), 1.e-10_rt)

          enddo
       enddo
    enddo

  end subroutine init_state


  subroutine print_nrhs(lo, hi, state, s_lo, s_hi) bind(C, name="print_nrhs")

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: s_lo(3), s_hi(3)
    integer, intent(inout) :: state(s_lo(1):s_hi(1), s_lo(2):s_hi(2), s_lo(3):s_hi(3), 1)

    integer :: ii, jj, kk

    do ii = lo(1), hi(1)
       do jj = lo(2), hi(2)
          do kk = lo(3), hi(3)
             print *, ' nrhs for ', ii, jj, kk, ' = ', state(ii,jj,kk,1)
          enddo
       enddo
    enddo

  end subroutine print_nrhs


  subroutine do_react(lo, hi, &
                      state, s_lo, s_hi, &
                      n_rhs, n_rhs_lo, n_rhs_hi) bind(C, name="do_react")

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: s_lo(3), s_hi(3)
    integer, intent(in) :: n_rhs_lo(3), n_rhs_hi(3)
    real(rt), intent(inout) :: state(s_lo(1):s_hi(1), s_lo(2):s_hi(2), s_lo(3):s_hi(3), p % n_plot_comps)
    integer,  intent(inout) :: n_rhs(n_rhs_lo(1):n_rhs_hi(1), n_rhs_lo(2):n_rhs_hi(2), n_rhs_lo(3):n_rhs_hi(3), 1)

    type (burn_t)   :: burn_state_in, burn_state_out
    integer         :: ii, jj, kk, j

    !$gpu

    do ii = lo(1), hi(1)
       do jj = lo(2), hi(2)
          do kk = lo(3), hi(3)

             burn_state_in % rho = state(ii, jj, kk, p % irho)
             burn_state_in % T = state(ii, jj, kk, p % itemp)
             do j = 1, nspec
                burn_state_in % xn(j) = state(ii, jj, kk, p % ispec_old + j - 1)
             enddo

             call normalize_abundances_burn(burn_state_in)

             ! the integrator doesn't actually care about the initial internal
             ! energy.
             burn_state_in % e = ZERO

             burn_state_in % i = ii
             burn_state_in % j = jj
             burn_state_in % k = kk

             call actual_burner(burn_state_in, burn_state_out, tmax, ZERO)

             do j = 1, nspec
                state(ii, jj, kk, p % ispec + j - 1) = burn_state_out % xn(j)
             enddo

             do j = 1, nspec
                ! an explicit loop is needed here to keep the GPU happy if running on GPU
                state(ii, jj, kk, p % irodot + j - 1) = &
                     (burn_state_out % xn(j) - burn_state_in % xn(j)) / tmax
             enddo

             state(ii, jj, kk, p % irho_hnuc) = &
                  state(ii, jj, kk, p % irho) * (burn_state_out % e - burn_state_in % e) / tmax

             n_rhs(ii, jj, kk, 1) = burn_state_out % n_rhs

          enddo
       enddo
    enddo

  end subroutine do_react

end module react_zones_module
