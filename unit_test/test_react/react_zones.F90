module react_zones_module

  use variables
  use network
  use eos_type_module
  use eos_module
  use burn_type_module
  use amrex_fort_module, only : rt => amrex_real
  use amrex_constants_module
  use extern_probin_module
  use actual_burner_module

  implicit none

contains

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

#if NAUX_NET > 0
             do j = 1, naux
                burn_state_in % aux(j) = state(ii, jj, kk, p % iaux_old + j - 1)
             end do
#endif

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

#if NAUX_NET > 0
             do j = 1, naux
                state(ii, jj, kk, p % iaux + j - 1) = burn_state_out % aux(j)
             end do
#endif


             state(ii, jj, kk, p % irho_hnuc) = &
                  state(ii, jj, kk, p % irho) * (burn_state_out % e - burn_state_in % e) / tmax

             n_rhs(ii, jj, kk, 1) = burn_state_out % n_rhs

          enddo
       enddo
    enddo

  end subroutine do_react

end module react_zones_module
