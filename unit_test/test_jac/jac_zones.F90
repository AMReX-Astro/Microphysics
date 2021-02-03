module react_zones_module

  use variables
  use network
  use eos_type_module
  use eos_module
  use burn_type_module
  use amrex_fort_module, only : rt => amrex_real
  use amrex_constants_module
  use extern_probin_module
  use actual_rhs_module
  use numerical_jac_module

  implicit none

contains

  subroutine do_jac(lo, hi, &
                    state, s_lo, s_hi) bind(C, name="do_jac")

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: s_lo(3), s_hi(3)
    real(rt), intent(inout) :: state(s_lo(1):s_hi(1), s_lo(2):s_hi(2), s_lo(3):s_hi(3), p % n_plot_comps)

    type (burn_t)   :: burn_state
    type (eos_t)    :: eos_state
    integer         :: ii, jj, kk, j, i, n
    real(rt)        :: jac_analytic(neqs, neqs), jac_numerical(neqs, neqs)

    !$gpu

    do ii = lo(1), hi(1)
       do jj = lo(2), hi(2)
          do kk = lo(3), hi(3)

             burn_state % rho = state(ii, jj, kk, p % irho)
             burn_state % T = state(ii, jj, kk, p % itemp)
             do j = 1, nspec
                burn_state % xn(j) = state(ii, jj, kk, p % ispec_old + j - 1)
             enddo

             call normalize_abundances_burn(burn_state)

             call burn_to_eos(burn_state, eos_state);
             call eos(eos_input_rt, eos_state);
             call eos_to_burn(eos_state, burn_state);

             ! the integrator doesn't actually care about the initial internal
             ! energy.
             burn_state % e = ZERO

             burn_state % i = ii
             burn_state % j = jj
             burn_state % k = kk

             burn_state % self_heat = .true.

             call actual_jac(burn_state, jac_analytic)

             ! we need to convert this to be in terms of X instead of Y
             do n = 1, nspec
                jac_analytic(n,:) = jac_analytic(n,:) * aion(n)
                jac_analytic(:,n) = jac_analytic(:,n) * aion_inv(n)
             enddo


             call numerical_jac(burn_state, jac_numerical)

             n = 0
             do j = 1, neqs
                do i = 1, neqs
                   state(ii, jj, kk, p % ijac + n) = abs(jac_numerical(i, j) - jac_analytic(i, j))
                   if (jac_analytic(i, j) /= 0.0_rt) then
                      state(ii, jj, kk, p % ijac + n) = state(ii, jj, kk, p % ijac + n) / &
                           abs(jac_analytic(i, j))
                   end if
                   n = n + 1
                end do
             end do

          enddo
       enddo
    enddo

  end subroutine do_jac

end module react_zones_module
