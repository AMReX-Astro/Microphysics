module react_utils_module

  use microphysics_type_module, only: rt, ZERO
  use variables
  use network
  use eos_type_module
  use eos_module
  use burn_type_module
  use extern_probin_module
  use util_module

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
       dlogT = ZERO
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


  subroutine get_state(state, s_lo, s_hi, ncomp, i, j, k, c, fval) bind(C, name="get_state")
    integer, intent(in) :: s_lo(3), s_hi(3)
    real(rt), intent(inout) :: state(s_lo(1):s_hi(1), s_lo(2):s_hi(2), s_lo(3):s_hi(3), ncomp)
    integer, intent(in), value :: ncomp
    integer, intent(in), value :: i, j, k, c
    real(rt), intent(inout) :: fval
    fval = state(i, j, k, c)
  end subroutine get_state


  subroutine set_state(state, s_lo, s_hi, ncomp, i, j, k, c, fval) bind(C, name="set_state")
    integer, intent(in) :: s_lo(3), s_hi(3)
    real(rt), intent(inout) :: state(s_lo(1):s_hi(1), s_lo(2):s_hi(2), s_lo(3):s_hi(3), ncomp)
    integer, intent(in), value :: ncomp
    integer, intent(in), value :: i, j, k, c
    real(rt), intent(in), value :: fval
    state(i, j, k, c) = fval
  end subroutine set_state


  subroutine get_number_equations(neqs, nspec_not_evolved) bind(C, name="get_number_equations")
    use network, only: nspec, nspec_evolve

    implicit none

    integer, intent(inout) :: neqs, nspec_not_evolved

    neqs = nspec_evolve + 2
    nspec_not_evolved = nspec - nspec_evolve
  end subroutine get_number_equations


  subroutine convert_to_molar(xspec) bind(C, name="sk_convert_to_molar")
    use network, only: aion_inv, nspec_evolve

    implicit none

    real(rt), intent(inout) :: xspec(nspec_evolve)
    integer :: i

    do i = 1, nspec_evolve
       xspec(i) = xspec(i) * aion_inv(i)
    enddo
  end subroutine convert_to_molar


  subroutine get_num_rpar_comps(number_rpar_comps) bind(C, name="sk_get_num_rpar_comps")
    use rpar_indices, only: n_rpar_comps

    implicit none

    integer, intent(inout) :: number_rpar_comps

    number_rpar_comps = n_rpar_comps
  end subroutine get_num_rpar_comps

  subroutine get_nspec_evolve(number_species_evolved) bind(C, name="sk_get_nspec_evolve")
    use network, only: nspec_evolve

    implicit none

    integer, intent(inout) :: number_species_evolved

    number_species_evolved = nspec_evolve
  end subroutine get_nspec_evolve

  subroutine get_csr_jac_rowcols(csr_row_count, csr_col_index) bind(C, name="sk_get_csr_jac_rowcols")

    use cvode_type_module, only: VODE_NEQS
    use network, only: NETWORK_SPARSE_JAC_NNZ, csr_jac_col_index, csr_jac_row_count

    implicit none

    integer, intent(inout) :: csr_row_count(VODE_NEQS+1)
    integer, intent(inout) :: csr_col_index(NETWORK_SPARSE_JAC_NNZ)

    csr_col_index = csr_jac_col_index
    csr_row_count = csr_jac_row_count

  end subroutine get_csr_jac_rowcols

  subroutine get_sparse_jac_nnz(number_nonzero) bind(C, name="sk_get_sparse_jac_nnz")

    use network, only: NETWORK_SPARSE_JAC_NNZ

    implicit none

    integer, intent(inout) :: number_nonzero

    number_nonzero = NETWORK_SPARSE_JAC_NNZ

  end subroutine get_sparse_jac_nnz

  subroutine get_store_jacobian(sjac) bind(C, name="sk_get_store_jacobian")

    use extern_probin_module, only: store_jacobian

    implicit none

    integer, intent(inout) :: sjac

    sjac = store_jacobian

  end subroutine get_store_jacobian

  subroutine get_num_steps_save_jacobian(sjac) bind(C, name="sk_get_num_steps_save_jacobian")

    use extern_probin_module, only: num_steps_save_jacobian

    implicit none

    integer, intent(inout) :: sjac

    sjac = num_steps_save_jacobian

  end subroutine get_num_steps_save_jacobian

end module react_utils_module
