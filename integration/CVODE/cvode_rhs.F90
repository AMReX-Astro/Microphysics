module cvode_rhs_module

  implicit none
  
contains
  
  ! The f_rhs routine provides the right-hand-side for the DVODE solver.
  ! This is a generic interface that calls the specific RHS routine in the
  ! network you're actually using.

  subroutine sk_f_rhs(time, y, ydot, rpar) bind(C, name="sk_f_rhs")

    !$acc routine seq
    
    use actual_network, only: aion, nspec_evolve
    use amrex_fort_module, only: rt => amrex_real
    use burn_type_module, only: burn_t, net_ienuc, net_itemp
    use amrex_constants_module, only: ZERO, ONE
    use actual_rhs_module, only: actual_rhs
    use extern_probin_module, only: renormalize_abundances, &
         integrate_temperature, integrate_energy
    use cvode_type_module, only: sk_clean_state, sk_renormalize_species, sk_update_thermodynamics, burn_to_vode, vode_to_burn, VODE_NEQS
    use cvode_rpar_indices, only: n_rpar_comps, irp_y_init, irp_t_sound

    implicit none

    real(rt), intent(INOUT) :: time, y(VODE_NEQS)
    real(rt), intent(INOUT) :: rpar(n_rpar_comps)
    real(rt), intent(INOUT) :: ydot(VODE_NEQS)

    type (burn_t) :: burn_state

    real(rt) :: limit_factor, t_sound, t_enuc

    !$gpu

    ! We are integrating a system of
    !
    ! y(1:nspec_evolve) = dX/dt
    ! y(net_itemp) = dT/dt
    ! y(net_ienuc) = denuc/dt

    ydot = ZERO

    ! Fix the state as necessary.

    call sk_clean_state(y, rpar)

    ! Renormalize the abundances as necessary.

    if (renormalize_abundances) then
       call sk_renormalize_species(y, rpar)
    endif

    ! Update the thermodynamics as necessary.

    call sk_update_thermodynamics(y, rpar)

    ! Call the specific network routine to get the RHS.

    call vode_to_burn(y, rpar, burn_state)

    burn_state % time = time
    call actual_rhs(burn_state)

    ! We integrate X, not Y
    burn_state % ydot(1:nspec_evolve) = &
         burn_state % ydot(1:nspec_evolve) * aion(1:nspec_evolve)

    ! Allow temperature and energy integration to be disabled.
    if (.not. integrate_temperature) then
       burn_state % ydot(net_itemp) = ZERO
    endif

    if (.not. integrate_energy) then
       burn_state % ydot(net_ienuc) = ZERO
    endif

    call burn_to_vode(burn_state, y, rpar, ydot = ydot)

  end subroutine sk_f_rhs



  ! Analytical Jacobian
  subroutine sk_analytic_jac(time, y, jac_mat, rpar) bind(C, name="sk_analytic_jac")

    !$acc routine seq
    
    use network, only: aion, aion_inv, nspec_evolve, NETWORK_SPARSE_JAC_NNZ
    use amrex_constants_module, only: ZERO
    use actual_rhs_module, only: actual_jac
    use burn_type_module, only: burn_t, net_ienuc, net_itemp
    use jacobian_sparsity_module, only: get_jac_entry, set_jac_entry, scale_jac_entry
    use cvode_type_module, only: vode_to_burn, burn_to_vode, VODE_NEQS
    use cvode_rpar_indices, only: n_rpar_comps, irp_y_init, irp_t_sound
    use amrex_fort_module, only: rt => amrex_real
    use extern_probin_module, only: integrate_temperature, integrate_energy

    implicit none

    real(rt), intent(INOUT) :: y(VODE_NEQS), rpar(n_rpar_comps), time

#ifdef REACT_SPARSE_JACOBIAN
    real(rt), intent(  OUT) :: jac_mat(NETWORK_SPARSE_JAC_NNZ)
#else
    real(rt), intent(  OUT) :: jac_mat(VODE_NEQS,VODE_NEQS)
#endif

    type (burn_t) :: state
    real(rt) :: limit_factor, t_sound, t_enuc, scratch
    integer :: n, i

    !$gpu

    ! Call the specific network routine to get the Jacobian.

    call vode_to_burn(y, rpar, state)
    state % time = time
    call actual_jac(state)

    ! We integrate X, not Y
    do n = 1, nspec_evolve
       do i = 1, VODE_NEQS
          call scale_jac_entry(state, n, i, aion(n))
          call scale_jac_entry(state, i, n, aion_inv(n))
       enddo
    enddo

    ! Allow temperature and energy integration to be disabled.
    if (.not. integrate_temperature) then
       scratch = ZERO
       do i = 1, VODE_NEQS
          call set_jac_entry(state, net_itemp, i, scratch)
       enddo
    endif

    if (.not. integrate_energy) then
       scratch = ZERO
       do i = 1, VODE_NEQS
          call set_jac_entry(state, net_ienuc, i, scratch)
       enddo
    endif

    call burn_to_vode(state, y, rpar, jac = jac_mat)

  end subroutine sk_analytic_jac


  subroutine sk_jac_times_vec(jac_mat, vec, jtv) bind(C, name="sk_jac_times_vec")

#ifdef REACT_SPARSE_JACOBIAN
    use network, only: NETWORK_SPARSE_JAC_NNZ, csr_jac_col_index, csr_jac_row_count
#endif

    use burn_type_module
    use cvode_type_module, only: VODE_NEQS
    use amrex_fort_module, only: rt => amrex_real
    use amrex_constants_module, only: ZERO    

    implicit none

#ifdef REACT_SPARSE_JACOBIAN
    real(rt), intent(IN   ) :: jac_mat(NETWORK_SPARSE_JAC_NNZ)
#else
    real(rt), intent(IN   ) :: jac_mat(VODE_NEQS,VODE_NEQS)
#endif
    real(rt), intent(IN   ) :: vec(VODE_NEQS)
    real(rt), intent(INOUT) :: jtv(VODE_NEQS)

    real(rt) :: scratch
    integer :: i, j

#ifdef REACT_SPARSE_JACOBIAN
    integer :: num_in_row, loc, col
#endif

    !$gpu

#ifdef REACT_SPARSE_JACOBIAN
    loc = 1
    do i = 1, VODE_NEQS
       jtv(i) = ZERO
       num_in_row = csr_jac_row_count(i+1) - csr_jac_row_count(i)
       do j = 1, num_in_row
          col = csr_jac_col_index(loc)
          jtv(i) = jtv(i) + jac_mat(loc) * vec(col)
          loc = loc + 1
       enddo
    enddo
#else
    do i = 1, VODE_NEQS
       jtv(i) = ZERO
       do j = 1, VODE_NEQS
          jtv(i) = jtv(i) + jac_mat(i,j) * vec(j)
       enddo
    enddo
#endif

  end subroutine sk_jac_times_vec

  
  subroutine sk_fill_csr_jac(jac_mat, csr_jac) bind(C, name="sk_fill_csr_jac")
    
    use cvode_type_module, only: VODE_NEQS
    use network, only: NETWORK_SPARSE_JAC_NNZ, csr_jac_col_index, csr_jac_row_count
    use amrex_fort_module, only: rt => amrex_real
    use amrex_constants_module, only: ZERO    

    implicit none

    real(rt), intent(IN   ) :: jac_mat(VODE_NEQS,VODE_NEQS)
    real(rt), intent(INOUT) :: csr_jac(NETWORK_SPARSE_JAC_NNZ)

    integer :: irow, icol, ninrow, j, loc

    !$gpu

    loc = 1
    do irow = 1, VODE_NEQS
       ninrow = csr_jac_row_count(irow+1) - csr_jac_row_count(irow)
       do j = 1, ninrow
          icol = csr_jac_col_index(loc)
          csr_jac(loc) = jac_mat(irow, icol)
          loc = loc + 1
       enddo
    enddo

  end subroutine sk_fill_csr_jac
  
  
  subroutine sk_full_jac(y, jac_mat, rpar, neq_total, ncells, neq_per_cell, nrpar_per_cell) bind(C, name="sk_full_jac")

    use amrex_fort_module, only: rt => amrex_real
    use amrex_constants_module, only: ZERO
    
    implicit none

    integer  :: neq_total, ncells, neq_per_cell, nrpar_per_cell
    real(rt) :: y(neq_total), jac_mat(neq_total, neq_total), rpar(ncells*nrpar_per_cell)

    integer  :: i, offset_y, offset_rpar, offset_y_end, offset_rpar_end
    real(rt) :: t = ZERO

    jac_mat(1:neq_total, 1:neq_total) = ZERO

    do i = 1, ncells
       offset_y = (i-1)*neq_per_cell + 1
       offset_y_end = i*neq_per_cell
       offset_rpar = (i-1)*nrpar_per_cell + 1
       offset_rpar_end = i*nrpar_per_cell       
       call sk_analytic_jac(t, y(offset_y:offset_y_end), &
                            jac_mat(offset_y:offset_y_end, offset_y:offset_y_end), &
                            rpar(offset_rpar:offset_rpar_end))
    enddo

  end subroutine sk_full_jac

end module cvode_rhs_module
