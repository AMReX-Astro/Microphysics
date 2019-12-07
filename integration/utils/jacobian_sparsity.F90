module jacobian_sparsity_module

#ifdef REACT_SPARSE_JACOBIAN
  use actual_network, only: NETWORK_SPARSE_JAC_NNZ
#endif

  use burn_type_module
  use amrex_fort_module, only : rt => amrex_real

  use amrex_fort_module, only : rt => amrex_real
  implicit none

contains

#ifdef REACT_SPARSE_JACOBIAN
  subroutine lookup_csr_jac_loc(row, col, csr_loc)

    !$acc routine seq

    use actual_network, only: csr_jac_col_index, csr_jac_row_count

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer, intent(in) :: row, col
    integer, intent(out) :: csr_loc

    integer :: num_in_row, row_start_loc, row_end_loc, i

    !$gpu

    ! Looks up the index into a CSR-formatted Jacobian
    ! matrix given row and col indices into the
    ! equivalent dense matrix.
    !
    ! Assumes the base in first element of CSR row count array is 1

    num_in_row = csr_jac_row_count(row+1) - csr_jac_row_count(row)
    row_start_loc = csr_jac_row_count(row)
    row_end_loc   = row_start_loc + num_in_row - 1

    csr_loc = -1
    do i = row_start_loc, row_end_loc
       if (csr_jac_col_index(i) == col) then
          csr_loc = i
          exit
       endif
    enddo
  end subroutine lookup_csr_jac_loc


  subroutine set_csr_jac_entry(csr_jac, row, col, val)

    !$acc routine seq

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    real(rt), intent(inout) :: csr_jac(NETWORK_SPARSE_JAC_NNZ)
    integer, intent(in) :: row, col
    real(rt), intent(in) :: val

    integer :: csr_loc

    !$gpu

    ! Get index into the CSR Jacobian
    call lookup_csr_jac_loc(row, col, csr_loc)

    ! Set value in CSR Jacobian only if row, col entry exists
    if (csr_loc /= -1) then
       csr_jac(csr_loc) = val
    else
#ifdef SPARSE_STOP_ON_OOB
       stop
#else
       return
#endif
    endif

  end subroutine set_csr_jac_entry


  subroutine scale_csr_jac_entry(csr_jac, row, col, val)

    !$acc routine seq

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    real(rt), intent(inout) :: csr_jac(NETWORK_SPARSE_JAC_NNZ)
    integer, intent(in) :: row, col
    real(rt), intent(in) :: val

    integer :: csr_loc

    !$gpu

    ! Get index into the CSR Jacobian
    call lookup_csr_jac_loc(row, col, csr_loc)

    ! Scale value in CSR Jacobian only if row, col entry exists
    if (csr_loc /= -1) then
       csr_jac(csr_loc) = csr_jac(csr_loc) * val
    else
#ifdef SPARSE_STOP_ON_OOB
       stop
#else
       return
#endif
    endif

  end subroutine scale_csr_jac_entry


  subroutine get_csr_jac_entry(csr_jac, row, col, val)

    !$acc routine seq

    use amrex_constants_module, only: ZERO

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    real(rt), intent(in) :: csr_jac(NETWORK_SPARSE_JAC_NNZ)
    integer, intent(in) :: row, col
    real(rt), intent(out) :: val

    integer :: csr_loc

    !$gpu

    ! Get index into the CSR Jacobian
    call lookup_csr_jac_loc(row, col, csr_loc)

    ! Get value from CSR Jacobian only if row, col entry exists
    ! otherwise return ZERO.
    if (csr_loc /= -1) then
       val = csr_jac(csr_loc)
    else
       val = ZERO
    endif

  end subroutine get_csr_jac_entry
#endif


  subroutine set_jac_entry(state, row, col, val)

    !$acc routine seq

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    type (burn_t), intent(inout) :: state
    integer, intent(in) :: row, col
    real(rt), intent(in) :: val

    !$gpu

#ifdef REACT_SPARSE_JACOBIAN
    call set_csr_jac_entry(state % sparse_jac, row, col, val)
#else
    state % jac(row, col) = val
#endif

  end subroutine set_jac_entry


  subroutine scale_jac_entry(state, row, col, val)

    !$acc routine seq

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    type (burn_t), intent(inout) :: state
    integer, intent(in) :: row, col
    real(rt), intent(in) :: val

    !$gpu

#ifdef REACT_SPARSE_JACOBIAN
    call scale_csr_jac_entry(state % sparse_jac, row, col, val)
#else
    state % jac(row, col) = state % jac(row, col) * val
#endif

  end subroutine scale_jac_entry  


  subroutine get_jac_entry(state, row, col, val)

    !$acc routine seq

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    type (burn_t), intent(in) :: state
    integer, intent(in) :: row, col
    real(rt), intent(out) :: val

    !$gpu

#ifdef REACT_SPARSE_JACOBIAN
    call get_csr_jac_entry(state % sparse_jac, row, col, val)
#else
    val = state % jac(row, col)
#endif

  end subroutine get_jac_entry


  subroutine set_jac_zero(state)

    !$acc routine seq

    use amrex_constants_module, only: ZERO

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    type (burn_t), intent(inout) :: state

    !$gpu

#ifdef REACT_SPARSE_JACOBIAN
    state % sparse_jac(:) = ZERO
#else
    state % jac(:,:) = ZERO
#endif

  end subroutine set_jac_zero

end module jacobian_sparsity_module

