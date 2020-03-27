#ifndef _jacobian_sparsity_H_
#define _jacobian_sparsity_H_

#include <AMReX_REAL.H>
#include <AMReX_Array.H>

#include <ArrayUtilities.H>
#include <network.H>

using namespace ArrayUtil;

#ifdef REACT_SPARSE_JACOBIAN
AMREX_GPU_HOST_DEVICE AMREX_INLINE
int lookup_csr_jac_loc(const int row, const int col)
{
    // return value: the flattened CSR location of (row, col)
    int csr_loc;

    int num_in_row, row_start_loc, row_end_loc;

    /*
    Looks up the index into a CSR-formatted Jacobian
    matrix given row and col indices into the
    equivalent dense matrix.

    Assumes the base in first element of CSR row count array is 1
    */

    num_in_row = csr_jac_row_count(row+1) - csr_jac_row_count(row);
    row_start_loc = csr_jac_row_count(row);
    row_end_loc   = row_start_loc + num_in_row - 1;

    csr_loc = -1;
    for (int i = row_start_loc; i <= row_end_loc; ++i) {
       if (csr_jac_col_index(i) == col) {
          csr_loc = i;
          break;
       }
    }

    return csr_loc;
}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
void set_csr_jac_entry(RealArray1D<NETWORK_SPARSE_JAC_NNZ>& csr_jac,
                       const int row, const int col, const Real val)
{
    // Get index into the CSR Jacobian
    const int csr_loc = lookup_csr_jac_loc(row, col, csr_loc);

    // Set value in CSR Jacobian only if row, col entry exists
    AMREX_ASSERT(csr_loc != -1);
    csr_jac(csr_loc) = val;
}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
void scale_csr_jac_entry(RealArray1D<NETWORK_SPARSE_JAC_NNZ>& csr_jac,
                         const int row, const int col, const Real val)
{
    // Get index into the CSR Jacobian
    const int csr_loc = lookup_csr_jac_loc(row, col);

    // Scale value in CSR Jacobian only if row, col entry exists
    AMREX_ASSERT(csr_loc != -1);
    csr_jac(csr_loc) = csr_jac(csr_loc) * val;
}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real get_csr_jac_entry(RealArray1D<NETWORK_SPARSE_JAC_NNZ> const& csr_jac,
                       const int row, const int col)
{
    // Get index into the CSR Jacobian
    const int csr_loc = lookup_csr_jac_loc(row, col);

    // Return value in CSR Jacobian only if row, col entry exists
    AMREX_ASSERT(csr_loc != -1);
    return csr_jac(csr_loc);
}

template <class T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void set_jac_entry(T& jac, const int row, const int col, const Real val)
{
#ifdef REACT_SPARSE_JACOBIAN
    set_csr_jac_entry(jac, row, col, val);
#else
    jac(row, col) = val;
#endif
}

template <class T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void scale_jac_entry(T& jac, const int row, const int col, const Real val)
{
#ifdef REACT_SPARSE_JACOBIAN
    scale_csr_jac_entry(jac, row, col, val);
#else
    jac(row, col) = jac(row, col) * val;
#endif
}

template <class T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real get_jac_entry(T& jac, const int row, const int col)
{
#ifdef REACT_SPARSE_JACOBIAN
    return get_csr_jac_entry(jac, row, col);
#else
    return jac(row, col);
#endif
}


template <class T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void set_jac_zero(T& jac)
{
#ifdef REACT_SPARSE_JACOBIAN
    for (int i = 0; i < NETWORK_SPARSE_JAC_NNZ; ++i)
        jac(i) = 0._rt;
#else
    for (int j = 0; j < Burn::njcols; ++j)
        for (int i = 0; i < Burn::njrows; ++i)
            jac(i,j) = 0._rt;
#endif
}

#endif
