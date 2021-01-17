#ifndef LINPACK_H
#define LINPACK_H

#include <AMReX_REAL.H>
#include <AMReX_Array.H>

#include <ArrayUtilities.H>

template <int num_eqs>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void dgesl (RArray2D& a1, IArray1D& pivot1, RArray1D& b1)
{
    auto const& a = reinterpret_cast<Array2D<Real, 0, num_eqs-1, 0, num_eqs-1>&>(a1);
    auto const& pivot = reinterpret_cast<Array1D<int, 0, num_eqs-1>&>(pivot1);
    auto& b = reinterpret_cast<Array1D<Real, 0, num_eqs-1>&>(b1);

    int nm1 = num_eqs - 1;

    // solve a * x = b
    // first solve l * y = b
    if (nm1 >= 1) {
        for (int k = 0; k < nm1; ++k) {
            int l = pivot(k) - 1;
            Real t = b(l);
            if (l != k) {
                b(l) = b(k);
                b(k) = t;
            }

            for (int j = k+1; j < num_eqs; ++j) {
                b(j) += t * a(j,k);
            }
        }
    }

    // now solve u * x = y
    for (int kb = 0; kb < num_eqs; ++kb) {

        int k = num_eqs - kb - 1;
        b(k) = b(k) / a(k,k);
        Real t = -b(k);
        for (int j = 0; j < k; ++j) {
            b(j) += t * a(j,k);
        }
    }

}



template <int num_eqs>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void dgefa (RArray2D& a1, IArray1D& pivot1, int& info)
{
    auto& a = reinterpret_cast<Array2D<Real, 0, num_eqs-1, 0, num_eqs-1>&>(a1);
    auto& pivot = reinterpret_cast<Array1D<int, 0, num_eqs-1>&>(pivot1);

    // dgefa factors a matrix by gaussian elimination.
    // a is returned in the form a = l * u where
    // l is a product of permutation and unit lower
    // triangular matrices and u is upper triangular.

    // gaussian elimination with partial pivoting

    info = 0;
    int nm1 = num_eqs - 1;

    Real t;

    if (nm1 >= 1) {

        for (int k = 0; k < nm1; ++k) {

            // find l = pivot index
            int l = k;
            Real dmax = std::abs(a(k,k));
            for (int i = k+1; i < num_eqs; ++i) {
                if (std::abs(a(i,k)) > dmax) {
                    l = i;
                    dmax = std::abs(a(i,k));
                }
            }

            pivot(k) = l + 1;

            // zero pivot implies this column already triangularized
            if (a(l,k) != 0.0e0_rt) {

                // interchange if necessary
                if (l != k) {
                    t = a(l,k);
                    a(l,k) = a(k,k);
                    a(k,k) = t;
                }

                // compute multipliers
                t = -1.0e0_rt / a(k,k);
                for (int j = k+1; j < num_eqs; ++j) {
                    a(j,k) *= t;
                }

                // row elimination with column indexing
                for (int j = k+1; j < num_eqs; ++j) {
                    t = a(l,j);
                    if (l != k) {
                        a(l,j) = a(k,j);
                        a(k,j) = t;
                    }
                    for (int i = k+1; i < num_eqs; ++i) {
                        a(i,j) += t * a(i,k);
                    }
                }
            }
            else {

                info = k+1;

            }

        }

    }

    pivot(num_eqs-1) = num_eqs;

    if (a(num_eqs-1,num_eqs-1) == 0.0e0_rt) {
        info = num_eqs;
    }

}

#endif
