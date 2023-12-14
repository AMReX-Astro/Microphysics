#ifndef LINPACK_H
#define LINPACK_H

#include <AMReX_REAL.H>
#include <AMReX_Array.H>

#include <ArrayUtilities.H>

template <int num_eqs>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void dgesl (RArray2D& a1, IArray1D& pivot1, RArray1D& b1)
{
    int nm1 = num_eqs - 1;

    // solve a * x = b
    // first solve l * y = b
    if (nm1 >= 1) {
        for (int k = 0; k < nm1; ++k) {
            int l = pivot1(k+1) - 1;
            Real t = b1(l+1);
            if (l != k) {
                b1(l+1) = b1(k+1);
                b1(k+1) = t;
            }

            for (int j = k+1; j < num_eqs; ++j) {
                b1(j+1) += t * a1(j+1,k+1);
            }
        }
    }

    // now solve u * x = y
    for (int kb = 0; kb < num_eqs; ++kb) {

        int k = num_eqs - kb - 1;
        b1(k+1) = b1(k+1) / a1(k+1,k+1);
        Real t = -b1(k+1);
        for (int j = 0; j < k; ++j) {
            b1(j+1) += t * a1(j+1,k+1);
        }
    }

}



template <int num_eqs>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void dgefa (RArray2D& a1, IArray1D& pivot1, int& info)
{
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
            Real dmax = std::abs(a1(k+1,k+1));
            for (int i = k+1; i < num_eqs; ++i) {
                if (std::abs(a1(i+1,k+1)) > dmax) {
                    l = i;
                    dmax = std::abs(a1(i+1,k+1));
                }
            }

            pivot1(k+1) = l + 1;

            // zero pivot implies this column already triangularized
            if (a1(l+1,k+1) != 0.0e0_rt) {

                // interchange if necessary
                if (l != k) {
                    t = a1(l+1,k+1);
                    a1(l+1,k+1) = a1(k+1,k+1);
                    a1(k+1,k+1) = t;
                }

                // compute multipliers
                t = -1.0e0_rt / a1(k+1,k+1);
                for (int j = k+1; j < num_eqs; ++j) {
                    a1(j+1,k+1) *= t;
                }

                // row elimination with column indexing
                for (int j = k+1; j < num_eqs; ++j) {
                    t = a1(l+1,j+1);
                    if (l != k) {
                        a1(l+1,j+1) = a1(k+1,j+1);
                        a1(k+1,j+1) = t;
                    }
                    for (int i = k+1; i < num_eqs; ++i) {
                        a1(i+1,j+1) += t * a1(i+1,k+1);
                    }
                }
            }
            else {

                info = k+1;

            }

        }

    }

    pivot1(num_eqs) = num_eqs;

    if (a1(num_eqs,num_eqs) == 0.0e0_rt) {
        info = num_eqs;
    }

}

#endif
