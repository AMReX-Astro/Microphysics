#ifndef _vode_linpack_H_
#define _vode_linpack_H_

#include <vode_type.H>

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void dgesl (RArray2D& a, IArray1D& pivot, RArray1D& b)
{

    int nm1 = VODE_NEQS - 1;

    // solve a * x = b
    // first solve l * y = b
    if (nm1 >= 1) {
        for (int k = 1; k <= nm1; ++k) {
            int l = pivot(k);
            Real t = b(l);
            if (l != k) {
                b(l) = b(k);
                b(k) = t;
            }

            for (int j = k+1; j <= VODE_NEQS; ++j) {
                b(j) += t * a(j,k);
            }
        }
    }

    // now solve u * x = y
    for (int kb = 1; kb <= VODE_NEQS; ++kb) {

        int k = VODE_NEQS + 1 - kb;
        b(k) = b(k) / a(k,k);
        Real t = -b(k);
        for (int j = 1; j <= k-1; ++j) {
            b(j) += t * a(j,k);
        }
    }

}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void dgefa (RArray2D& a, IArray1D& pivot, int& info)
{

    // dgefa factors a matrix by gaussian elimination.
    // a is returned in the form a = l * u where
    // l is a product of permutation and unit lower
    // triangular matrices and u is upper triangular.

    // gaussian elimination with partial pivoting

    info = 0;
    int nm1 = VODE_NEQS - 1;

    Real t;

    if (nm1 >= 1) {

        for (int k = 1; k <= nm1; ++k) {

            // find l = pivot index
            int l = k;
            Real dmax = std::abs(a(k,k));
            for (int i = k+1; i <= VODE_NEQS; ++i) {
                if (std::abs(a(i,k)) > dmax) {
                    l = i;
                    dmax = std::abs(a(i,k));
                }
            }

            pivot(k) = l;

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
                for (int j = k+1; j <= VODE_NEQS; ++j) {
                    a(j,k) *= t;
                }

                // row elimination with column indexing
                for (int j = k+1; j <= VODE_NEQS; ++j) {
                    t = a(l,j);
                    if (l != k) {
                        a(l,j) = a(k,j);
                        a(k,j) = t;
                    }
                    for (int i = k+1; i <= VODE_NEQS; ++i) {
                        a(i,j) += t * a(i,k);
                    }
                }
            }
            else {

                info = k;

            }

        }

    }

    pivot(VODE_NEQS) = VODE_NEQS;

    if (a(VODE_NEQS,VODE_NEQS) == 0.0e0_rt) {
        info = VODE_NEQS;
    }

}

#endif
