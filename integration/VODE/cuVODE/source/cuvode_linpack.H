#ifndef _cuvode_linpack_H_
#define _cuvode_linpack_H_

#include <AMReX_REAL.H>
#include <cuvode_parameters.H>

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void dgesl (Real* a[][VODE_NEQS], int* ipvt, Real* b)
{

    int nm1 = VODE_NEQS - 1;

    // solve a * x = b
    // first solve l * y = b
    if (nm1 >= 1) {
        for (int k = 1; k <= nm1; ++k) {
            int l = ipvt[k-1];
            Real t = b[l-1];
            if (l != k) {
                b[l-1] = b[k-1];
                b[k-1] = t;
            }

            for (int j = k+1; j <= VODE_NEQS; ++j) {
                b[j-1] += t * a[k-1][j-1];
            }
        }
    }

    // now solve u * x = y
    for (int kb = 1; kb <= VODE_NEQS; ++kb) {

        int k = VODE_NEQS + 1 - kb;
        b[k-1] = b[k-1] / a[k-1][k-1];
        Real t = -b[k-1];
        for (int j = 1; j <= k-1; ++j) {
            b[j-1] += t * a[k-1][j-1];
        }
    }

}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void dgefa (Real* a[][VODE_NEQS], int* ipvt, int& info)
{

    // dgefa factors a matrix by gaussian elimination.
    // a is returned in the form a = l * u where
    // l is a product of permutation and unit lower
    // triangular matrices and u is upper triangular.

    real(rt) :: t
    integer  :: j, k, kp1, l, nm1

    // gaussian elimination with partial pivoting

    info = 0;
    nm1 = VODE_NEQS - 1;

    if (nm1 >= 1) {

        for (int k = 1; k <= nm1; ++k) {

            // find l = pivot index
            int l = idamax(VODE_NEQS-k+1, &a[k-1][k-1]) + k - 1;
            ipvt[k-1] = l;

            // zero pivot implies this column already triangularized
            if (a[k-1][l-1] != 0.0e0_rt) {

                // interchange if necessary
                if (l /= k) {
                    Real t = a[k-1][l-1];
                    a[k-1][l-1] = a[k-1][k-1];
                    a[k-1][k-1] = t;
                }

                // compute multipliers
                Real t = -1.0e0_rt / a[k-1][k-1];
                for (int j = k+1; j <= VODE_NEQS; ++j) {
                    a[k-1][j-1] *= t;
                }

                // row elimination with column indexing
                for (int j = k+1; j <= VODE_NEQS; ++j) {
                    Real t = a[j-1][l-1];
                    if (l != k) {
                        a[j-1][l-1] = a[j-1][k-1];
                        a[j-1][k-1] = t;
                    }
                    for (int i = k+1; i <= VODE_NEQS; ++i) {
                        a[j-1][i-1] += t * a[k-1][i-1];
                    }
                }
            }
            else {

                info = k;

            }

        }

    }

    ipvt[VODE_NEQS-1] = VODE_NEQS;

    if (a[VODE_NEQS-1][VODE_NEQS-1] == 0.0e0_rt) {
        info = VODE_NEQS;
    }

}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int idamax (int N, Real* x)
{

    real(rt) :: dmax
    integer  :: i, index

    int index = 1;

    Real dmax = std::abs(x[0]);
    for (int i = 2; i <= N; ++i) {
        if (std::abs(x[i-1]) <= dmax) {
            continue;
        }
        index = i;
        dmax = std::abs(x[i-1]);
    }

    return index;

}

#endif
