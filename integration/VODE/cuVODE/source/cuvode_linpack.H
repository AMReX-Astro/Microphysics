#ifndef _cuvode_linpack_H_
#define _cuvode_linpack_H_

#include <AMReX_REAL.H>
#include <cuvode_parameters.H>

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void dgesl (Real* a, int* ipvt, Real* b)
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
                b[j-1] += t * a[(k-1) * VODE_NEQS + j-1];
            }
        }
    }

    // now solve u * x = y
    for (int kb = 1; kb <= VODE_NEQS; ++kb) {

        int k = VODE_NEQS + 1 - kb;
        b[k-1] = b[k-1] / a[(k-1) * VODE_NEQS + k-1];
        Real t = -b[k-1];
        for (int j = 1; j <= k-1; ++j) {
            b[j-1] += t * a[(k-1) * VODE_NEQS + j-1];
        }
    }

}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int idamax (int N, Real* x)
{

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


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void dgefa (Real* a, int* ipvt, int& info)
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
            int l = idamax(VODE_NEQS-k+1, &a[(k-1) * VODE_NEQS + k-1]) + k - 1;
            ipvt[k-1] = l;

            // zero pivot implies this column already triangularized
            if (a[(k-1) * VODE_NEQS + l-1] != 0.0e0_rt) {

                // interchange if necessary
                if (l != k) {
                    t = a[(k-1) * VODE_NEQS + l-1];
                    a[(k-1) * VODE_NEQS + l-1] = a[(k-1) * VODE_NEQS + k-1];
                    a[(k-1) * VODE_NEQS + k-1] = t;
                }

                // compute multipliers
                t = -1.0e0_rt / a[(k-1) * VODE_NEQS + k-1];
                for (int j = k+1; j <= VODE_NEQS; ++j) {
                    a[(k-1) * VODE_NEQS + j-1] *= t;
                }

                // row elimination with column indexing
                for (int j = k+1; j <= VODE_NEQS; ++j) {
                    t = a[(j-1) * VODE_NEQS + l-1];
                    if (l != k) {
                        a[(j-1) * VODE_NEQS + l-1] = a[(j-1) * VODE_NEQS + k-1];
                        a[(j-1) * VODE_NEQS + k-1] = t;
                    }
                    for (int i = k+1; i <= VODE_NEQS; ++i) {
                        a[(j-1) * VODE_NEQS + i-1] += t * a[(k-1) * VODE_NEQS + i-1];
                    }
                }
            }
            else {

                info = k;

            }

        }

    }

    ipvt[VODE_NEQS-1] = VODE_NEQS;

    if (a[(VODE_NEQS-1) * VODE_NEQS + VODE_NEQS-1] == 0.0e0_rt) {
        info = VODE_NEQS;
    }

}

#endif
