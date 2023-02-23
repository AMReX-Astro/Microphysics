#ifndef HYBRJ_QFORM_H
#define HYBRJ_QFORM_H

#include <hybrj_type.H>

// this subroutine proceeds from the computed qr factorization of
// an m by n matrix a to accumulate the m by m orthogonal matrix
// q from its factored form.
//
//  argonne national laboratory. minpack project. march 1980.
//  burton s. garbow, kenneth e. hillstrom, jorge j. more

// changes:
// originally, this routine supported an mxn matrix, but we assume m=n

template<int neqs>
void qform(amrex::Array2D<Real, 1, neqs, 1, neqs>& q,
           amrex::Array1D<Real, 1, neqs>& wa) {

    // q is an n by n array. on input the full lower trapezoid in
    // the n columns of q contains the factored form.
    // on output q has been accumulated into a square matrix.

    // zero out upper triangle of q in the first n

    if (neqs >= 2) {
        for (int j = 2; j <= neqs; ++j) {
            int jm1 = j - 1;
            for (int i = 1; i <= jm1; ++i) {
                q(i,j) = 0.0_rt;
            }
        }
    }

    // accumulate q from its factored form.

    for (int l = 1; l <= neqs; ++l) {
        int k = neqs - l + 1;
        for (int i = k; i <= neqs; ++i) {
            wa(i) = q(i,k);
            q(i,k) = 0.0_rt;
        }
        q(k,k) = 1.0_rt;
        if (wa(k) != 0.0_rt) {
            for (int j = k; j <= neqs; ++j) {
                Real sum = 0.0_rt;
                for (int i = k; i <= neqs; ++i) {
                    sum += q(i,j) * wa(i);
                }
                Real temp = sum / wa(k);
                for (int i = k; i <= neqs; ++i) {
                    q(i,j) += -temp * wa(i);
                }
            }
        }
    }
}
#endif
