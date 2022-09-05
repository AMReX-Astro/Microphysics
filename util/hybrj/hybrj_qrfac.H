#ifndef HYBRJ_QRFAC_H
#define HYBRJ_QRFAC_H

#include <hybrj_type.H>
#include <hybrj_enorm.H>

// this subroutine uses householder transformations with column
// pivoting (optional) to compute a qr factorization of the
// m by n matrix a. that is, qrfac determines an orthogonal
// matrix q, a permutation matrix p, and an upper trapezoidal
// matrix r with diagonal elements of nonincreasing magnitude,
// such that a*p = q*r. the householder transformation for
// column k, k = 1,2,...,min(m,n), is of the form
//                        t
//        i - (1/u(k))*u*u
//
// where u has zeros in the first k-1 positions. the form of
// this transformation and the method of pivoting first
// appeared in the corresponding linpack subroutine.

// changes from the Fortran version:
//
// * removed the ability to pivot
// * assume input is square
// * written to use our hybrj_t type

// hybrj_t.fjac is an n by n array. on input a contains the matrix for
// which the qr factorization is to be computed. on output
// the strict upper trapezoidal part of a contains the strict
// upper trapezoidal part of r, and the lower trapezoidal
// part of a contains a factored form of q (the non-trivial
// elements of the u vectors described above).
//
// the following are assigned to different work arrays in the hybrj type
//
// rdiag is an output array of length n which contains the
// diagonal elements of r.
//
// acnorm is an output array of length n which contains the
// norms of the corresponding columns of the input matrix a.
// if this information is not needed, then acnorm can coincide
// with rdiag.
//
// wa is a work array of length n. if pivot is false, then wa
// can coincide with rdiag.
//
//
// argonne national laboratory. minpack project. march 1980.
// burton s. garbow, kenneth e. hillstrom, jorge j. more

template<int neqs>
void qrfac(Array2D<Real, 1, neqs, 1, neqs>& a,
           Array1D<Real, 1, neqs>& rdiag,
           Array1D<Real, 1, neqs>& acnorm,
           Array1D<Real, 1, neqs>& wa) {

    // temporary storage for a column of the Jacobian -- the Fortran
    // code would pass to enorm something like a(i, j) to compute the
    // norm of column j from row i to the last row.

    amrex::Array1D<Real, 1, neqs> tmp;

    // compute the initial column norms and initialize several arrays.

    for (int j = 1; j <= neqs; ++j) {
        for (int irow = 1; irow <= neqs; ++irow) {
            tmp(irow) = a(irow, j);
        }
        acnorm(j) = enorm(neqs, tmp);
        rdiag(j) = acnorm(j);
        wa(j) = rdiag(j);
    }


    // reduce a to r with householder transformations.

    for (int j = 1; j <= neqs; ++j) {

        // compute the householder transformation to reduce the
        // j-th column of a to a multiple of the j-th unit vector.

        // we want the norm of a(j,j:n)
        for (int irow = j; irow <= neqs; ++irow) {
            tmp(irow-j+1) = a(irow,j);
        }
        Real ajnorm = enorm(neqs-j+1, tmp);
        if (ajnorm != 0.0_rt) {
            if (a(j,j) < 0.0_rt) {
                ajnorm = -ajnorm;
            }
            for (int irow = j; irow <= neqs; ++irow) {
                a(irow,j) /= ajnorm;
            }
            a(j,j) += 1.0_rt;

            // apply the transformation to the remaining columns
            // and update the norms.

            int jp1 = j + 1;
            if (neqs >= jp1) {
               for (int k = jp1; k <= neqs; ++k) {
                   Real sum = 0.0_rt;
                   for (int i = j; i <= neqs; ++i) {
                       sum += a(i,j) * a(i,k);
                   }
                   Real temp = sum / a(j,j);
                   for (int i = j; i <= neqs; ++i) {
                       a(i,k) += -temp * a(i,j);
                   }
               }
            }
        }
        rdiag(j) = -ajnorm;
    }

}
#endif
