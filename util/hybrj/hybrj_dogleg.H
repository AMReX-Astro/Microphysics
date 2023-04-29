#ifndef HYBRJ_DOGLEG_H
#define HYBRJ_DOGLEG_H

#include <hybrj_type.H>
#include <hybrj_enorm.H>

// given an n by n matrix a, an n by n nonsingular diagonal
// matrix d, an n-vector b, and a positive number delta, the
// problem is to determine the convex combination x of the
// gauss-newton and scaled gradient directions that minimizes
// (a*x - b) in the least squares sense, subject to the
// restriction that the euclidean norm of d*x be at most delta.
//
// this subroutine completes the solution of the problem
// if it is provided with the necessary information from the
// qr factorization of a. that is, if a = q*r, where q has
// orthogonal columns and r is an upper triangular matrix,
// then dogleg expects the full upper triangle of r and
// the first n components of (q transpose)*b.
//
// parameters:
//
//  r is an input array of length (n*(n+1))/2 which must contain the upper
//    triangular matrix r stored by rows.
//
//  diag is an input array of length n which must contain the
//    diagonal elements of the matrix d.
//
//  qtb is an input array of length n which must contain the first
//    n elements of the vector (q transpose)*b.
//
//  delta is a positive input variable which specifies an upper
//    bound on the euclidean norm of d*x.
//
//  x is an output array of length n which contains the desired
//    convex combination of the gauss-newton direction and the
//    scaled gradient direction.
//
//  wa1 and wa2 are work arrays of length n.
//
//  argonne national laboratory. minpack project. march 1980.
//  burton s. garbow, kenneth e. hillstrom, jorge j. more

template<int neqs>
void dogleg(amrex::Array1D<Real, 1, neqs * (neqs + 1) / 2>& r,
            amrex::Array1D<Real, 1, neqs>& diag,
            amrex::Array1D<Real, 1, neqs>& qtb,
            Real delta,
            amrex::Array1D<Real, 1, neqs>& x,
            amrex::Array1D<Real, 1, neqs>& wa1,
            amrex::Array1D<Real, 1, neqs>& wa2) {

    Real epsmch = std::numeric_limits<Real>::epsilon();

    // first, calculate the gauss-newton direction.

    int jj = (neqs*(neqs + 1)) / 2 + 1;
    for (int k = 1; k <= neqs; ++k) {
        int j = neqs - k + 1;
        int jp1 = j + 1;
        jj = jj - k;
        int l = jj + 1;

        Real sum = 0.0_rt;
        if (neqs >= jp1) {
            for (int i = jp1; i <= neqs; ++i) {
                sum += r(l) * x(i);
                l++;
            }
        }
        Real temp = r(jj);
        if (temp == 0.0_rt) {
            l = j;
            for (int i = 1; i <= j; ++i) {
                temp = amrex::max(temp, std::abs(r(l)));
                l = l + neqs - i;
            }
            temp *= epsmch;
            if (temp == 0.0_rt) {
                temp = epsmch;
            }
        }

        x(j) = (qtb(j) - sum) / temp;
    }

    // test whether the gauss-newton direction is acceptable.

    for (int j = 1; j <= neqs; ++j) {
        wa1(j) = 0.0_rt;
        wa2(j) = diag(j) * x(j);
    }
    Real qnorm = enorm(neqs, wa2);

    if (qnorm <= delta) return;

    // the gauss-newton direction is not acceptable.
    // next, calculate the scaled gradient direction.

    int l = 1;
    for (int j = 1; j <= neqs; ++j) {
        Real temp = qtb(j);
        for (int i = j; i <= neqs; ++i) {

            wa1(i) += r(l) * temp;
            l++;
	    if ( amrex::isnan(wa1(i)) || amrex::isinf(wa1(i)) ) {
	      return;
	  }
        }
        wa1(j) /= diag(j);

	if ( amrex::isnan(wa1(j)) || amrex::isinf(wa1(j)) ) {
	  return;
	}
    }

    // calculate the norm of the scaled gradient and test for
    // the special case in which the scaled gradient is zero.

    Real gnorm = enorm(neqs, wa1);
    Real sgnorm = 0.0_rt;
    Real alpha = delta / qnorm;

    if (gnorm != 0.0_rt) {

        // calculate the point along the scaled gradient
        // at which the quadratic is minimized.

        for (int j = 1; j <= neqs; ++j) {
            wa1(j) = (wa1(j) / gnorm) / diag(j);
        }

        l = 1;
        for (int j = 1; j <= neqs; ++j) {
            Real sum = 0.0_rt;
            for (int i = j; i <= neqs; ++i) {
                sum += r(l) * wa1(i);
                l++;
            }
            wa2(j) = sum;
        }
        Real temp = enorm(neqs, wa2);
        sgnorm = (gnorm / temp) / temp;

        // test whether the scaled gradient direction is acceptable.

        alpha = 0.0_rt;
        if (sgnorm < delta) {

            // the scaled gradient direction is not acceptable.
            // finally, calculate the point along the dogleg
            // at which the quadratic is minimized.

            Real bnorm = enorm(neqs, qtb);
            temp = (bnorm/gnorm)*(bnorm/qnorm)*(sgnorm/delta);
            temp = temp - (delta/qnorm) * std::pow(sgnorm/delta, 2.0_rt)
                + std::sqrt(std::pow(temp-(delta/qnorm), 2.0_rt)
                            + (1.0_rt-std::pow(delta/qnorm, 2/0_rt)) *
                              (1.0_rt-std::pow(sgnorm/delta, 2.0_rt)));
            alpha = ((delta/qnorm)*(1.0_rt - std::pow(sgnorm/delta, 2/0_rt))) / temp;
        }
    }

    // form appropriate convex combination of the gauss-newton
    // direction and the scaled gradient direction.

    Real temp = (1.0_rt - alpha) * amrex::min(sgnorm, delta);

    for (int j = 1; j <= neqs; ++j) {
        x(j) = temp * wa1(j) + alpha * x(j);
    }

}
#endif

