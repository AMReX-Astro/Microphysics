#ifndef HYBRJ_TYPE_H
#define HYBRJ_TYPE_H

#include <AMReX_REAL.H>
#include <AMReX_Array.H>

using namespace amrex;

// changes compared to Minpack's hybrj:
//
// 1. we assume that the Jacobian is square, of shape n x n, hence
//    no need for ldfjac
//
// 2. nprint was removed

// maxfev is a positive integer input variable. termination
// occurs when the number of calls to fcn with iflag = 1
// has reached maxfev.
const int maxfev = 1000;

// factor is a positive input variable used in determining the
// initial step bound. this bound is set to the product of
// factor and the euclidean norm of diag*x if nonzero, or else
// to factor itself. in most cases factor should lie in the
// interval (.1,100.). 100. is a generally recommended value.
const Real factor = 100.0_rt;

template<int neqs>
struct hybrj_t {

    static_assert(neqs > 0);

    // n is a positive integer input variable set to the number
    // of functions and variables.
    const int n = neqs;

    // x is an array of length n. on input x must contain
    // an initial estimate of the solution vector. on output x
    // contains the final estimate of the solution vector.
    amrex::Array1D<Real, 1, neqs> x;

    // fvec is an output array of length n which contains
    // the functions evaluated at the output x.
    amrex::Array1D<Real, 1, neqs> fvec;

    // fjac is an output n by n array which contains the
    // orthogonal matrix q produced by the qr factorization
    // of the final approximate jacobian.
    amrex::Array2D<Real, 1, neqs, 1, neqs> fjac;

    // xtol is a nonnegative input variable. termination
    // occurs when the relative error between two consecutive
    // iterates is at most xtol.
    Real xtol;

    // diag is an array of length n. if mode = 1 (see
    // below), diag is internally set. if mode = 2, diag
    // must contain positive entries that serve as
    // multiplicative scale factors for the variables.
    amrex::Array1D<Real, 1, neqs> diag;

    // mode is an integer input variable. if mode = 1, the
    // variables will be scaled internally. if mode = 2,
    // the scaling is specified by the input diag. other
    // values of mode are equivalent to mode = 1.
    int mode;

    // info is an integer output variable. if the user has
    // terminated execution, info is set to the (negative)
    // value of iflag. see description of fcn. otherwise,
    // info is set as follows.
    //
    // info = 0   improper input parameters.
    //
    // info = 1   relative error between two consecutive iterates
    //            is at most xtol.
    //
    // info = 2   number of calls to fcn with iflag = 1 has
    //            reached maxfev.
    //
    // info = 3   xtol is too small. no further improvement in
    //            the approximate solution x is possible.
    //
    // info = 4   iteration is not making good progress, as
    //            measured by the improvement from the last
    //            five jacobian evaluations.
    //
    // info = 5   iteration is not making good progress, as
    //            measured by the improvement from the last
    //            ten iterations.
    int info;

    // nfev is an integer output variable set to the number of
    // calls to fcn with iflag = 1.
    int nfev;

    // njev is an integer output variable set to the number of
    // calls to fcn with iflag = 2.
    int njev;

    // r is an output array of length lr which contains the
    // upper triangular matrix produced by the qr factorization
    // of the final approximate jacobian, stored rowwise.
    amrex::Array1D<Real, 1, neqs * (neqs + 1) / 2> r;

    // qtf is an output array of length n which contains
    // the vector (q transpose)*fvec.
    amrex::Array1D<Real, 1, neqs> qtf;

    // wa1, wa2, wa3, and wa4 are work arrays of length n.
    amrex::Array1D<Real, 1, neqs> wa1;
    amrex::Array1D<Real, 1, neqs> wa2;
    amrex::Array1D<Real, 1, neqs> wa3;
    amrex::Array1D<Real, 1, neqs> wa4;
};

#endif
