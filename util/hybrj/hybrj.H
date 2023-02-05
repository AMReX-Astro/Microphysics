#ifndef HYBRJ_H
#define HYBRJ_H

#include <hybrj_type.H>

#include <hybrj_qrfac.H>
#include <hybrj_qform.H>
#include <hybrj_dogleg.H>
#include <hybrj_enorm.H>
#include <hybrj_r1mpyq.H>
#include <hybrj_r1updt.H>
#include <AMReX_GpuUtility.H>

template<int neqs, typename T>
void hybrj(hybrj_t<neqs>& hj,
           const T& data,
           std::function<void(amrex::Array1D<Real, 1, neqs>&,
                              amrex::Array1D<Real, 1, neqs>&,
                              const T&, int&)> fcn,
           std::function<void(amrex::Array1D<Real, 1, neqs>&,
                              amrex::Array2D<Real, 1, neqs, 1, neqs>&,
                              const T&, int&)> jcn) {

    // the purpose of hybrj is to find a zero of a system of
    // n nonlinear functions in n variables by a modification
    // of the powell hybrid method. the user must provide a
    // subroutine which calculates the functions and the jacobian.
    //
    // fcn jcn are the names of the user-supplied subroutine which
    // calculates the functions and the jacobian.
    //
    // void fcn(Array1D<Real, 1, neqs>& x, Array1D<Real, 1, neqs>& fvec, int& iflag)
    // void jcn(Array1D<Real, 1, neqs>& x, Array2D<Real, 1, neqs, 1, neqs> fjac, int& iflag)
    //
    // the value of iflag should not be changed by fcn unless
    // the user wants to terminate execution of hybrj.
    // in this case set iflag to a negative integer.

    bool finished = false;

    Real epsmch = std::numeric_limits<Real>::epsilon();

    hj.info = 0;
    int iflag = 0;
    hj.nfev = 0;
    hj.njev = 0;

    Real delta = 0.0_rt;
    Real xnorm = 0.0_rt;

    // check the input parameters for errors.

    if (hj.xtol < 0.0_rt) {
        amrex::Error("error: xtol must be > 0");
    }

    if (hj.mode == 2) {
        for (int j = 1; j <= neqs; ++j) {
            if (hj.diag(j) <= 0.0_rt) {
                finished = true;
                break;
            }
        }
    }

    // evaluate the function at the starting point
    // and calculate its norm.
    Real fnorm = 0.0_rt;

    iflag = 1;

    fcn(hj.x, hj.fvec, data, iflag);
    hj.nfev = 1;
    if (iflag < 0) {
        // user requested abort
        finished = true;
    } else {
        fnorm = enorm(neqs, hj.fvec);
    }

    // initialize iteration counter and monitors.

    int iter = 1;
    int ncsuc = 0;
    int ncfail = 0;
    int nslow1 = 0;
    int nslow2 = 0;

    // beginning of the outer loop.

    while (! finished) {

        bool jeval = true;

        // calculate the jacobian matrix.

        jcn(hj.x, hj.fjac, data, iflag);
        hj.njev++;
        if (iflag < 0) {
            finished = true;
            break;
        }

        // compute the qr factorization of the jacobian.

        qrfac(hj.fjac, hj.wa1, hj.wa2, hj.wa3);

        // on the first iteration and if mode is 1, scale according
        // to the norms of the columns of the initial jacobian.

        if (iter == 1) {
            if (hj.mode != 2) {
                for (int j = 1; j <= hj.n; ++j) {
                    hj.diag(j) = hj.wa2(j);
                    if (hj.wa2(j) == 0.0_rt) {
                        hj.diag(j) = 1.0_rt;
                    }
                }
            }

            // on the first iteration, calculate the norm of the scaled x
            // and initialize the step bound delta.

            for (int j = 1; j <= neqs; ++j) {
                hj.wa3(j) = hj.diag(j) * hj.x(j);
            }
            xnorm = enorm(hj.n, hj.wa3);
            delta =  factor * xnorm;
            if (delta == 0.0_rt) {
                delta = factor;
            }
        }

        // form (q transpose)*fvec and store in qtf.

        for (int i = 1; i <= hj.n; ++i) {
            hj.qtf(i) = hj.fvec(i);
        }
        for (int j = 1; j <= neqs; ++j) {
            if (hj.fjac(j,j) != 0.0_rt) {
                Real sum = 0.0_rt;
                for (int i = j; i <= hj.n; ++i) {
                    sum += hj.fjac(i,j) * hj.qtf(i);
                }
                Real temp = -sum / hj.fjac(j,j);
                for (int i = j; i <= hj.n; ++i) {
                    hj.qtf(i) += hj.fjac(i,j) * temp;
                }
            }
        }

        // copy the triangular factor of the qr factorization into r.

        bool sing = false;

        for (int j = 1; j <= hj.n; ++j) {
            int l = j;
            int jm1 = j - 1;
            if (jm1 >= 1) {
                for (int i = 1; i <= jm1; ++i) {
                    hj.r(l) = hj.fjac(i,j);
                    l = l + neqs - i;
                }
            }
            hj.r(l) = hj.wa1(j);
            if (hj.wa1(j) == 0.0_rt) {
                sing = true;
            }
        }

        // accumulate the orthogonal factor in fjac.

        qform(hj.fjac, hj.wa1);

        // rescale if necessary.

        if (hj.mode != 2) {
            for (int j = 1; j <= hj.n; ++j) {
                hj.diag(j) = amrex::max(hj.diag(j), hj.wa2(j));
            }
        }

        // beginning of the inner loop.

        while (true) {

            // determine the direction p.

            dogleg(hj.r, hj.diag, hj.qtf, delta, hj.wa1, hj.wa2, hj.wa3);

            // store the direction p and x + p. calculate the norm of p.

	    bool invalid = false;

	    for (int i = 1; i < neqs; ++i) {
	      if ( amrex::isinf(hj.wa1(i)) || amrex::isnan(hj.wa1(i)) ){
		invalid = true;
	      }
	    }

	    if (invalid) {
	      finished = true;
	      break;
	    }

            for (int j = 1; j <= neqs; ++j) {
                hj.wa1(j) = -hj.wa1(j);
                hj.wa2(j) = hj.x(j) + hj.wa1(j);
                hj.wa3(j) = hj.diag(j) * hj.wa1(j);
            }

            Real pnorm = enorm(hj.n, hj.wa3);

            // on the first iteration, adjust the initial step bound.

            if (iter == 1) {
                delta = amrex::min(delta, pnorm);
            }

            // evaluate the function at x + p and calculate its norm.

            fcn(hj.wa2, hj.wa4, data, iflag);
            hj.nfev++;

            if (iflag < 0) {
                finished = true;
                break;
            }

            Real fnorm1 = enorm(hj.n, hj.wa4);

            // compute the scaled actual reduction.

            Real actred = -1.0_rt;
            if (fnorm1 < fnorm) {
                actred = 1.0_rt - std::pow(fnorm1/fnorm, 2.0_rt);
            }

            // compute the scaled predicted reduction.

            int l = 1;
            for (int i = 1; i <= hj.n; ++i) {
                Real sum = 0.0_rt;
                for (int j = i; j <= neqs; ++j) {
                    sum += hj.r(l) * hj.wa1(j);
                    l += 1;
                }
                hj.wa3(i) = hj.qtf(i) + sum;
            }
            Real temp = enorm(hj.n, hj.wa3);
            Real prered = 0.0_rt;
            if (temp < fnorm) {
                prered = 1.0_rt - std::pow(temp/fnorm, 2.0_rt);
            }

            // compute the ratio of the actual to the predicted
            // reduction.

            Real ratio = 0.0_rt;
            if (prered > 0.0_rt) {
                ratio = actred / prered;
            }

            // update the step bound.

            if (ratio < 0.1_rt) {
                ncsuc = 0;
                ncfail += 1;
                delta = 0.5_rt * delta;
            } else {
                ncfail = 0;
                ncsuc += 1;
                if (ratio >= 0.5_rt || ncsuc > 1) {
                    delta = amrex::max(delta, pnorm / 0.5_rt);
                }
                if (std::abs(ratio-1.0_rt) <= 0.1_rt) {
                    delta = pnorm / 0.5_rt;
                }
            }

            // test for successful iteration.

            if (ratio >= 1.e-4_rt) {

                // successful iteration. update x, fvec, and their norms.

                for (int j = 1; j <= hj.n; ++j) {
                    hj.x(j) = hj.wa2(j);
                    hj.wa2(j) = hj.diag(j) * hj.x(j);
                    hj.fvec(j) = hj.wa4(j);
                }

                xnorm = enorm(hj.n, hj.wa2);
                fnorm = fnorm1;
                iter++;
            }

            // determine the progress of the iteration.

            nslow1++;
            if (actred >= 1.e-3_rt) {
                nslow1 = 0;
            }
            if (jeval) {
                nslow2++;
            }
            if (actred >= 0.1_rt) {
                nslow2 = 0;
            }

            // test for convergence.

            if (delta <= hj.xtol * xnorm || fnorm == 0.0_rt) {
                hj.info = 1;
            }
            if (hj.info != 0) {
                finished = true;
                break;
            }

            // tests for termination and stringent tolerances.

            if (hj.nfev >= maxfev) {
                hj.info = 2;
            }
            if (0.1_rt * amrex::max(0.1_rt*delta, pnorm) <= epsmch * xnorm) {
                hj.info = 3;
            }
            if (nslow2 == 5) {
                hj.info = 4;
            }
            if (nslow1 == 10) {
                hj.info = 5;
            }

            if (hj.info != 0) {
                finished = true;
                break;
            }

            // criterion for recalculating jacobian.

            if (ncfail == 2) {
                break;
            }

            // calculate the rank one modification to the jacobian
            // and update qtf if necessary.
	    if (pnorm == 0.0_rt) {
	      break;
	    }

            for (int j = 1; j <= hj.n; ++j) {
                Real sum = 0.0_rt;
                for (int i = 1; i <= hj.n; ++i) {
                    sum += hj.fjac(i,j) * hj.wa4(i);
                }
                hj.wa2(j) = (sum - hj.wa3(j)) / pnorm;
                hj.wa1(j) = hj.diag(j) * ((hj.diag(j) * hj.wa1(j)) / pnorm);
                if (ratio >= 1.e-4_rt) {
                    hj.qtf(j) = sum;
                }
            }

            // compute the qr factorization of the updated jacobian.

            r1updt(hj.r, hj.wa1, hj.wa2, hj.wa3, sing);
            r1mpyq(hj.fjac, hj.wa2, hj.wa3);
            r1mpyq(hj.qtf, hj.wa2, hj.wa3);

            // end of the inner loop.

            jeval = false;
        }

        // end of the outer loop

        if (finished) {
            break;
        }
    }


    // termination, either normal or user imposed.

    if (iflag < 0) {
        hj.info = iflag;
    }
    iflag = 0;

}
#endif
