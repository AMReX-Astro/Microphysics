#ifndef ROCK4_H
#define ROCK4_H

#include <rock4_type.H>
#include <burn_type.H>
#ifdef STRANG
#include <integrator_type_strang.H>
#include <integrator_rhs_strang.H>
#endif
#ifdef SIMPLIFIED_SDC
#include <integrator_type_simplified_sdc.H>
#include <integrator_rhs_simplified_sdc.H>
#endif
#include <initial_timestep.H>
#include <rock4_coefficients.H>

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mdegre(rock4_t& rstate)
{

    constexpr amrex::Array1D<int, 1, 50> ms = {
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
    22, 24, 26, 28, 30, 32, 34, 36, 38, 41, 44, 47, 50, 53, 56, 59, 63, 67, 71,
    76, 81, 86, 92, 98, 105, 112, 120, 129, 138, 148};


    // Find the optimal degree.
    // MP(1): pointer which select the degree in ms(i)\1,2,..
    //        such that mdeg<=ms(i).
    // MP(2): pointer which gives the corresponding position
    //        of a_1 in the data recf for the selected degree.

    rstate.mp(2) = 1;

    // Find the degree

    for (int i = 1; i <= 50; ++i) {
        if ((ms(i) / rstate.mdeg) >= 1) {
            rstate.mdeg = ms(i);
            rstate.mp(1) = i;
            return;
        }
        rstate.mp(2) += ms(i)*2 - 1;
    }
    return;
}


template <typename BurnT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real rfstep (BurnT& state, rock4_t& rstate)
{

    // Solut. at t+h by an explicit (mdeg+4)-stages formula.

    int mz = rstate.mp(1);
    int mr = rstate.mp(2);

    // First stage

    Real temp1 = rstate.h * recf(mr);
    Real ci1 = rstate.t + temp1;
    Real ci2 = rstate.t + temp1;
    Real ci3 = rstate.t;

    Real temp2, temp3, temp4, temp5;

    for (int i = 1; i <= INT_NEQS; ++i) {
        rstate.yjm2(i) = rstate.yn(i);
        rstate.yjm1(i) = rstate.yn(i) + temp1 * rstate.fn(i);
        if (rstate.mdeg < 2) {
            rstate.y(i) = rstate.yjm1(i);
        }
    }

    // Stage for j=2..mdeg.

    for (int i = 2; i <= rstate.mdeg; ++i) {
        temp1 = rstate.h * recf(mr+2*(i-2)+1);
        temp3 = -recf(mr+2*(i-2)+2);
        temp2 = 1.0_rt - temp3;
        // we want to use rstate.yjm1 as input and store ydot in the original y
        // so we simply need to copy rstate.yjm1 to y here and use y_backup
        // to store the RHS and then copy it after the call
        for (int j = 1; j <= INT_NEQS; ++j) {
            rstate.y(j) = rstate.yjm1(j);
        }
        rhs(ci1, state, rstate, rstate.y_backup);
        for (int j = 1; j <= INT_NEQS; ++j) {
            rstate.y(j) = rstate.y_backup(j);
        }
        ci1 = temp1 + temp2 * ci2 + temp3 * ci3;
        for (int j = 1; j <= INT_NEQS; ++j) {
            rstate.y(j) = temp1 * rstate.y(j) + temp2 * rstate.yjm1(j) + temp3 * rstate.yjm2(j);

            // Shift the value "y" for the next stage.
            if (i < rstate.mdeg) {
                rstate.yjm2(j) = rstate.yjm1(j);
                rstate.yjm1(j) = rstate.y(j);
            }
        }
        ci3 = ci2;
        ci2 = ci1;
    }

    // The finishing procedure (4-stage method)

    // stage 1
    temp1 = rstate.h * fpa(mz, 1);
    rhs(ci1, state, rstate, rstate.yjm1);
    for (int j = 1; j <= INT_NEQS; ++j) {
        rstate.yjm3(j) = rstate.y(j) + temp1 * rstate.yjm1(j);
    }

    // stage 2
    ci2 = ci1 + temp1;
    temp1 = rstate.h * fpa(mz,2);
    temp2 = rstate.h * fpa(mz,3);

    // call the RHS with yjm3 as input and store the output in yjm2
    for (int i = 1; i <= INT_NEQS; ++i) {
        rstate.y_backup(i) = rstate.y(i);
        rstate.y(i) = rstate.yjm3(i);
    }
    rhs(ci2, state, rstate, rstate.yjm2);
    for (int i = 1; i <= INT_NEQS; ++i) {
        rstate.y(i) = rstate.y_backup(i);
    }

    for (int j = 1; j <= INT_NEQS; ++j) {
        rstate.yjm4(j) = rstate.y(j) + temp1 * rstate.yjm1(j) + temp2 * rstate.yjm2(j);
    }

    // stage 3
    ci2 = ci1 + temp1 + temp2;
    temp1 = rstate.h * fpa(mz,4);
    temp2 = rstate.h * fpa(mz,5);
    temp3 = rstate.h * fpa(mz,6);

    // call the RHS with yjm4 as input and store the output in yjm3
    for (int i = 1; i <= INT_NEQS; ++i) {
        // this is still backed-up: rstate.y_backup(i) = rstate.y(i);
        rstate.y(i) = rstate.yjm4(i);
    }
    rhs(ci2, state, rstate, rstate.yjm3);
    for (int i = 1; i <= INT_NEQS; ++i) {
        rstate.y(i) = rstate.y_backup(i);
    }

    for (int j = 1; j <= INT_NEQS; ++j) {
        rstate.fnt(j) = rstate.y(j) + temp1 * rstate.yjm1(j) + temp2 * rstate.yjm2(j) + temp3 * rstate.yjm3(j);
    }

    // stage 4
    ci2 = ci1 + temp1 + temp2 + temp3;
    temp1 = rstate.h * fpb(mz,1);
    temp2 = rstate.h * fpb(mz,2);
    temp3 = rstate.h * fpb(mz,3);
    temp4 = rstate.h * fpb(mz,4);

    // call the RHS with fnt as the input and store the output in yjm4
    for (int i = 1; i <= INT_NEQS; ++i) {
        // this is still backed-up: rstate.y_backup(i) = rstate.y(i);
        rstate.y(i) = rstate.fnt(i);
    }
    rhs(ci2, state, rstate, rstate.yjm4);
    for (int i = 1; i <= INT_NEQS; ++i) {
        rstate.y(i) = rstate.y_backup(i);
    }

    for (int j = 1; j <= INT_NEQS; ++j) {
        rstate.y(j) = rstate.y(j) + temp1 * rstate.yjm1(j) + temp2 * rstate.yjm2(j)
                                  + temp3 * rstate.yjm3(j) + temp4 * rstate.yjm4(j);
    }

    // Error evaluation (embedded method of order 3).
    temp1 = rstate.h * fpbe(mz,1) - temp1;
    temp2 = rstate.h * fpbe(mz,2) - temp2;
    temp3 = rstate.h * fpbe(mz,3) - temp3;
    temp4 = rstate.h * fpbe(mz,4) - temp4;
    temp5 = rstate.h * fpbe(mz,5);

    rhs(rstate.t + rstate.h, state, rstate, rstate.fnt);


    Real err{};

    for (int j = 1; j <= INT_NEQS; ++j) {
        Real wt{};
        if (j <= NumSpec) {
            wt = rstate.rtol_spec * std::abs(rstate.yn(j)) + rstate.atol_spec;
        } else {
            wt = rstate.rtol_enuc * std::abs(rstate.yn(j)) + rstate.atol_enuc;
        }

        err += std::pow((temp1 * rstate.yjm1(j) +
                         temp2 * rstate.yjm2(j) +
                         temp3 * rstate.yjm3(j) +
                         temp4 * rstate.yjm4(j) +
                         temp5 * rstate.fnt(j)) / wt, 2);
    }

    err = std::sqrt(err / static_cast<Real>(INT_NEQS));

    return err;
}


template <typename BurnT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real rockfrho (BurnT& state, rock4_t& rstate)
{


    // Rockfrho compute eigmax, a close upper bound of the spectral
    // radius of the Jacobian matrix using a power method
    // (J.N. Franklin (matrix theory).  The algorithm used is a small
    // change (initial vector and stopping criteria) of that of
    // Sommeijer-Shampine-Verwer, implemented in RKC.

    // note: in the original code, the following names were used:
    //   z = yjm1
    //  fz = yjm2

    constexpr int maxiter{50};
    constexpr Real safe{1.2_rt};

    const Real sqrtu = std::sqrt(UROUND);

    // The initial vectors for the power method are yn and
    // yn+c*f(v_n), where vn=f(yn) a perturbation of yn (if first
    // step) or a perturbation of the last computed eigenvector
    // (otherwise)

    if (rstate.nsteps == 0) {
        for (int i = 1; i <= INT_NEQS; ++i) {
            rstate.yjm2(i) = rstate.fn(i);
        }
        // call the RHS with yjm2 as input and yjm1 as output
        for (int i = 1; i <= INT_NEQS; ++i) {
            rstate.y_backup(i) = rstate.y(i);
            rstate.y(i) = rstate.yjm2(i);
        }
        rhs(rstate.t, state, rstate, rstate.yjm1);
        for (int i = 1; i <= INT_NEQS; ++i) {
            rstate.y(i) = rstate.y_backup(i);
        }
        rstate.nfesig++;
    } else {
        for (int i = 1; i <= INT_NEQS; ++i) {
            rstate.yjm1(i) = rstate.sprad(i);
        }
    }

    // Perturbation.

    Real ynor{};
    Real znor{};

    for (int i = 1; i <= INT_NEQS; ++i) {
        ynor += std::pow(rstate.yn(i), 2);
        znor += std::pow(rstate.yjm1(i), 2);
    }
    ynor = std::sqrt(ynor);
    znor = std::sqrt(znor);

    // Normalization of the vector z so that the difference z-yn lie
    // in a circle around yn (i.e has a constant modulus).

    Real dzyn{};

    if (ynor != 0.0_rt && znor != 0.0_rt) {
        dzyn = ynor * sqrtu;
        Real quot = dzyn / znor;
        for (int i = 1; i <= INT_NEQS; ++i) {
            rstate.yjm1(i) = rstate.yn(i) + rstate.yjm1(i) * quot;
        }
    } else if (ynor != 0.0_rt) {
        dzyn = ynor * sqrtu;
        for (int i = 1; i <= INT_NEQS; ++i) {
            rstate.yjm1(i) = rstate.yn(i) + rstate.yn(i) * sqrtu;
        }
    } else if (znor != 0.0_rt) {
        dzyn = UROUND;
        Real quot = dzyn / znor;
        for (int i = 1; i <= INT_NEQS; ++i) {
            rstate.yjm1(i) = rstate.yjm1(i) * quot;
        }
    } else {
        dzyn = UROUND;
        for (int i = 1; i <= INT_NEQS; ++i) {
            rstate.yjm1(i) = dzyn;
        }
    }

    // Start the power method.

    Real eigmax{};

    int ntest{};
    int nind{};
    int ind{};

    for (int iter = 0; iter < maxiter; ++iter) {
        // call RHS with yjm1 as input and store the output as yjm2

        for (int i = 1; i <= INT_NEQS; ++i) {
            rstate.y_backup(i) = rstate.y(i);
            rstate.y(i) = rstate.yjm1(i);
        }
        rhs(rstate.t, state, rstate, rstate.yjm2);
        for (int i = 1; i <= INT_NEQS; ++i) {
            rstate.y(i) = rstate.y_backup(i);
        }

        rstate.nfesig++;

        Real dfzfn{};
        for (int i = 1; i <= INT_NEQS; ++i) {
            dfzfn += std::pow(rstate.yjm2(i) - rstate.fn(i), 2);
        }
        dfzfn = std::sqrt(dfzfn);

        Real eigmaxo = eigmax;
        eigmax = dfzfn / dzyn;
        eigmax = safe * eigmax;

        // The stopping criteria is based on a relative error between
        // two succesive estimation ``eigmax'' of the spectral radius.

        if (iter >= 1 && std::abs(eigmax-eigmaxo) <= (eigmax*0.05_rt)) {
            // The last eigenvector is stored.
            for (int i = 1; i <= INT_NEQS; ++i) {
                rstate.sprad(i) = rstate.yjm1(i) - rstate.yn(i);
            }
            return eigmax;
        }

        // The next z is defined by z_new=yn+coef*(fz-fn) where coef
        // is choosen so that norm(z_new-yn)=norm(z_old-yn)

        if (dfzfn != 0.0_rt) {
            Real quot = dzyn / dfzfn;
            for (int i = 1; i <= INT_NEQS; ++i) {
                rstate.yjm1(i) = rstate.yn(i) + (rstate.yjm2(i) - rstate.fn(i)) * quot;
            }
        } else {
            // The new z is defined by an arbitrary perturbation of
            // the current approximation of the eigenvector.
            ntest = 0;
            nind = INT_NEQS;
            ind = 1 + iter+1 % nind;
            if (rstate.yjm1(ind) != rstate.yn(ind) || ntest == 10) {
                rstate.yjm1(ind) = rstate.yn(ind) - (rstate.yjm1(ind) - rstate.yn(ind));
            } else {
                nind = INT_NEQS + ind;
                ntest++;
            }
        }
    }

    // if we made it here, we failed, and the caller will abort
    rstate.idid = -3;
    return 0.0_rt;
}


template <typename BurnT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int rockfcore (BurnT& state, rock4_t& rstate)
{

    // Core integrator for ROCK4.

    // initializations

    rstate.nfe = 0;
    rstate.nsteps = 0;
    rstate.naccpt = 0;
    rstate.nreject = 0;
    rstate.nfesig = 0;
    rstate.max_stages = 0;
    rstate.max_rho = 0;
    rstate.min_rho = std::numeric_limits<decltype(rstate.min_rho)>::max();

    rstate.idid = 1;


    Real facmax{5.0_rt};
    Real told{0.0_rt};
    Real err{0.0_rt};

    int nrej{0};
    int nrho{0};
    int mdego{0};

    Real hp{rstate.h};

    bool last{false};
    bool reject{false};

    // initialization of the integration step.

    for (int i = 1; i <= INT_NEQS; ++i) {
        rstate.yn(i) = rstate.y(i);
    }
    rhs(rstate.t, state, rstate, rstate.fn);   // TODO: I think this can be removed -- we already filled for finding rstate.h

    rstate.nfe++;
    Real errp = err;
    Real eigmax{};

    while (true) {

        // Step-size is adjusted.

        if (1.01_rt * rstate.h >= std::abs(rstate.tout - rstate.t)) {
            rstate.h = std::abs(rstate.tout - rstate.t);
            last = true;
        }

#if 0
        if (rstate.h < 10.0_rt * UROUND) {
            rstate.idid = -2;
            return rstate.idid;
        }
#endif

        // Spectral radius.

        if (nrho == 0) {
            // Computed internally by rockfrho
            eigmax = rockfrho(state, rstate);
            if (rstate.idid == -3) {
                return rstate.idid;
            }
            rstate.max_rho = std::max(rstate.max_rho, static_cast<int>(std::round(eigmax)));
            rstate.min_rho = std::min(rstate.min_rho, static_cast<int>(std::round(eigmax)));
        }

        // The number of stages.

        rstate.mdeg = static_cast<int>(std::sqrt((3.0_rt + rstate.h * eigmax) / 0.353_rt)) + 1;

        if (rstate.mdeg > 152) {
            rstate.h = 0.8_rt * (std::pow(152.0_rt, 2) * 0.353_rt - 3.0_rt) / eigmax;
            rstate.mdeg = 152;
            last = false;
        }

        rstate.mdeg = std::max(rstate.mdeg, 5) - 4;

        if (rstate.mdeg != mdego) {
            mdegre(rstate);
        }

        if (rstate.mdeg + 4 > rstate.max_stages) {
            rstate.max_stages = rstate.mdeg + 4;
        }

        // Computation of an integration step.

        err = rfstep(state, rstate);

        mdego = rstate.mdeg;
        rstate.nsteps++;
        rstate.nfe += rstate.mdeg + 4;

        // Error control procedure.

        Real fac = std::pow(1.0_rt / err, 0.25_rt);

        if (errp != 0.0_rt && ! reject) {
            Real facp = std::pow(errp, 0.25_rt) * fac * fac * (rstate.h / hp);
            fac = std::min(fac, facp);
        }

        if (reject) {
            facmax = 1.0;
        }

        fac = std::min(facmax, std::max(0.1_rt, 0.8_rt * fac));
        Real hnew = rstate.h * fac;

        // Accepted step.
        if (err < 1.0_rt) {
            rstate.naccpt++;
            facmax = 2.0_rt;
            rstate.t += rstate.h;
            if (reject) {
                hnew = std::min(hnew, rstate.h);
                if (rstate.tout < rstate.t) {
                    hnew = std::max(hnew, rstate.h);
                }
                reject = false;
                nrej=0;
            }
            hp = rstate.h;
            rstate.h = hnew;
            nrho++;
            nrho = nrho+1 % 10;
            if (last) {
                return rstate.idid;
            } else {
                for (int i = 1; i <= INT_NEQS; ++i) {
                    rstate.yn(i) = rstate.y(i);
                    rstate.fn(i) = rstate.fnt(i);
                }
            }
            rstate.nfe++;
            errp = err;
            continue;
        }

        // Rejected step.
        rstate.nreject++;
        reject = true;
        last = false;
        rstate.h = 0.8_rt * hnew;
        if (rstate.nsteps == 0) {
            rstate.h = 0.1e0 * rstate.h;
        }
        if (told == rstate.t) {
            nrej++;
            if (nrej == 10) rstate.h = 1.0e-5;   // this may not be enough for us?
        }
        told = rstate.t;
        // The spectral radius is recomputed after a step failure
        if (nrho != 0) {
            nrho = 0;
        } else {
            nrho = 1;
        }
        continue;
    }

    return rstate.idid;
}


template <typename BurnT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int rock4 (BurnT& state, rock4_t& rstate)
{

    const Real rmax = 0.1_rt;
    const Real rmin = 10.0_rt * UROUND;

    // initializations

    // get the initial timestep

    rhs(rstate.t, state, rstate, rstate.fn);

    rstate.h = initial_dt(state, rstate, rstate.fn);

    std::cout << "initial timestep = " << rstate.h << std::endl;

    // Test the initial step size and tolerances.

#if 0
    if (rstate.h < 10.0_rt * UROUND) {
        rstate.idid = -1;
        return rstate.idid;
    }
#endif

    bool valid = true;

    if ((rtol_spec > rmax) || (rtol_enuc > rmax) ||
        (rtol_spec < rmin) || (rtol_spec < rmin)) {
        valid = false;
    }

    if (atol_spec < 0.0_rt || atol_enuc < 0.0_rt) {
        valid = false;
    }

    if (! valid) {
        rstate.idid = -1;  // TODO: remove idid from rstate
        return rstate.idid;
    }

    // call the integrator

    rockfcore(state, rstate);

    return rstate.idid;
}

#endif
