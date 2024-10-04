#ifndef ROCK_H
#define ROCK_H

#include <rock_type.H>
#include <burn_type.H>
#ifdef STRANG
#include <integrator_type_strang.H>
#include <integrator_rhs_strang.H>
#endif
#ifdef SIMPLIFIED_SDC
#include <integrator_type_simplified_sdc.H>
#include <integrator_rhs_simplified_sdc.H>
#endif
#include <rock2_coefficients.H>
#include <rock4_coefficients.H>
#include <rkc_util.H>
#include <circle_theorem.H>
#include <integrator_data.H>

template <typename RockT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mdegre2(RockT& rstate)
{

    // this is the ROCK2 version

    constexpr amrex::Array1D<int, 1, 46> ms = {
        1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
        22, 24, 26, 28, 30, 33, 36, 39, 43, 47, 51, 56, 61, 66, 72, 78, 85, 93,
        102, 112, 123, 135, 148, 163, 180, 198};

    // Find the optimal degree.
    // MP(1): pointer which select the degree in ms(i)\1,2,..
    //        such that mdeg<=ms(i).
    // MP(2): pointer which gives the corresponding position
    //        of a_1 in the data recf for the selected degree.

    rstate.mp(2) = 1;

    // Find the degree

    for (int i = 1; i <= 46; ++i) {
        if ((ms(i) / rstate.mdeg) >= 1) {
            rstate.mdeg = ms(i);
            rstate.mp(1) = i;
            return;
        }
        rstate.mp(2) += ms(i)*2 - 1;
    }
}


template <typename RockT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mdegre4(RockT& rstate)
{

    // this is the ROCK4 version

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
}


template <typename BurnT, typename RockT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real r2step (BurnT& state, RockT& rstate)
{

    // Solut. at t+h by an explicit (mdeg+2)-stages formula.

    int mz = rstate.mp(1);
    int mr = rstate.mp(2);

    // First stage

    Real temp1 = rstate.h * rock2_coeffs::recf(mr);
    Real ci1 = rstate.t + temp1;
    Real ci2 = rstate.t + temp1;
    Real ci3 = rstate.t;

    Real temp2, temp3;

    for (int i = 1; i <= INT_NEQS; ++i) {
        rstate.yjm2(i) = rstate.yn(i);
        rstate.yjm1(i) = rstate.yn(i) + temp1 * rstate.fn(i);
        if (rstate.mdeg < 2) {
            rstate.y(i) = rstate.yjm1(i);
        }
    }

    // Stage for j=2..mdeg.

    for (int i = 2; i <= rstate.mdeg; ++i) {
        temp1 = rstate.h * rock2_coeffs::recf(mr+2*(i-2)+1);
        temp3 = -rock2_coeffs::recf(mr+2*(i-2)+2);
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

    // The finishing procedure

    temp1 = rstate.h * rock2_coeffs::fp1(mz);
    temp2 = rstate.h * rock2_coeffs::fp2(mz);
    rhs(ci1, state, rstate, rstate.yjm2);
    for (int j = 1; j <= INT_NEQS; ++j) {
        rstate.yjm1(j) = rstate.y(j) + temp1 * rstate.yjm2(j);
    }

    ci1 += temp1;

    // call the RHS with yjm1 as the input and store the output in y
    for (int i = 1; i <= INT_NEQS; ++i) {
        rstate.y(i) = rstate.yjm1(i);
    }
    rhs(ci1, state, rstate, rstate.y_backup);
    for (int i = 1; i <= INT_NEQS; ++i) {
        rstate.y(i) = rstate.y_backup(i);
    }

    Real err{};

    for (int j = 1; j <= INT_NEQS; ++j) {
        temp3 = temp2 * (rstate.y(j) - rstate.yjm2(j));
        rstate.y(j) = rstate.yjm1(j) + temp1 * rstate.y(j) + temp3;
        Real wt{};
        if (j <= NumSpec) {
            wt = rstate.rtol_spec * std::abs(rstate.y(j)) + rstate.atol_spec;
        } else {
            wt = rstate.rtol_enuc * std::abs(rstate.y(j)) + rstate.atol_enuc;
        }

        err += amrex::Math::powi<2>(temp3 / wt);
    }

    err = std::sqrt(err / static_cast<Real>(INT_NEQS));

    return err;
}

template <typename BurnT, typename RockT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real r4step (BurnT& state, RockT& rstate)
{

    // Solut. at t+h by an explicit (mdeg+4)-stages formula.

    int mz = rstate.mp(1);
    int mr = rstate.mp(2);

    // First stage

    Real temp1 = rstate.h * rock4_coeffs::recf(mr);
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
        temp1 = rstate.h * rock4_coeffs::recf(mr+2*(i-2)+1);
        temp3 = -rock4_coeffs::recf(mr+2*(i-2)+2);
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
    temp1 = rstate.h * rock4_coeffs::fpa(mz, 1);
    rhs(ci1, state, rstate, rstate.yjm1);
    for (int j = 1; j <= INT_NEQS; ++j) {
        rstate.yjm3(j) = rstate.y(j) + temp1 * rstate.yjm1(j);
    }

    // stage 2
    ci2 = ci1 + temp1;
    temp1 = rstate.h * rock4_coeffs::fpa(mz,2);
    temp2 = rstate.h * rock4_coeffs::fpa(mz,3);

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
    temp1 = rstate.h * rock4_coeffs::fpa(mz,4);
    temp2 = rstate.h * rock4_coeffs::fpa(mz,5);
    temp3 = rstate.h * rock4_coeffs::fpa(mz,6);

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
    temp1 = rstate.h * rock4_coeffs::fpb(mz,1);
    temp2 = rstate.h * rock4_coeffs::fpb(mz,2);
    temp3 = rstate.h * rock4_coeffs::fpb(mz,3);
    temp4 = rstate.h * rock4_coeffs::fpb(mz,4);

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
    temp1 = rstate.h * rock4_coeffs::fpbe(mz,1) - temp1;
    temp2 = rstate.h * rock4_coeffs::fpbe(mz,2) - temp2;
    temp3 = rstate.h * rock4_coeffs::fpbe(mz,3) - temp3;
    temp4 = rstate.h * rock4_coeffs::fpbe(mz,4) - temp4;
    temp5 = rstate.h * rock4_coeffs::fpbe(mz,5);

    rhs(rstate.t + rstate.h, state, rstate, rstate.fnt);


    Real err{};

    for (int j = 1; j <= INT_NEQS; ++j) {
        Real wt{};
        if (j <= NumSpec) {
            wt = rstate.rtol_spec * std::abs(rstate.y(j)) + rstate.atol_spec;
        } else {
            wt = rstate.rtol_enuc * std::abs(rstate.y(j)) + rstate.atol_enuc;
        }

        err += amrex::Math::powi<2>((temp1 * rstate.yjm1(j) +
                                     temp2 * rstate.yjm2(j) +
                                     temp3 * rstate.yjm3(j) +
                                     temp4 * rstate.yjm4(j) +
                                     temp5 * rstate.fnt(j)) / wt);
    }

    err = std::sqrt(err / static_cast<Real>(INT_NEQS));

    return err;
}


template <typename BurnT, typename RockT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int rockcore (BurnT& state, RockT& rstate)
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

    while (rstate.nsteps < integrator_rp::ode_max_steps) {

        // Step-size is adjusted.

        if (1.01_rt * rstate.h >= std::abs(rstate.tout - rstate.t)) {
            rstate.h = std::abs(rstate.tout - rstate.t);
            last = true;
        }

        // did we underflow?

        if (rstate.t + rstate.h == rstate.t) {
            return IERR_DT_UNDERFLOW;
        }

        // Spectral radius.

        if (nrho == 0) {
            if (integrator_rp::use_circle_theorem) {
                circle_theorem_sprad(rstate.t, state, rstate, eigmax);
            } else {
                int ierr = rkcrho(state, rstate, rstate.tout, eigmax);
                if (ierr != 0) {
                    return IERR_SPRAD_CONVERGENCE;
                }
            }
            rstate.max_rho = std::max(rstate.max_rho, static_cast<int>(std::round(eigmax)));
            rstate.min_rho = std::min(rstate.min_rho, static_cast<int>(std::round(eigmax)));
        }

        // The number of stages
        if (rock_order == 2) {
            rstate.mdeg = static_cast<int>(std::sqrt((1.5_rt + rstate.h * eigmax) / 0.811_rt)) + 1;

            if (rstate.mdeg > 200) {
                rstate.h = 0.8_rt * (amrex::Math::powi<2>(200.0_rt) * 0.811_rt - 1.5_rt) / eigmax;
                rstate.mdeg = 200;
                last = false;
            }

            rstate.mdeg = std::max(rstate.mdeg, 3) - 2;

            if (rstate.mdeg != mdego) {
                mdegre2(rstate);
            }

            if (rstate.mdeg + 2 > rstate.max_stages) {
                rstate.max_stages = rstate.mdeg + 2;
            }

        } else {
            rstate.mdeg = static_cast<int>(std::sqrt((3.0_rt + rstate.h * eigmax) / 0.353_rt)) + 1;

            if (rstate.mdeg > 152) {
                rstate.h = 0.8_rt * (amrex::Math::powi<2>(152.0_rt) * 0.353_rt - 3.0_rt) / eigmax;
                rstate.mdeg = 152;
                last = false;
            }

            rstate.mdeg = std::max(rstate.mdeg, 5) - 4;

            if (rstate.mdeg != mdego) {
                mdegre4(rstate);
            }

            if (rstate.mdeg + 4 > rstate.max_stages) {
                rstate.max_stages = rstate.mdeg + 4;
            }
        }

        // Computation of an integration step.

        if (rock_order == 2) {
            err = r2step(state, rstate);
            rstate.nfe += rstate.mdeg + 1;
        } else {
            err = r4step(state, rstate);
            rstate.nfe += rstate.mdeg + 4;
        }

        mdego = rstate.mdeg;
        rstate.nsteps++;

        // Error control procedure.

        Real fac;
        if (rock_order == 2) {
            fac = std::sqrt(1.0_rt / err);

            if (errp != 0.0_rt && ! reject) {
                Real facp = std::sqrt(errp) * fac * fac * (rstate.h / hp);
                fac = std::min(fac, facp);
            }

        } else {
            fac = std::pow(1.0_rt / err, 0.25_rt);

            if (errp != 0.0_rt && ! reject) {
                Real facp = std::pow(errp, 0.25_rt) * fac * fac * (rstate.h / hp);
                fac = std::min(fac, facp);
            }
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
            int rho_cycle = rock_order == 2 ? 25 : 10;
            nrho = nrho+1 % rho_cycle;
            if (last) {
                return IERR_SUCCESS;
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

    }

    return IERR_TOO_MANY_STEPS;
}


template <typename BurnT, typename RockT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int rock (BurnT& state, RockT& rstate)
{

    using namespace integrator_rp;

    const Real rmax = 0.1_rt;
    const Real rmin = 10.0_rt * UROUND;

    // initializations

    // get the initial timestep
    for (int i = 1; i <= INT_NEQS; ++i) {
        rstate.yn(i) = rstate.y(i);
    }

    rhs(rstate.t, state, rstate, rstate.fn);

    Real sprad{};
    if (integrator_rp::use_circle_theorem) {
        circle_theorem_sprad(rstate.t, state, rstate, sprad);
    } else {
        int ierr = rkcrho(state, rstate, rstate.tout, sprad);
        if (ierr != 0) {
            return IERR_SPRAD_CONVERGENCE;
        }
    }

    //rstate.h = initial_dt(state, rstate, rstate.fn);
    rstate.h = rkc_init_dt(state, rstate, rstate.tout, sprad);

    // Test the initial step size and tolerances.


    bool valid = true;

    if ((rtol_spec > rmax) || (rtol_enuc > rmax) ||
        (rtol_spec < rmin) || (rtol_enuc < rmin)) {
        valid = false;
    }

    if (atol_spec < 0.0_rt || atol_enuc < 0.0_rt) {
        valid = false;
    }

    if (! valid) {
        return IERR_BAD_INPUTS;
    }

    // call the integrator

    int ierr = rockcore(state, rstate);

    return ierr;
}

#endif
