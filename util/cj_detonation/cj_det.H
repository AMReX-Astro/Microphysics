#ifndef CJ_DET_H
#define CJ_DET_H

AMREX_INLINE
void
adiabat(const eos_t eos_state_fuel, eos_t& eos_state_ash, const Real q, int& istatus) {


    const Real tol = 1.e-8_rt;
    const int max_iter = 50;

    // we want to zero e_1 + q - e_2 + 0.5*(p_1 + p_2)*(v_1 - v_2)
    // where v = 1/rho
    //
    // we have a rho_2 (the ash input), so we need to find the T that
    // makes this all work
    //
    // we assume that we come in with a reasonable guess for T

    bool converged = false;

    int iter = 0;

    while (!converged && iter < max_iter) {

        eos(eos_input_rt, eos_state_ash);

        Real f = eos_state_fuel.e + q - eos_state_ash.e +
            0.5_rt * (eos_state_fuel.p + eos_state_ash.p) *
            (1.0_rt / eos_state_fuel.rho - 1.0_rt / eos_state_ash.rho);

        Real dfdT = -eos_state_ash.dedT + 0.5_rt * eos_state_ash.dpdT *
            (1.0_rt / eos_state_fuel.rho - 1.0_rt / eos_state_ash.rho);

        Real dT = -f / dfdT;

        if (std::abs(dT) < tol * eos_state_ash.T) {
            converged = true;
        }

        eos_state_ash.T += dT;

        iter += 1;

    }

    if (!converged) {
        istatus = -1;
    } else {
        istatus = 0;
    }

}

AMREX_INLINE
void
cj_cond(const eos_t eos_state_fuel, eos_t& eos_state_ash, const Real q) {


    const Real tol = 1.e-8_rt;
    const int max_iter = 50;

    // iterate, picking the density that corresponds to the CJ point

    eos(eos_input_rt, eos_state_ash);

    Real drho = 1.e30_rt;

    // this is the density we find from the tangent point to the
    // Hugoniot curve

    eos_state_ash.rho = eos_state_fuel.rho *
        (1.0_rt + (eos_state_ash.p - eos_state_fuel.p) /
         (eos_state_ash.gam1 * eos_state_ash.p));

    int istatus;

    int iter = 0;
    bool converged = false;
    while (!converged && iter < max_iter) {

        Real rho_old = eos_state_ash.rho;

        adiabat(eos_state_fuel, eos_state_ash, q, istatus);

        // this is the density we find from the tangent point to the
        // Hugoniot curve

        eos_state_ash.rho = eos_state_fuel.rho *
            (1.0_rt + (eos_state_ash.p - eos_state_fuel.p) /
             (eos_state_ash.gam1 * eos_state_ash.p));

        Real drho = eos_state_ash.rho - rho_old;

        if (std::abs(drho) < tol * eos_state_ash.rho) {
            converged = true;
        }

       iter += 1;
    }


    if (!converged || istatus == -1) {
        amrex::Error("CJ did not converge");
    }
}

#endif
