// VODE integrator and Strang integrator glue for the Loom burn kernel
// (see README.md), for any rhs.H-framework network (NEW_NETWORK_IMPLEMENTATION),
// helmholtz, STRANG, no NSE, no NONAKA_PLOT, HIP branches: no
// ALLOW_JACOBIAN_CACHING (no jac_save), AMREX_USE_GPU (no printing, no
// xn_in in state_backup_t), no ALLOW_GPU_PRINTF.
//
// Translated:
//   integration/VODE/vode_type.H        dvode_t (POD of flat arrays, by pointer)
//   integration/VODE/vode_dvode.H       dvode
//   integration/VODE/vode_dvhin.H       dvhin
//   integration/VODE/vode_dvstep.H      advance_nordsieck, retract_nordsieck, dvstep
//   integration/VODE/vode_dvnlsd.H      dvnlsd
//   integration/VODE/vode_dvjac.H       dvjac (analytic + numerical branch)
//   integration/VODE/vode_dvjust.H      dvjust
//   integration/VODE/vode_dvset.H       dvset
//   integration/integrator_data.H       integrator_errors
//   integration/integrator_type.H       update_thermodynamics
//   integration/integrator_type_strang.H integrator_to_burn, burn_to_integrator, clean_state
//   integration/integrator_rhs_strang.H rhs, jac
//   integration/integrator_setup_strang.H state_backup_t, integrator_setup,
//                                       integrator_backup, integrator_cleanup
//   interfaces/burn_type.H              normalize_abundances_burn
//   integration/utils/numerical_jacobian.H numerical_jac (NOT reachable from
//                                       VODE: dvjac's jacobian_type == 2
//                                       branch differences rhs() inline;
//                                       ported on request, validated separately)
//
// Linear algebra: RHS::dgefa / RHS::dgesl (no pivoting) from network.H, so
// integrator.linalg_do_pivoting is not read and dvode_t has no pivot array.
//
// Control flow: every `return` inside a VODE loop became a flag + `break`
// (result code kept in a local), preserving each exit's side effects.
// Loop bodies that index arrays with runtime-dependent indices carry
// MPL_ASSUME bounds (no-op on the host).
#ifndef MPL_VODE_H
#define MPL_VODE_H

#include "support.H"
#include "burn.H"
#include "eos.H"
#include "network.H"

namespace mpl {

// ---------------------------------------------------------------------------
// runtime parameters (README rule 8)
// ---------------------------------------------------------------------------

// integrator_rp::* read on the burn path, plus network_rp::small_x
// (normalize_abundances_burn). bool parameters are stored as int (0/1).
struct IntegratorParams {
    Real integrator_rtol_spec;
    Real integrator_rtol_enuc;
    Real integrator_atol_spec;
    Real integrator_atol_enuc;
    Real integrator_retry_rtol_spec;
    Real integrator_retry_rtol_enuc;
    Real integrator_retry_atol_spec;
    Real integrator_retry_atol_enuc;
    int integrator_jacobian;
    int integrator_use_burn_retry;
    int integrator_retry_swap_jacobian;
    int integrator_scale_system;
    int integrator_call_eos_in_rhs;
    int integrator_integrate_energy;
    int integrator_use_number_densities;
    int integrator_subtract_internal_energy;
    int integrator_do_species_clip;
    int integrator_renormalize_abundances;
    int integrator_do_corrector_validation;
    int integrator_correct_jacobian_for_const_e;
    int integrator_ode_max_steps;
    Real integrator_ode_max_dt;
    Real integrator_SMALL_X_SAFE;
    Real integrator_MAX_TEMP;
    Real integrator_react_boost;
    Real integrator_species_failure_tolerance;
    Real integrator_X_reject_buffer;
    Real network_small_x;
};

// ---------------------------------------------------------------------------
// integration/integrator_data.H, integration/VODE/vode_type.H constants
// ---------------------------------------------------------------------------

// enum integrator_errors
constexpr int IERR_SUCCESS = 1;
constexpr int IERR_BAD_INPUTS = -1;
constexpr int IERR_DT_UNDERFLOW = -2;
constexpr int IERR_SPRAD_CONVERGENCE = -3;
constexpr int IERR_TOO_MANY_STEPS = -4;
constexpr int IERR_TOO_MUCH_ACCURACY_REQUESTED = -5;
constexpr int IERR_CORRECTOR_CONVERGENCE = -6;
constexpr int IERR_LU_DECOMPOSITION_ERROR = -7;
constexpr int IERR_BAD_STATE_IN_CORRECTOR = -8;

// std::numeric_limits<double>::epsilon()
constexpr Real UROUND = 2.220446049250313080847263336181640625e-16;
constexpr Real CCMXJ = 0.2e0;
constexpr Real HMIN = 0.0;
constexpr Real vode_increase_change_factor = 4.0;
constexpr Real vode_decrease_change_factor = 0.25;
constexpr Real vode_final_state_species_failure_tolerance_factor = 1.5;
constexpr int VODE_MAXORD = 5;
constexpr int VODE_LMAX = VODE_MAXORD + 1;

// network_properties.H aion[] / aion_inv[] are local const arrays where
// needed (README rule 7), filled from the generated MPL_NET_*_LIST macros.

// Nordsieck history yh(i, j), i = 1..neqs, j = 1..VODE_LMAX
#define MPL_YH(i, j) MPL_IDX2((i), (j), (mpl::VODE_LMAX + 1))
#define MPL_YH_SIZE ((mpl::neqs + 1) * (mpl::VODE_LMAX + 1))

// from integration/VODE/vode_type.H:38 dvode_t<int_neqs = 14>
// Array1D<.., 1, N> -> a[N + 1] (slot 0 unused); jac -> MPL_JAC (network.H);
// yh -> MPL_YH. No jac_save (caching compiled out), no pivot (RHS::dgefa).
struct dvode_t {
    Real CONP{};
    Real CRATE{};
    Real DRC{};
    Real ETA{};
    Real ETAMAX{};
    Real H{};
    Real HSCAL{};
    Real PRL1{};
    Real HMXI{};
    Real RC{};
    Real RL1{};
    Real tn{};

    int n_rhs{};
    int n_jac{};
    int n_step{};

    short ICF{};
    short IPUP{};
    short JCUR{};
    short L{};
    short NEWH{};
    short NEWQ{};
    short NQ{};
    short NQWAIT{};

    int NSLJ{};
    int NSLP{};

    short jacobian_type{};

    Real el[VODE_LMAX + 1]{};
    Real tau[VODE_LMAX + 1]{};
    Real tq[5 + 1]{};

    Real rtol_spec{};
    Real atol_spec{};
    Real rtol_enuc{};
    Real atol_enuc{};

    Real t{};
    Real tout{};

    Real y[neqs + 1]{};

    Real jac[MPL_JAC_SIZE]{};

    Real yh[MPL_YH_SIZE]{};

    Real ewt[neqs + 1]{};
    Real savf[neqs + 1]{};

    Real acor[neqs + 1]{};
};

// ---------------------------------------------------------------------------
// interfaces/burn_type.H
// ---------------------------------------------------------------------------

// from interfaces/burn_type.H:338 normalize_abundances_burn
MPL_DEVICE void normalize_abundances_burn(burn_t* state, const IntegratorParams* ip)
{
    Real sum = 0.0;
    for (int n = 0; n < NumSpec; ++n) {
        state->xn[n] = math::max(ip->network_small_x, math::min(1.0, state->xn[n]));
        sum += state->xn[n];
    }
    for (int n = 0; n < NumSpec; ++n) {
        state->xn[n] /= sum;
    }
}

// ---------------------------------------------------------------------------
// integration/integrator_type_strang.H, integration/integrator_type.H
// ---------------------------------------------------------------------------

// from integration/integrator_type_strang.H:10 integrator_to_burn
MPL_DEVICE void integrator_to_burn(const dvode_t* vs, burn_t* state, const IntegratorParams* ip)
{
    for (int n = 1; n <= NumSpec; ++n) {
        state->xn[n-1] = vs->y[n];
    }

    state->e = vs->y[net_ienuc];

    if (ip->integrator_scale_system) {
        state->e *= state->e_scale;
    }

    state->time = vs->t;
}

// from integration/integrator_type_strang.H:30 burn_to_integrator
MPL_DEVICE void burn_to_integrator(const burn_t* state, dvode_t* vs, const IntegratorParams* ip)
{
    for (int n = 1; n <= NumSpec; ++n) {
        vs->y[n] = state->xn[n-1];
    }

    vs->y[net_ienuc] = state->e;

    if (ip->integrator_scale_system) {
        vs->y[net_ienuc] /= state->e_scale;
    }
}

// from integration/integrator_type_strang.H:47 clean_state
MPL_DEVICE void clean_state(const Real time, burn_t* state, dvode_t* vs, const IntegratorParams* ip)
{
    (void) time;
    (void) state;

    if (ip->integrator_do_species_clip) {
        for (int n = 1; n <= NumSpec; ++n) {
            vs->y[n] = math::clamp(vs->y[n], ip->integrator_SMALL_X_SAFE, 1.0);
        }
    }

    if (ip->integrator_use_number_densities) {
        for (int n = 1; n <= NumSpec; ++n) {
            vs->y[n] = math::max(vs->y[n], ip->integrator_SMALL_X_SAFE);
        }
    }

    if (ip->integrator_renormalize_abundances) {
        Real sum = 0.0;

        for (int n = 1; n <= NumSpec; ++n) {
            sum += vs->y[n];
        }

        for (int n = 1; n <= NumSpec; ++n) {
            vs->y[n] /= sum;
        }
    }
}

// from integration/integrator_type.H:34 update_thermodynamics (no AUX_THERMO)
MPL_DEVICE void update_thermodynamics(burn_t* state, const dvode_t* vs, const IntegratorParams* ip,
                                      const EosParams* ep, EosTables tab)
{
    for (int n = 1; n <= NumSpec; ++n) {
        state->xn[n-1] = vs->y[n];
    }
    state->e = vs->y[net_ienuc];

    if (ip->integrator_scale_system) {
        state->e *= state->e_scale;
    }

    if (ip->integrator_call_eos_in_rhs) {
        eos(eos_input_re, state, ep, tab);
    }

    if (state->T_fixed > 0.0) {
        state->T = state->T_fixed;
    }
}

// ---------------------------------------------------------------------------
// integration/integrator_rhs_strang.H
// ---------------------------------------------------------------------------

// from integration/integrator_rhs_strang.H:24 rhs (ydot[1..neqs];
// in_jacobian only feeds NONAKA_PLOT and is dropped)
MPL_DEVICE void rhs(const Real time, burn_t* state, dvode_t* vs, Real* ydot,
                    const IntegratorParams* ip, const EosParams* ep, const NetParams* np, EosTables tab)
{
    const Real aion[NumSpec] = {MPL_NET_AION_LIST};

    clean_state(time, state, vs, ip);

    update_thermodynamics(state, vs, ip, ep, tab);

    if (state->T <= ep->EOSData_mintemp || state->T >= ip->integrator_MAX_TEMP) {

        for (int n = 1; n <= neqs; ++n) {
            ydot[n] = 0.0;
        }

        return;
    }

    state->time = time;

    actual_rhs(state, np, ydot);

    if (!ip->integrator_use_number_densities) {
        for (int n = 1; n <= NumSpec; ++n) {
            ydot[n] *= aion[n-1];
        }
    }

    if (ip->integrator_scale_system) {
        ydot[net_ienuc] /= state->e_scale;
    }

    if (!ip->integrator_integrate_energy) {
        ydot[net_ienuc] = 0.0;
    }

    if (ip->integrator_react_boost > 0.0) {
        for (int n = 1; n <= neqs; ++n) {
            ydot[n] *= ip->integrator_react_boost;
        }
    }
}

// from integration/integrator_rhs_strang.H:114 jac (pd[MPL_JAC(i, j)])
MPL_DEVICE void jac(const Real time, burn_t* state, dvode_t* vs, Real* pd,
                    const IntegratorParams* ip, const EosParams* ep, const NetParams* np, EosTables tab)
{
    (void) time;

    const Real aion[NumSpec] = {MPL_NET_AION_LIST};
    const Real aion_inv[NumSpec] = {MPL_NET_AION_INV_LIST};

    if (state->T <= ep->EOSData_mintemp || state->T >= ip->integrator_MAX_TEMP) {
        // pd.zero()
        for (int i = 1; i <= neqs; ++i) {
            for (int j = 1; j <= neqs; ++j) {
                pd[MPL_JAC(i, j)] = 0.0;
            }
        }
        return;
    }

    integrator_to_burn(vs, state, ip);

    actual_jac(state, np, pd);

    if (!ip->integrator_use_number_densities) {
        for (int j = 1; j <= NumSpec; ++j) {
            for (int i = 1; i <= neqs; ++i) {
                pd[MPL_JAC(j, i)] *= aion[j-1];
                pd[MPL_JAC(i, j)] *= aion_inv[j-1];
            }
        }
    }

    if (ip->integrator_correct_jacobian_for_const_e) {

        eos_re_extra_t eos_state{};
        eos_state.rho = state->rho;
        eos_state.T = state->T;
        eos_state.e = state->e;
        for (int n = 0; n < NumSpec; n++) {
            eos_state.xn[n] = state->xn[n];
        }

        eos(eos_input_re, &eos_state, ep, tab);

        eos_xderivs_t eos_xderivs{};
        composition_derivatives(&eos_state, &eos_xderivs);

        for (int jcol = 1; jcol <= NumSpec; ++jcol) {
            for (int irow = 1; irow <= neqs; ++irow) {
                pd[MPL_JAC(irow, jcol)] -= eos_xderivs.dedX[jcol-1] * pd[MPL_JAC(irow, net_ienuc)];
            }
        }
    }

    if (ip->integrator_scale_system) {
        for (int j = 1; j <= neqs; ++j) {
            pd[MPL_JAC(net_ienuc, j)] /= state->e_scale;
        }

        for (int i = 1; i <= neqs; ++i) {
            pd[MPL_JAC(i, net_ienuc)] *= state->e_scale;
        }
    }

    if (ip->integrator_react_boost > 0.0) {
        // pd.mul(react_boost)
        for (int i = 1; i <= neqs; ++i) {
            for (int j = 1; j <= neqs; ++j) {
                pd[MPL_JAC(i, j)] *= ip->integrator_react_boost;
            }
        }
    }

    if (!ip->integrator_integrate_energy) {
        for (int j = 1; j <= neqs; ++j) {
            pd[MPL_JAC(net_ienuc, j)] = 0.0;
        }
    }
}

// ---------------------------------------------------------------------------
// integration/utils/numerical_jacobian.H (not on the VODE path, see header)
// ---------------------------------------------------------------------------

// from integration/utils/numerical_jacobian.H:40 numerical_jac
// jac_info.h -> h; jac[MPL_JAC(i, j)].
MPL_DEVICE void numerical_jac(const burn_t* state, const Real h, Real* jac,
                              const IntegratorParams* ip, const EosParams* ep, const NetParams* np,
                              EosTables tab)
{
    const Real aion[NumSpec] = {MPL_NET_AION_LIST};

    Real ydotp[neqs + 1] = {};
    Real ydotm[neqs + 1] = {};

    burn_t state_delp = *state;

    actual_rhs(state, np, ydotm);

    for (int q = 1; q <= NumSpec; q++) {
        ydotm[q] *= aion[q-1];
    }

    Real fac = 0.0;
    Real w = 0.0;
    for (int i = 1; i <= NumSpec; i++) {
        w = 1.0 / (ip->integrator_rtol_spec * math::abs(state->xn[i-1]) + ip->integrator_atol_spec);
        fac += (ydotm[i] * w) * (ydotm[i] * w);
    }

    w = 1.0 / (ip->integrator_rtol_enuc * math::abs(state->e) + ip->integrator_atol_enuc);
    fac += (ydotm[net_ienuc] * w) * (ydotm[net_ienuc] * w);
    fac = math::sqrt(fac / neqs);

    Real r0 = 1000.0 * math::abs(h) * UROUND * neqs * fac;
    if (r0 == 0.0) {
        r0 = 1.0;
    }

    for (int n = 1; n <= NumSpec; n++) {
        Real yj = state_delp.xn[n-1];

        w = ip->integrator_rtol_spec * math::abs(yj) + ip->integrator_atol_spec;

        Real dy = math::max(math::sqrt(UROUND) * math::abs(yj), r0 * w);

        state_delp.xn[n-1] += dy;

        actual_rhs(&state_delp, np, ydotp);

        for (int q = 1; q <= NumSpec; q++) {
            ydotp[q] *= aion[q-1];
        }

        for (int m = 1; m <= neqs; m++) {
            jac[MPL_JAC(m, n)] = (ydotp[m] - ydotm[m]) / dy;
        }

        state_delp.xn[n-1] = yj;
    }

    w = ip->integrator_rtol_enuc * math::abs(state->T) + ip->integrator_atol_enuc;

    Real dy = math::max(math::sqrt(UROUND) * math::abs(state->T), r0 * w);

    state_delp.T += dy;

    if (state_delp.T <= ep->EOSData_mintemp || state_delp.T >= ip->integrator_MAX_TEMP) {
        // jac.zero()
        for (int i = 1; i <= neqs; ++i) {
            for (int j = 1; j <= neqs; ++j) {
                jac[MPL_JAC(i, j)] = 0.0;
            }
        }
        return;
    }

    eos(eos_input_rt, &state_delp, ep, tab);

    actual_rhs(&state_delp, np, ydotp);

    for (int q = 1; q <= NumSpec; q++) {
        ydotp[q] *= aion[q-1];
    }

    for (int m = 1; m <= neqs; m++) {
        jac[MPL_JAC(m, net_ienuc)] = (ydotp[m] - ydotm[m]) / dy;
    }

    eos_re_extra_t eos_state{};
    eos_state.rho = state->rho;
    eos_state.T = state->T;
    eos_state.e = state->e;
    for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = state->xn[n];
    }

    eos(eos_input_re, &eos_state, ep, tab);

    eos_xderivs_t eos_xderivs{};
    composition_derivatives(&eos_state, &eos_xderivs);

    for (int m = 1; m <= neqs; m++) {
        jac[MPL_JAC(m, net_ienuc)] *= 1.0 / eos_state.cv;
    }

    for (int n = 1; n <= NumSpec; n++) {
        for (int m = 1; m <= neqs; m++) {
            jac[MPL_JAC(m, n)] -= eos_xderivs.dedX[n-1] * jac[MPL_JAC(m, net_ienuc)];
        }
    }

    if (ip->integrator_scale_system) {
        for (int n = 1; n <= neqs; ++n) {
            jac[MPL_JAC(net_ienuc, n)] /= state->e_scale;
        }

        for (int m = 1; m <= neqs; ++m) {
            jac[MPL_JAC(m, net_ienuc)] *= state->e_scale;
        }
    }

    if (ip->integrator_react_boost > 0.0) {
        for (int i = 1; i <= neqs; ++i) {
            for (int j = 1; j <= neqs; ++j) {
                jac[MPL_JAC(i, j)] *= ip->integrator_react_boost;
            }
        }
    }

    if (!ip->integrator_integrate_energy) {
        for (int n = 1; n <= neqs; ++n) {
            jac[MPL_JAC(net_ienuc, n)] = 0.0;
        }
    }
}

// ---------------------------------------------------------------------------
// integration/VODE
// ---------------------------------------------------------------------------

// from integration/VODE/vode_dvset.H:7 dvset
MPL_DEVICE void dvset(dvode_t* vs)
{
    constexpr Real CORTES = 0.1e0;

    const Real FLOTL = vs->L;
    const int NQM1 = vs->NQ - 1;
    const int NQM2 = vs->NQ - 2;

    for (int i = 3; i <= vs->L; ++i) {
        MPL_ASSUME(3 <= i && i <= VODE_LMAX);
        vs->el[i] = 0.0;
    }

    vs->el[1] = 1.0;
    vs->el[2] = 1.0;
    Real ALPH0 = -1.0;
    Real AHATN0 = -1.0;
    Real HSUM = vs->H;
    Real RXI = 1.0;
    Real RXIS = 1.0;

    if (vs->NQ != 1) {
        for (int j = 1; j <= NQM2; ++j) {
            MPL_ASSUME(1 <= j && j <= VODE_MAXORD - 2);
            HSUM += vs->tau[j];
            RXI = vs->H / HSUM;
            ALPH0 -= 1.0 / (j+1);
            for (int iback = 1; iback <= j+1; ++iback) {
                const int i = (j + 3) - iback;
                MPL_ASSUME(2 <= i && i <= VODE_LMAX);
                vs->el[i] += vs->el[i-1] * RXI;
            }
        }

        ALPH0 -= 1.0 / vs->NQ;
        RXIS = -vs->el[2] - ALPH0;
        MPL_ASSUME(1 <= NQM1 && NQM1 <= VODE_LMAX);
        HSUM += vs->tau[NQM1];
        RXI = vs->H / HSUM;
        AHATN0 = -vs->el[2] - RXI;
        for (int iback = 1; iback <= vs->NQ; ++iback) {
            const int i = (vs->NQ + 2) - iback;
            MPL_ASSUME(2 <= i && i <= VODE_LMAX);
            vs->el[i] += vs->el[i-1] * RXIS;
        }
    }

    // Loom's __builtin_assume takes plain values, not member paths: copy first
    const int La = vs->L;
    MPL_ASSUME(1 <= La && La <= VODE_LMAX);
    const int NQa = vs->NQ;
    MPL_ASSUME(1 <= NQa && NQa <= VODE_MAXORD);

    const Real T1 = 1.0 - AHATN0 + ALPH0;
    const Real T2 = 1.0 + vs->NQ * T1;
    vs->tq[2] = math::abs(ALPH0 * T2 / T1);
    vs->tq[5] = math::abs(T2 / (vs->el[La] * RXI / RXIS));

    if (vs->NQWAIT == 1) {
        const Real CNQM1 = RXIS / vs->el[La];
        const Real T3 = ALPH0 + 1.0 / vs->NQ;
        const Real T4 = AHATN0 + RXI;
        Real ELP = T3 / (1.0 - T4 + T3);
        vs->tq[1] = math::abs(ELP / CNQM1);
        HSUM += vs->tau[NQa];
        RXI = vs->H / HSUM;
        const Real T5 = ALPH0 - 1.0 / (vs->NQ+1);
        const Real T6 = AHATN0 - RXI;
        ELP = T2 / (1.0 - T6 + T5);
        vs->tq[3] = math::abs(ELP * RXI * (FLOTL + 1.0) * T5);
    }

    vs->tq[4] = CORTES * vs->tq[2];
}

// from integration/VODE/vode_dvjust.H:7 dvjust
MPL_DEVICE void dvjust(int IORD, burn_t* state, dvode_t* vs)
{
    (void) state;

    if ((vs->NQ == 2) && (IORD != 1)) {
        return;
    }

    const int NQM1 = vs->NQ - 1;
    const int NQM2 = vs->NQ - 2;

    Real XI = 0.0;

    if (IORD != 1) {

        // Order decrease.
        for (int j = 1; j <= VODE_LMAX; ++j) {
            vs->el[j] = 0.0;
        }
        vs->el[3] = 1.0;
        Real HSUM = 0.0;
        for (int j = 1; j <= NQM2; ++j) {
            MPL_ASSUME(1 <= j && j <= VODE_MAXORD - 2);
            HSUM += vs->tau[j];
            XI = HSUM / vs->HSCAL;
            for (int iback = 1; iback <= j+1; ++iback) {
                const int i = (j + 4) - iback;
                MPL_ASSUME(2 <= i && i <= VODE_LMAX);
                vs->el[i] = vs->el[i] * XI + vs->el[i-1];
            }
        }

        const int Lb = vs->L;
        MPL_ASSUME(1 <= Lb && Lb <= VODE_LMAX);
        for (int j = 3; j <= vs->NQ; ++j) {
            MPL_ASSUME(3 <= j && j <= VODE_LMAX);
            for (int i = 1; i <= neqs; ++i) {
                vs->yh[MPL_YH(i, j)] -= vs->yh[MPL_YH(i, Lb)] * vs->el[j];
            }
        }

    }
    else {

        // Order increase.
        for (int j = 1; j <= VODE_LMAX; ++j) {
            vs->el[j] = 0.0;
        }

        vs->el[3] = 1.0;
        Real ALPH0 = -1.0;
        Real ALPH1 = 1.0;
        Real PROD = 1.0;
        Real XIOLD = 1.0;
        Real HSUM = vs->HSCAL;

        if (vs->NQ != 1) {
            for (int j = 1; j <= NQM1; ++j) {
                MPL_ASSUME(1 <= j && j <= VODE_MAXORD - 1);
                HSUM += vs->tau[j+1];
                XI = HSUM / vs->HSCAL;
                PROD *= XI;
                ALPH0 -= 1.0 / (j + 1);
                ALPH1 += 1.0 / XI;
                for (int iback = 1; iback <= j+1; ++iback) {
                    int i = (j + 4) - iback;
                    MPL_ASSUME(2 <= i && i <= VODE_LMAX);
                    vs->el[i] = vs->el[i] * XIOLD + vs->el[i-1];
                }
                XIOLD = XI;
            }
        }

        Real T1 = (-ALPH0 - ALPH1) / PROD;
        const int Lc = vs->L;
        MPL_ASSUME(1 <= Lc && Lc <= VODE_LMAX - 1);
        for (int i = 1; i <= neqs; ++i) {
            vs->yh[MPL_YH(i, Lc+1)] = T1 * vs->yh[MPL_YH(i, VODE_LMAX)];
        }

        for (int j = 3; j <= vs->NQ + 1; ++j) {
            MPL_ASSUME(3 <= j && j <= VODE_LMAX);
            for (int i = 1; i <= neqs; ++i) {
                vs->yh[MPL_YH(i, j)] += vs->el[j] * vs->yh[MPL_YH(i, Lc+1)];
            }
        }
    }
}

// from integration/VODE/vode_dvstep.H:10 advance_nordsieck
MPL_DEVICE void advance_nordsieck(burn_t* state, dvode_t* vs)
{
    (void) state;

    for (int k = vs->NQ; k >= 1; --k) {
        MPL_ASSUME(1 <= k && k <= VODE_MAXORD);
        for (int j = k; j <= vs->NQ; ++j) {
            MPL_ASSUME(1 <= j && j <= VODE_LMAX - 1);
            for (int i = 1; i <= neqs; ++i) {
                vs->yh[MPL_YH(i, j)] += vs->yh[MPL_YH(i, j+1)];
            }
        }
    }
}

// from integration/VODE/vode_dvstep.H:31 retract_nordsieck
MPL_DEVICE void retract_nordsieck(burn_t* state, dvode_t* vs)
{
    (void) state;

    for (int k = vs->NQ; k >= 1; --k) {
        MPL_ASSUME(1 <= k && k <= VODE_MAXORD);
        for (int j = k; j <= vs->NQ; ++j) {
            MPL_ASSUME(1 <= j && j <= VODE_LMAX - 1);
            for (int i = 1; i <= neqs; ++i) {
                vs->yh[MPL_YH(i, j)] -= vs->yh[MPL_YH(i, j+1)];
            }
        }
    }
}

// rescale yh(:, 2..L) by successive powers of ETA (the block repeated in
// dvstep at five places; same statements, same order)
MPL_DEVICE void dvstep_rescale_history(dvode_t* vs)
{
    Real R = 1.0;

    for (int j = 2; j <= vs->L; ++j) {
        MPL_ASSUME(2 <= j && j <= VODE_LMAX);
        R *= vs->ETA;
        for (int i = 1; i <= neqs; ++i) {
            vs->yh[MPL_YH(i, j)] *= R;
        }
    }
}

// from integration/VODE/vode_dvjac.H:16 dvjac (no ALLOW_JACOBIAN_CACHING;
// NEW_NETWORK_IMPLEMENTATION: RHS::dgefa)
MPL_DEVICE void dvjac(int& IERPJ, burn_t* state, dvode_t* vs,
                      const IntegratorParams* ip, const EosParams* ep, const NetParams* np, EosTables tab)
{
    IERPJ = 0;

    if (vs->jacobian_type == 1) {

        vs->n_jac += 1;

        vs->NSLJ = vs->n_step;

        vs->JCUR = 1;

        jac(vs->tn, state, vs, vs->jac, ip, ep, np, tab);

    }
    else {

        vs->n_jac += 1;

        vs->NSLJ = vs->n_step;

        vs->JCUR = 1;

        Real fac = 0.0;
        for (int i = 1; i <= neqs; ++i) {
            fac += (vs->savf[i] * vs->ewt[i]) * (vs->savf[i] * vs->ewt[i]);
        }
        fac = math::sqrt(fac / neqs);

        Real R0 = 1000.0 * math::abs(vs->H) * UROUND * neqs * fac;
        if (R0 == 0.0) {
            R0 = 1.0;
        }

        for (int j = 1; j <= neqs; ++j) {
            const Real yj = vs->y[j];

            const Real R = math::max(math::sqrt(UROUND) * math::abs(yj), R0 / vs->ewt[j]);
            vs->y[j] += R;
            fac = 1.0 / R;

            rhs(vs->tn, state, vs, vs->acor, ip, ep, np, tab);
            for (int i = 1; i <= neqs; ++i) {
                vs->jac[MPL_JAC(i, j)] = (vs->acor[i] - vs->savf[i]) * fac;
            }

            vs->y[j] = yj;
        }

        vs->n_rhs += neqs;
    }

    const Real hrl1 = vs->H * vs->RL1;
    const Real con = -hrl1;

    // jac.mul(con); jac.add_identity()
    for (int i = 1; i <= neqs; ++i) {
        for (int j = 1; j <= neqs; ++j) {
            vs->jac[MPL_JAC(i, j)] *= con;
        }
    }
    for (int i = 1; i <= neqs; ++i) {
        vs->jac[MPL_JAC(i, i)] += 1.0;
    }

    int IER = rhs_dgefa(vs->jac);

    if (IER != 0) {
        IERPJ = 1;
    }
}

// from integration/VODE/vode_dvnlsd.H:11 dvnlsd
// The three `return ACNRM;` exits inside the loops set `failed` and break.
MPL_DEVICE Real dvnlsd(int& NFLAG, burn_t* state, dvode_t* vs,
                       const IntegratorParams* ip, const EosParams* ep, const NetParams* np, EosTables tab)
{
    Real ACNRM = 1.e10;

    constexpr Real CCMAX = 0.3e0;
    constexpr Real CRDOWN = 0.3e0;
    constexpr Real RDIV = 2.0e0;
    constexpr int MAXCOR = 3;
    constexpr int MSBP = 20;

    Real DEL = 0.0;
    int M = 0;

    if (NFLAG == 0) {
        vs->ICF = 0;
    }
    if (NFLAG == -2) {
        vs->IPUP = 1;
    }

    vs->DRC = math::abs(vs->RC - 1.0);
    if (vs->DRC > CCMAX || vs->n_step >= vs->NSLP + MSBP) {
        vs->IPUP = 1;
    }

    bool converged = false;
    bool failed = false;

    while (true) {

        M = 0;
        Real DELP = 0.0;

        for (int i = 1; i <= neqs; ++i) {
            vs->y[i] = vs->yh[MPL_YH(i, 1)];
        }

        rhs(vs->tn, state, vs, vs->savf, ip, ep, np, tab);
        vs->n_rhs += 1;

        if (vs->IPUP == 1) {

            int IERPJ = 0;
            dvjac(IERPJ, state, vs, ip, ep, np, tab);

            vs->IPUP = 0;
            vs->RC = 1.0;
            vs->DRC = 0.0;
            vs->CRATE = 1.0;
            vs->NSLP = vs->n_step;

            if (IERPJ != 0) {
                NFLAG = -1;
                vs->ICF = 2;
                vs->IPUP = 1;
                failed = true;
                break;  // stock: return ACNRM
            }
        }

        for (int i = 1; i <= neqs; ++i) {
            vs->acor[i] = 0.0;
        }

        // Corrector iteration loop.
        while (true) {

            for (int i = 1; i <= neqs; ++i) {
                vs->y[i] = (vs->RL1 * vs->H) * vs->savf[i] -
                           (vs->RL1 * vs->yh[MPL_YH(i, 2)] + vs->acor[i]);
            }

            rhs_dgesl(vs->jac, vs->y);

            if (vs->RC != 1.0) {
                const Real CSCALE = 2.0 / (1.0 + vs->RC);
                for (int i = 1; i <= neqs; ++i) {
                    vs->y[i] *= CSCALE;
                }
            }

            DEL = 0.0;
            for (int i = 1; i <= neqs; ++i) {
                DEL += (vs->y[i] * vs->ewt[i]) * (vs->y[i] * vs->ewt[i]);
            }
            DEL = math::sqrt(DEL / neqs);

            for (int i = 1; i <= neqs; ++i) {
                vs->acor[i] += vs->y[i];
            }

            for (int i = 1; i <= neqs; ++i) {
                vs->y[i] = vs->yh[MPL_YH(i, 1)] + vs->acor[i];
            }

            if (ip->integrator_do_corrector_validation && !ip->integrator_use_number_densities) {
                const Real thresh = ip->integrator_species_failure_tolerance;
                bool fail = false;
                for (int i = 1; i <= NumSpec; ++i) {
                    // loop with break: Loom can't bound i from the loop alone
                    MPL_ASSUME(1 <= i && i <= NumSpec);
                    if (vs->y[i] < -thresh) {
                        fail = true;
                        break;
                    }
                }

                if (fail) {
                    NFLAG = -1;
                    vs->ICF = 2;
                    vs->IPUP = 1;
                    failed = true;
                    break;  // stock: return ACNRM
                }
            }

            if (M != 0) {
                vs->CRATE = math::max(CRDOWN * vs->CRATE, DEL / DELP);
            }

            const Real DCON = DEL * math::min(1.0, vs->CRATE) / vs->tq[4];
            if (DCON <= 1.0) {
                converged = true;
                break;
            }

            M += 1;
            if (M == MAXCOR) {
                break;
            }

            if (M >= 2 && DEL > RDIV * DELP) {
                break;
            }

            DELP = DEL;
            rhs(vs->tn, state, vs, vs->savf, ip, ep, np, tab);
            vs->n_rhs += 1;
        }

        if (failed) {
            break;
        }

        if (converged) {
            break;
        }

        if (vs->JCUR == 1) {
            NFLAG = -1;
            vs->ICF = 2;
            vs->IPUP = 1;
            failed = true;
            break;  // stock: return ACNRM
        }

        vs->ICF = 1;
        vs->IPUP = 1;
    }

    if (failed) {
        return ACNRM;
    }

    NFLAG = 0;
    vs->JCUR = 0;
    vs->ICF = 0;
    if (M == 0) {
        ACNRM = DEL;
    }

    if (M > 0) {
        ACNRM = 0.0;
        for (int i = 1; i <= neqs; ++i) {
            ACNRM += (vs->acor[i] * vs->ewt[i]) * (vs->acor[i] * vs->ewt[i]);
        }
        ACNRM = math::sqrt(ACNRM / neqs);
    }

    return ACNRM;
}

// from integration/VODE/vode_dvstep.H:52 dvstep (STRANG constraint block)
// The four `return kflag;` exits inside the step loop set `finished` and break.
MPL_DEVICE int dvstep(burn_t* state, dvode_t* vs,
                      const IntegratorParams* ip, const EosParams* ep, const NetParams* np, EosTables tab)
{
    constexpr int KFC = -3;
    constexpr int KFH = -7;
    constexpr int MXNCF = 10;
    constexpr Real ADDON = 1.0e-6;
    constexpr Real BIAS1 = 6.0e0;
    constexpr Real BIAS2 = 6.0e0;
    constexpr Real BIAS3 = 10.0e0;
    constexpr Real ETACF = 0.25e0;
    constexpr Real ETAMIN = 0.1e0;
    constexpr Real ETAMXF = 0.2e0;
    constexpr Real ETAMX2 = 10.0e0;
    constexpr Real ETAMX3 = 10.0e0;
    constexpr Real ONEPSM = 1.00001e0;
    constexpr Real THRESH = 1.5e0;

    Real DSM = 0.0;
    int NCF = 0;
    int NFLAG = 0;

    int kflag = 0;
    Real TOLD = vs->tn;
    NCF = 0;
    vs->JCUR = 0;
    NFLAG = 0;

    if (vs->NEWH != 0) {

        if (vs->NEWQ < vs->NQ) {
            dvjust(-1, state, vs);
            vs->NQ = vs->NEWQ;
            vs->L = static_cast<short>(vs->NQ + 1);
            vs->NQWAIT = vs->L;
        }
        else if (vs->NEWQ > vs->NQ) {
            dvjust(1, state, vs);
            vs->NQ = vs->NEWQ;
            vs->L = static_cast<short>(vs->NQ + 1);
            vs->NQWAIT = vs->L;
        }

        dvstep_rescale_history(vs);

        vs->H = vs->HSCAL * vs->ETA;
        vs->HSCAL = vs->H;
        vs->RC = vs->RC * vs->ETA;
    }

    Real y_save[neqs + 1] = {};
    for (int i = 1; i <= neqs; ++i) {
        y_save[i] = vs->y[i];
    }

    bool finished = false;

    while (true) {

        vs->tn += vs->H;

        advance_nordsieck(state, vs);

        dvset(vs);

        vs->RL1 = 1.0 / vs->el[2];
        vs->RC *= (vs->RL1 / vs->PRL1);
        vs->PRL1 = vs->RL1;

        const Real ACNRM = dvnlsd(NFLAG, state, vs, ip, ep, np, tab);

        if (NFLAG != 0) {

            NCF += 1;
            vs->ETAMAX = 1.0;
            vs->tn = TOLD;

            retract_nordsieck(state, vs);

            if (math::abs(vs->H) <= HMIN * ONEPSM) {
                kflag = -2;
                finished = true;
                break;  // stock: return kflag
            }

            if (NCF == MXNCF) {
                kflag = -2;
                finished = true;
                break;  // stock: return kflag
            }

            vs->ETA = ETACF;
            vs->ETA = math::max(vs->ETA, HMIN / math::abs(vs->H));

            dvstep_rescale_history(vs);

            vs->H = vs->HSCAL * vs->ETA;
            vs->HSCAL = vs->H;
            vs->RC *= vs->ETA;

            continue;
        }

        bool valid_update = true;

        for (int i = 1; i <= NumSpec; ++i) {
            // loop with break: Loom can't bound i from the loop alone
            MPL_ASSUME(1 <= i && i <= NumSpec);

            if (math::abs(y_save[i]) > ip->integrator_X_reject_buffer * vs->atol_spec &&
                math::abs(vs->y[i]) > ip->integrator_X_reject_buffer * vs->atol_spec &&
                (math::abs(vs->y[i]) > vode_increase_change_factor * math::abs(y_save[i]) ||
                 math::abs(vs->y[i]) < vode_decrease_change_factor * math::abs(y_save[i]))) {
                valid_update = false;
                break;
            }

            if (vs->y[i] < -ip->integrator_species_failure_tolerance) {
                valid_update = false;
                break;
            }

            if (!ip->integrator_use_number_densities) {
                if (vs->y[i] > 1.0 + ip->integrator_species_failure_tolerance) {
                    valid_update = false;
                    break;
                }
            }
        }

        DSM = ACNRM / vs->tq[2];
        if (DSM <= 1.0 && valid_update) {

            kflag = 0;
            vs->n_step += 1;
            for (int iback = 1; iback <= vs->NQ; ++iback) {
                const int i = vs->L - iback;
                MPL_ASSUME(1 <= i && i <= VODE_LMAX - 1);
                vs->tau[i+1] = vs->tau[i];
            }

            vs->tau[1] = vs->H;
            for (int j = 1; j <= vs->L; ++j) {
                MPL_ASSUME(1 <= j && j <= VODE_LMAX);
                for (int i = 1; i <= neqs; ++i) {
                    vs->yh[MPL_YH(i, j)] += vs->el[j] * vs->acor[i];
                }
            }

            vs->NQWAIT -= 1;
            if ((vs->L != VODE_LMAX) && (vs->NQWAIT == 1)) {
                for (int i = 1; i <= neqs; ++i) {
                    vs->yh[MPL_YH(i, VODE_LMAX)] = vs->acor[i];
                }
                vs->CONP = vs->tq[5];
            }

            if (vs->ETAMAX != 1.0) {
                break;  // on to the order/step selection below
            }

            if (vs->NQWAIT < 2) {
                vs->NQWAIT = 2;
            }

            vs->NEWQ = vs->NQ;
            vs->NEWH = 0;
            vs->ETA = 1.0;
            vs->ETAMAX = ETAMX3;
            if (vs->n_step <= 10) {
                vs->ETAMAX = ETAMX2;
            }

            Real R = 1.0 / vs->tq[2];
            for (int i = 1; i <= neqs; ++i) {
                vs->acor[i] *= R;
            }

            finished = true;
            break;  // stock: return kflag
        }

        kflag -= 1;
        NFLAG = -2;
        vs->tn = TOLD;

        retract_nordsieck(state, vs);

        if (math::abs(vs->H) <= HMIN * ONEPSM) {
            kflag = -1;
            finished = true;
            break;  // stock: return kflag
        }

        vs->ETAMAX = 1.0;

        if (kflag > KFC) {

            Real FLOTL = vs->L;
            vs->ETA = 1.0 / (math::pow(BIAS2 * DSM, 1.0 / FLOTL) + ADDON);
            vs->ETA = math::max(math::max(vs->ETA, HMIN / math::abs(vs->H)), ETAMIN);
            if ((kflag <= -2) && (vs->ETA > ETAMXF)) {
                vs->ETA = ETAMXF;
            }

            dvstep_rescale_history(vs);

            vs->H = vs->HSCAL * vs->ETA;
            vs->HSCAL = vs->H;
            vs->RC = vs->RC * vs->ETA;

            continue;
        }

        if (kflag == KFH) {
            kflag = -1;
            finished = true;
            break;  // stock: return kflag
        }

        if (vs->NQ != 1) {
            vs->ETA = math::max(ETAMIN, HMIN / math::abs(vs->H));
            dvjust(-1, state, vs);
            vs->L = vs->NQ;
            vs->NQ -= 1;
            vs->NQWAIT = vs->L;

            dvstep_rescale_history(vs);

            vs->H = vs->HSCAL * vs->ETA;
            vs->HSCAL = vs->H;
            vs->RC = vs->RC * vs->ETA;

            continue;
        }

        vs->ETA = math::max(ETAMIN, HMIN / math::abs(vs->H));
        vs->H *= vs->ETA;
        vs->HSCAL = vs->H;
        vs->tau[1] = vs->H;
        rhs(vs->tn, state, vs, vs->savf, ip, ep, np, tab);
        vs->n_rhs += 1;
        for (int i = 1; i <= neqs; ++i) {
            vs->yh[MPL_YH(i, 2)] = vs->H * vs->savf[i];
        }

        vs->NQWAIT = 10;
    }

    if (finished) {
        return kflag;
    }

    bool already_set_eta = false;

    Real FLOTL = vs->L;
    const Real ETAQ = 1.0 / (math::pow(BIAS2 * DSM, 1.0 / FLOTL) + ADDON);
    if (vs->NQWAIT == 0) {
        vs->NQWAIT = 2;

        Real ETAQM1 = 0.0;

        if (vs->NQ != 1) {
            Real DDN = 0.0;
            const int Le = vs->L;
            MPL_ASSUME(1 <= Le && Le <= VODE_LMAX);
            for (int i = 1; i <= neqs; ++i) {
                DDN += (vs->yh[MPL_YH(i, Le)] * vs->ewt[i]) * (vs->yh[MPL_YH(i, Le)] * vs->ewt[i]);
            }
            DDN = math::sqrt(DDN / neqs) / vs->tq[1];
            ETAQM1 = 1.0 / (math::pow(BIAS1 * DDN, 1.0 / (FLOTL - 1.0)) + ADDON);
        }

        Real ETAQP1 = 0.0;

        if (vs->L != VODE_LMAX) {
            // std::pow(Real, short L): the short promotes to double
            Real CNQUOT = (vs->tq[5] / vs->CONP) * math::pow(vs->H / vs->tau[2], static_cast<Real>(vs->L));
            for (int i = 1; i <= neqs; ++i) {
                vs->savf[i] = vs->acor[i] - CNQUOT * vs->yh[MPL_YH(i, VODE_LMAX)];
            }
            Real DUP = 0.0;
            for (int i = 1; i <= neqs; ++i) {
                DUP += (vs->savf[i] * vs->ewt[i]) * (vs->savf[i] * vs->ewt[i]);
            }
            DUP = math::sqrt(DUP / neqs) / vs->tq[3];
            ETAQP1 = 1.0 / (math::pow(BIAS3 * DUP, 1.0 / (FLOTL + 1.0)) + ADDON);
        }

        if (ETAQ < ETAQP1) {
            if (ETAQP1 > ETAQM1) {
                vs->ETA = ETAQP1;
                vs->NEWQ = static_cast<short>(vs->NQ + 1);
                for (int i = 1; i <= neqs; ++i) {
                    vs->yh[MPL_YH(i, VODE_LMAX)] = vs->acor[i];
                }
            }
            else {
                vs->ETA = ETAQM1;
                vs->NEWQ = static_cast<short>(vs->NQ - 1);
            }
            already_set_eta = true;
        }

        if (ETAQ < ETAQM1 && !already_set_eta) {
            vs->ETA = ETAQM1;
            vs->NEWQ = static_cast<short>(vs->NQ - 1);
            already_set_eta = true;
        }
    }

    if (!already_set_eta) {
        vs->ETA = ETAQ;
        vs->NEWQ = vs->NQ;
    }

    if (vs->ETA >= THRESH && vs->ETAMAX != 1.0) {
        vs->ETA = math::min(vs->ETA, vs->ETAMAX);
        vs->ETA = vs->ETA / math::max(1.0, math::abs(vs->H) * vs->HMXI * vs->ETA);
        vs->NEWH = 1;
        vs->ETAMAX = ETAMX3;
        if (vs->n_step <= 10) {
            vs->ETAMAX = ETAMX2;
        }
        Real R = 1.0 / vs->tq[2];
        for (int i = 1; i <= neqs; ++i) {
            vs->acor[i] *= R;
        }
        return kflag;
    }

    vs->NEWQ = vs->NQ;
    vs->NEWH = 0;
    vs->ETA = 1.0;
    vs->ETAMAX = ETAMX3;
    if (vs->n_step <= 10) {
        vs->ETAMAX = ETAMX2;
    }
    Real R = 1.0 / vs->tq[2];
    for (int i = 1; i <= neqs; ++i) {
        vs->acor[i] *= R;
    }

    return kflag;
}

// from integration/VODE/vode_dvhin.H:15 dvhin
MPL_DEVICE void dvhin(burn_t* state, dvode_t* vs, Real& H0, int& NITER, int& IER,
                      const IntegratorParams* ip, const EosParams* ep, const NetParams* np, EosTables tab)
{
    constexpr Real PT1 = 0.1e0;

    NITER = 0;
    const Real TDIST = math::abs(vs->tout - vs->t);
    const Real TROUND = UROUND * math::max(math::abs(vs->t), math::abs(vs->tout));

    if (TDIST < 2.0 * TROUND) {
        IER = -1;
        return;
    }

    const Real HLB = 100.0 * TROUND;

    Real HUB = PT1 * TDIST;

    for (int i = 1; i <= neqs; ++i) {
        const Real atol = i <= NumSpec ? vs->atol_spec : vs->atol_enuc;
        const Real DELYI = PT1 * math::abs(vs->yh[MPL_YH(i, 1)]) + atol;
        const Real AFI = math::abs(vs->yh[MPL_YH(i, 2)]);
        if (AFI * HUB > DELYI) {
            HUB = DELYI / AFI;
        }
    }

    int iter = 0;
    Real HG = math::sqrt(HLB * HUB);

    bool do_iterations = true;
    if (HUB < HLB) {
        H0 = HG;
        do_iterations = false;
    }

    Real hnew = 0.0;

    if (do_iterations) {

        while (true) {

            const Real H = math::copysign(HG, vs->tout - vs->t);
            for (int i = 1; i <= neqs; ++i) {
                vs->y[i] = vs->yh[MPL_YH(i, 1)] + H * vs->yh[MPL_YH(i, 2)];
            }

            const Real t1 = vs->t + H;

            rhs(t1, state, vs, vs->acor, ip, ep, np, tab);

            for (int i = 1; i <= neqs; ++i) {
                vs->acor[i] = (vs->acor[i] - vs->yh[MPL_YH(i, 2)]) / H;
            }

            Real YDDNRM = 0.0;
            for (int i = 1; i <= neqs; ++i) {
                YDDNRM += (vs->acor[i] * vs->ewt[i]) * (vs->acor[i] * vs->ewt[i]);
            }
            YDDNRM = math::sqrt(YDDNRM / neqs);

            if (YDDNRM * HUB * HUB > 2.0) {
                hnew = math::sqrt(2.0 / YDDNRM);
            }
            else {
                hnew = math::sqrt(HG * HUB);
            }
            iter += 1;

            if (iter >= 4) {
                break;
            }

            const Real HRAT = hnew / HG;
            if ((HRAT > 0.5) && (HRAT < 2.0)) {
                break;
            }

            if ((iter >= 2) && (hnew > 2.0 * HG)) {
                hnew = HG;
                break;
            }

            HG = hnew;
        }

        H0 = hnew * 0.5;
        H0 = math::clamp(H0, HLB, HUB);
    }

    H0 = math::copysign(H0, vs->tout - vs->t);
    NITER = iter;
    IER = 0;
}

// from integration/VODE/vode_dvode.H:23 dvode
// The `return`s inside the step loop set `result` and break.
MPL_DEVICE int dvode(burn_t* state, dvode_t* vs,
                     const IntegratorParams* ip, const EosParams* ep, const NetParams* np, EosTables tab)
{
    Real H0 = 0.0;
    Real S = 0.0;
    int IER = 0;
    int NITER = 0;

    int istate = IERR_SUCCESS;

    vs->HMXI = 1.0 / ip->integrator_ode_max_dt;

    if (vs->tout == vs->t) {
        return istate;
    }

    vs->tn = vs->t;

    vs->n_step = 0;
    vs->n_jac = 0;
    vs->NSLJ = 0;

    Real f_init[neqs + 1] = {};

    rhs(vs->t, state, vs, f_init, ip, ep, np, tab);

    for (int i = 1; i <= neqs; ++i) {
        vs->yh[MPL_YH(i, 2)] = f_init[i];
    }

    vs->n_rhs = 1;

    for (int i = 1; i <= neqs; ++i) {
        vs->yh[MPL_YH(i, 1)] = vs->y[i];
    }

    vs->NQ = 1;
    vs->H = 1.0;

    for (int i = 1; i <= NumSpec; ++i) {
        vs->ewt[i] = vs->rtol_spec * math::abs(vs->yh[MPL_YH(i, 1)]) + vs->atol_spec;
        vs->ewt[i] = 1.0 / vs->ewt[i];
    }
    vs->ewt[NumSpec+1] = vs->rtol_enuc * math::abs(vs->yh[MPL_YH(NumSpec+1, 1)]) + vs->atol_enuc;
    vs->ewt[NumSpec+1] = 1.0 / vs->ewt[NumSpec+1];

    H0 = 0.0;
    dvhin(state, vs, H0, NITER, IER, ip, ep, np, tab);
    vs->n_rhs += NITER;

    if (IER != 0) {
        istate = -3;
        return istate;
    }

    H0 = math::copysign(math::min(math::abs(H0), ip->integrator_ode_max_dt), H0);

    vs->H = H0;
    for (int i = 1; i <= neqs; ++i) {
        vs->yh[MPL_YH(i, 2)] *= H0;
    }

    vs->NQ = 1;
    vs->NEWQ = 1;
    vs->L = 2;
    vs->tau[1] = vs->H;
    vs->PRL1 = 1.0;
    vs->RC = 0.0;
    vs->ETAMAX = 1.0e4;
    vs->NQWAIT = 2;
    vs->HSCAL = vs->H;
    vs->NEWH = 0;
    vs->NSLP = 0;
    vs->IPUP = 1;

    bool skip_loop_start = true;

    int result = IERR_SUCCESS;

    // Single-exit form of the stock `while (true)` loop with its several
    // `return`s and a `continue`: every exit sets `done` after the same side
    // effects, and the code after it is guarded by `if (!done)`. A loop with
    // multiple exits under divergent control fails AMDGPU branch lowering
    // (TARGET/034 masked_region_exits_by_cfg).
    bool done = false;

    while (!done) {

        bool too_many_steps = false;

        if (!skip_loop_start) {

            if (vs->n_step >= ip->integrator_ode_max_steps) {
                for (int i = 1; i <= neqs; ++i) {
                    vs->y[i] = vs->yh[MPL_YH(i, 1)];
                }

                vs->t = vs->tn;

                too_many_steps = true;
            }

            if (!too_many_steps) {
                for (int i = 1; i <= NumSpec; ++i) {
                    vs->ewt[i] = vs->rtol_spec * math::abs(vs->yh[MPL_YH(i, 1)]) + vs->atol_spec;
                    vs->ewt[i] = 1.0 / vs->ewt[i];
                }
                vs->ewt[NumSpec+1] = vs->rtol_enuc * math::abs(vs->yh[MPL_YH(NumSpec+1, 1)]) + vs->atol_enuc;
                vs->ewt[NumSpec+1] = 1.0 / vs->ewt[NumSpec+1];
            }
        }
        else {
            skip_loop_start = false;
        }

        if (too_many_steps) {
            result = IERR_TOO_MANY_STEPS;
            done = true;  // stock: return IERR_TOO_MANY_STEPS
        }

        if (!done) {
            Real TOLSF = 0.0;
            for (int i = 1; i <= neqs; ++i) {
                TOLSF += (vs->yh[MPL_YH(i, 1)] * vs->ewt[i]) * (vs->yh[MPL_YH(i, 1)] * vs->ewt[i]);
            }
            TOLSF = UROUND * math::sqrt(TOLSF / neqs);

            if (TOLSF > 1.0) {

                if (vs->n_step == 0) {
                    result = IERR_TOO_MUCH_ACCURACY_REQUESTED;
                    done = true;  // stock: return
                }
                else {
                    for (int i = 1; i <= neqs; ++i) {
                        vs->y[i] = vs->yh[MPL_YH(i, 1)];
                    }

                    vs->t = vs->tn;

                    result = IERR_TOO_MUCH_ACCURACY_REQUESTED;
                    done = true;  // stock: return
                }
            }
        }

        if (!done) {
            int kflag = dvstep(state, vs, ip, ep, np, tab);

            if (kflag == -1) {
                for (int i = 1; i <= neqs; ++i) {
                    vs->y[i] = vs->yh[MPL_YH(i, 1)];
                }

                vs->t = vs->tn;
                result = IERR_DT_UNDERFLOW;
                done = true;  // stock: return
            }
            else if (kflag == -2) {
                for (int i = 1; i <= neqs; ++i) {
                    vs->y[i] = vs->yh[MPL_YH(i, 1)];
                }

                vs->t = vs->tn;
                result = IERR_CORRECTOR_CONVERGENCE;
                done = true;  // stock: return
            }
            else if (!((vs->tn - vs->tout) * vs->H < 0.0)) {
                // stock: `if ((tn - tout) * H < 0.0) continue;` then finish

                const int Ld = vs->L;
                MPL_ASSUME(1 <= Ld && Ld <= VODE_LMAX);

                for (int i = 1; i <= neqs; ++i) {
                    vs->y[i] = vs->yh[MPL_YH(i, Ld)];
                }

                S = (vs->tout - vs->tn) / vs->H;

                for (int jb = 1; jb <= vs->NQ; ++jb) {
                    const int j = vs->NQ - jb;
                    MPL_ASSUME(0 <= j && j <= VODE_LMAX - 1);
                    for (int i = 1; i <= neqs; ++i) {
                        vs->y[i] = vs->yh[MPL_YH(i, j+1)] + S * vs->y[i];
                    }
                }

                vs->t = vs->tout;

                result = IERR_SUCCESS;
                done = true;  // stock: return IERR_SUCCESS
            }
        }
    }

    return result;
}

// ---------------------------------------------------------------------------
// integration/integrator_setup_strang.H
// ---------------------------------------------------------------------------

// from integration/integrator_setup_strang.H:19 state_backup_t (AMREX_USE_GPU: no xn_in)
struct state_backup_t {
    Real T_in{};
    Real e_in{};
};

// from integration/integrator_setup_strang.H:27 integrator_setup
// (stock returns the value-initialized IntegratorT; here the caller passes a
// value-initialized dvode_t)
MPL_DEVICE void integrator_setup(burn_t* state, Real dt, bool is_retry, dvode_t* vs,
                                 const IntegratorParams* ip, const EosParams* ep, EosTables tab)
{
    if (!is_retry) {
        vs->atol_spec = ip->integrator_atol_spec;
        vs->atol_enuc = ip->integrator_atol_enuc;

        vs->rtol_spec = ip->integrator_rtol_spec;
        vs->rtol_enuc = ip->integrator_rtol_enuc;
    } else {
        vs->atol_spec = ip->integrator_retry_atol_spec > 0 ?
            ip->integrator_retry_atol_spec : ip->integrator_atol_spec;
        vs->atol_enuc = ip->integrator_retry_atol_enuc > 0 ?
            ip->integrator_retry_atol_enuc : ip->integrator_atol_enuc;

        vs->rtol_spec = ip->integrator_retry_rtol_spec > 0 ?
            ip->integrator_retry_rtol_spec : ip->integrator_rtol_spec;
        vs->rtol_enuc = ip->integrator_retry_rtol_enuc > 0 ?
            ip->integrator_retry_rtol_enuc : ip->integrator_rtol_enuc;
    }

    if (is_retry && ip->integrator_retry_swap_jacobian) {
        vs->jacobian_type = (ip->integrator_jacobian == 1) ? 2 : 1;
    } else {
        vs->jacobian_type = static_cast<short>(ip->integrator_jacobian);
    }

    state->success = true;

    vs->t = 0.0;
    vs->tout = dt;

    for (int n = 0; n < SVAR; ++n) {
        state->ydot_a[n] = 0;
    }

    eos(eos_input_rt, state, ep, tab);

    state->e_scale = state->e;

    if (ip->integrator_scale_system) {
        vs->atol_enuc /= state->e_scale;
    }

    burn_to_integrator(state, vs, ip);
}

// from integration/integrator_setup_strang.H:101 integrator_backup
MPL_DEVICE state_backup_t integrator_backup(const burn_t* state)
{
    state_backup_t state_save{};

    state_save.T_in = state->T;
    state_save.e_in = state->e;

    return state_save;
}

// from integration/integrator_setup_strang.H:120 integrator_cleanup
// (no printing on GPU)
MPL_DEVICE void integrator_cleanup(dvode_t* vs, burn_t* state, int istate,
                                   const state_backup_t* state_save, Real dt,
                                   const Real final_state_species_failure_tolerance_factor,
                                   const IntegratorParams* ip)
{
    (void) dt;

    integrator_to_burn(vs, state, ip);

    if (ip->integrator_subtract_internal_energy) {
        state->e -= state_save->e_in;
    }

    if (!ip->integrator_use_number_densities) {
        normalize_abundances_burn(state, ip);
    }

    state->n_rhs = vs->n_rhs;
    state->n_jac = vs->n_jac;
    state->n_step = vs->n_step;

    state->error_code = static_cast<short>(istate);

    if (istate != IERR_SUCCESS) {
        state->success = false;
    }

    const Real final_state_species_failure_tolerance =
        final_state_species_failure_tolerance_factor * ip->integrator_species_failure_tolerance;

    for (int n = 1; n <= NumSpec; ++n) {
        if (vs->y[n] < -final_state_species_failure_tolerance) {
            state->success = false;
        }

        if (!ip->integrator_use_number_densities) {
            if (vs->y[n] > 1.0 + ip->integrator_species_failure_tolerance) {
                state->success = false;
            }
        }
    }
}

} // namespace mpl

#endif
