// Generic reaction-network machinery for the Loom burn kernel (see README.md):
// networks/rhs.H (NEW_NETWORK_IMPLEMENTATION: rhs, jac, construct_rate,
// rhs_term, jac_term, apply_density_scaling, apply_screening,
// apply_branching, fill_additional_rates, intermediate rates, is_rate_used,
// dgefa, dgesl), for STRANG, no NEUTRINOS, HIP branches.
//
// Network-independent: every network-specific value comes from the generated
// network_data.H (rd_* = rhs_data fields, net_intermediate_index, aion_inv,
// zion, screen factors), and every network-specific function from hooks the
// including network.H defines BEFORE including this file:
//   struct NetParams; RatesParams net_rates_params(const NetParams&)
//   template <int rate> rate_t<Real> net_analytical_rate_r(const burn_t&, const NetYArray1D& y,
//                                                            const tf_t<Real>&, const RatesParams&)
//   template <int rate> rate_t<Dual> net_analytical_rate_d(same, tf_t<Dual>)
//       actual_network.H evaluate_analytical_rate
//   template <int rate> void postprocess_rate_r(const burn_t&, const NetYArray1D& y, rate_t<Real>* rates,
//                                               rate_t<Real>* rates1, rate_t<Real>* rates2,
//                                               rate_t<Real>* rates3)   (and _d with rate_t<Dual>)
//       actual_network.H postprocess_rate (rates by pointer: it writes their
//       fields, which a mutable struct reference does not import)
//   template <int spec> Real ener_gener_rate(const burn_t&, const Real dydt)
//       actual_network.H ener_gener_rate
// NetYArray1D (rhs_state.y, Array1D<Real, 1, NumSpec>) and JacNetArray2D come
// from the generated burn.H, so network.H can use them.
// rhs_state_t is not materialized: its fields are burn_state.rho / y_e /
// eta (the stock code copies them from burn_state), y(1..NumSpec), tf and
// pstate, passed separately.
//
// Not ported: network.use_tables (rate_tab_t, ttab, rattab, tabulate_rates,
// evaluate_tabulated_rate). Only the analytic branch of construct_rate is
// translated; the host launcher must reject use_tables != 0
// (NetParams::network_use_tables documents the requirement). NEUTRINOS and
// NSE are refused at compile time.
//
// Entry points (called by the integrator module):
//   actual_rhs(const burn_t& state, const NetParams& params, Array1D<Real, 1, neqs>& ydot)
//       RHS::rhs<neqs>: ydot(1..neqs), Y/s and erg/g/s, exactly as stock
//       (before integrator scaling).
//   actual_jac(const burn_t& state, const NetParams& params, JacNetArray2D& jac)
//       RHS::jac: jac(i, j), i, j = 1..neqs.
//   rhs_dgefa(JacNetArray2D& a) -> info     RHS::dgefa (LU without pivoting)
//   rhs_dgesl(const JacNetArray2D& a, Array1D<Real, 1, neqs>& b)   RHS::dgesl
//
// Arrays keep the stock 1-based accessors (support.H Array1D / Array2D);
// JacNetArray2D = MathArray2D<jac_t, 1, neqs, 1, neqs>.
//
// Loops over rates/species: stock amrex::constexpr_for lambdas become
// template-struct recursion (README rule 5) in the stock order: rates
// 1..NumRates, within a rate species 1..NumSpec, Jacobian rows spec1 then
// columns spec2. Every rhs_data value used is bound to a constexpr scalar
// (README rule 6).
#ifndef MPL_NETWORK_RHS_H
#define MPL_NETWORK_RHS_H

#include "support.H"
#include "dual.H"
#include "network_data.H"
#include "burn.H"
#include "rates.H"
#include "screen.H"

#if MPL_GEN_NEUTRINOS
#error "network_rhs.H: NEUTRINOS is defined, but neutrino cooling (sneut5) is not ported; regenerate network_data.H with USE_NEUTRINOS=FALSE"
#endif
#if MPL_GEN_NSE || MPL_GEN_NSE_NET
#error "network_rhs.H: NSE / NSE_NET builds are not ported"
#endif

namespace mpl {


// from networks/rhs.H:35 factorial (as a loop; only folded)
MPL_DEVICE constexpr int net_factorial(int n)
{
    int r = 1;
    for (int k = 2; k <= n; ++k) {
        r = r * k;
    }
    return r;
}

// from networks/rhs.H:49 is_rate_used<species, rate>()
MPL_DEVICE constexpr int is_rate_used(int species, int rate)
{
    const int sA = rd_species(rate, 0);
    const int sB = rd_species(rate, 1);
    const int sC = rd_species(rate, 2);
    const int sD = rd_species(rate, 3);
    const int sE = rd_species(rate, 4);
    const int sF = rd_species(rate, 5);
    int used = 0;
    if (sA == species || sB == species || sC == species ||
        sD == species || sE == species || sF == species) {
        if (sA > NumSpec || sB > NumSpec || sC > NumSpec ||
            sD > NumSpec || sE > NumSpec || sF > NumSpec) {
            used = 0;
        } else {
            used = 1;
        }
    }
    return used;
}

// from networks/rhs.H:358 density_exponent_forward<rate>()
MPL_DEVICE constexpr int density_exponent_forward(int rate)
{
    int exponent = 0;
    if (rd_species(rate, 0) >= 0) { exponent += rd_exponent(rate, 0); }
    if (rd_species(rate, 1) >= 0) { exponent += rd_exponent(rate, 1); }
    if (rd_species(rate, 2) >= 0) { exponent += rd_exponent(rate, 2); }
    if (exponent > 0) { exponent -= 1; }
    return exponent;
}

// from networks/rhs.H:384 density_exponent_reverse<rate>()
MPL_DEVICE constexpr int density_exponent_reverse(int rate)
{
    int exponent = 0;
    if (rd_species(rate, 3) >= 0) { exponent += rd_exponent(rate, 3); }
    if (rd_species(rate, 4) >= 0) { exponent += rd_exponent(rate, 4); }
    if (rd_species(rate, 5) >= 0) { exponent += rd_exponent(rate, 5); }
    if (exponent > 0) { exponent -= 1; }
    return exponent;
}

// networks/rhs.H:1251 "We cannot have a zero-sized array, so just set the
// array size to 1 in that case."
constexpr int net_intermediate_array_size = net_num_intermediate > 0 ? net_num_intermediate : 1;

template <typename T>
using net_inter_t = Array1D<rate_t<T>, 1, net_intermediate_array_size>;

// ---------------------------------------------------------------------------
// networks/rhs.H:410 apply_density_scaling, :450 apply_screening,
// :548 apply_branching, :1129 construct_rate
// rhs_state_t is not materialized: burn_state (rho, y_e, eta), y, tf and
// pstate are passed separately (see the header comment).
// ---------------------------------------------------------------------------

// from networks/rhs.H:410 apply_density_scaling (number_t = Real)
template <int rate>
MPL_DEVICE void apply_density_scaling_r(const Real rho, rate_t<Real>* rates)
{
    constexpr int forward_exponent = density_exponent_forward(rate);
    constexpr int reverse_exponent = density_exponent_reverse(rate);

    static_assert(forward_exponent <= 3);
    static_assert(reverse_exponent <= 3);

    Real density_term_forward = 1.0;
    Real density_term_reverse = 1.0;

    if constexpr (forward_exponent == 1) {
        density_term_forward = rho;
    } else if constexpr (forward_exponent == 2) {
        density_term_forward = rho * rho;
    } else if constexpr (forward_exponent == 3) {
        density_term_forward = rho * rho * rho;
    }

    if constexpr (reverse_exponent == 1) {
        density_term_reverse = rho;
    } else if constexpr (reverse_exponent == 2) {
        density_term_reverse = rho * rho;
    } else if constexpr (reverse_exponent == 3) {
        density_term_reverse = rho * rho * rho;
    }

    rates->fr *= density_term_forward;
    rates->rr *= density_term_reverse;
}

// from networks/rhs.H:410 apply_density_scaling (number_t = autodiff::dual)
template <int rate>
MPL_DEVICE void apply_density_scaling_d(const Real rho, rate_t<Dual>* rates)
{
    constexpr int forward_exponent = density_exponent_forward(rate);
    constexpr int reverse_exponent = density_exponent_reverse(rate);

    static_assert(forward_exponent <= 3);
    static_assert(reverse_exponent <= 3);

    Real density_term_forward = 1.0;
    Real density_term_reverse = 1.0;

    if constexpr (forward_exponent == 1) {
        density_term_forward = rho;
    } else if constexpr (forward_exponent == 2) {
        density_term_forward = rho * rho;
    } else if constexpr (forward_exponent == 3) {
        density_term_forward = rho * rho * rho;
    }

    if constexpr (reverse_exponent == 1) {
        density_term_reverse = rho;
    } else if constexpr (reverse_exponent == 2) {
        density_term_reverse = rho * rho;
    } else if constexpr (reverse_exponent == 3) {
        density_term_reverse = rho * rho * rho;
    }

    // dual *= Real: compound assignment through a full dual temporary
    rates->fr = dual_mul(rates->fr, dual_const(density_term_forward));
    rates->rr = dual_mul(rates->rr, dual_const(density_term_reverse));
}

#if MPL_GEN_SCREENING
// from networks/rhs.H:450 apply_screening (number_t = Real)
template <int rate>
MPL_DEVICE void apply_screening_r(const plasma_state_t<Real>& pstate, rate_t<Real>* rates,
                                  const RatesParams& rp)
{
    constexpr int screen_fwd = rd_screen_forward_reaction(rate);
    constexpr int screen_rev = rd_screen_reverse_reaction(rate);
    constexpr int sA = rd_species(rate, 0);
    constexpr int sB = rd_species(rate, 1);
    constexpr int eA = rd_exponent(rate, 0);
    constexpr int eB = rd_exponent(rate, 1);
    constexpr int eC = rd_exponent(rate, 2);

    // stock: `if constexpr (no screening) return;` before the three blocks
    if constexpr (!(screen_fwd == 0 && screen_rev == 0)) {

        if constexpr (eA == 1 && eB == 1 && eC == 0) {
            // Forward reaction is A + B, screen using these two species
            constexpr int Z1 = net_zion_i(sA);
            constexpr int A1 = net_aion_i(sA);
            constexpr int Z2 = net_zion_i(sB);
            constexpr int A2 = net_aion_i(sB);

            const screen_factors_t scn_fac = screen_factor<Z1, A1, Z2, A2>();

            const Real sc = actual_screen_r(pstate, scn_fac, rp);

            if constexpr (screen_fwd == 1) {
                rates->fr *= sc;
            }
            if constexpr (screen_rev == 1) {
                rates->rr *= sc;
            }
        }

        if constexpr (eA == 2 && eB == 0 && eC == 0) {
            // Forward reaction is A + A, screen using just this species
            constexpr int Z1 = net_zion_i(sA);
            constexpr int A1 = net_aion_i(sA);

            const screen_factors_t scn_fac = screen_factor<Z1, A1, Z1, A1>();

            const Real sc = actual_screen_r(pstate, scn_fac, rp);

            if constexpr (screen_fwd == 1) {
                rates->fr *= sc;
            }
            if constexpr (screen_rev == 1) {
                rates->rr *= sc;
            }
        }

        if constexpr (eA == 3 && eB == 0 && eC == 0) {
            // triple alpha: A + A, then A + X with X = (2 Z1, 2 A1)
            constexpr int Z1 = net_zion_i(sA);
            constexpr int A1 = net_aion_i(sA);

            const screen_factors_t scn_fac1 = screen_factor<Z1, A1, Z1, A1>();

            const Real sc1 = actual_screen_r(pstate, scn_fac1, rp);

            const screen_factors_t scn_fac2 = screen_factor<Z1, A1, 2 * Z1, 2 * A1>();

            const Real sc2 = actual_screen_r(pstate, scn_fac2, rp);

            const Real sc = sc1 * sc2;

            if constexpr (screen_fwd == 1) {
                rates->fr *= sc;
            }
            if constexpr (screen_rev == 1) {
                rates->rr *= sc;
            }
        }
    }
}

// from networks/rhs.H:450 apply_screening (number_t = autodiff::dual)
template <int rate>
MPL_DEVICE void apply_screening_d(const plasma_state_t<Dual>& pstate, rate_t<Dual>* rates,
                                  const RatesParams& rp)
{
    constexpr int screen_fwd = rd_screen_forward_reaction(rate);
    constexpr int screen_rev = rd_screen_reverse_reaction(rate);
    constexpr int sA = rd_species(rate, 0);
    constexpr int sB = rd_species(rate, 1);
    constexpr int eA = rd_exponent(rate, 0);
    constexpr int eB = rd_exponent(rate, 1);
    constexpr int eC = rd_exponent(rate, 2);

    if constexpr (!(screen_fwd == 0 && screen_rev == 0)) {

        if constexpr (eA == 1 && eB == 1 && eC == 0) {
            constexpr int Z1 = net_zion_i(sA);
            constexpr int A1 = net_aion_i(sA);
            constexpr int Z2 = net_zion_i(sB);
            constexpr int A2 = net_aion_i(sB);

            const screen_factors_t scn_fac = screen_factor<Z1, A1, Z2, A2>();

            const Dual sc = actual_screen_d(pstate, scn_fac, rp);

            if constexpr (screen_fwd == 1) {
                rates->fr = dual_mul(rates->fr, sc);
            }
            if constexpr (screen_rev == 1) {
                rates->rr = dual_mul(rates->rr, sc);
            }
        }

        if constexpr (eA == 2 && eB == 0 && eC == 0) {
            constexpr int Z1 = net_zion_i(sA);
            constexpr int A1 = net_aion_i(sA);

            const screen_factors_t scn_fac = screen_factor<Z1, A1, Z1, A1>();

            const Dual sc = actual_screen_d(pstate, scn_fac, rp);

            if constexpr (screen_fwd == 1) {
                rates->fr = dual_mul(rates->fr, sc);
            }
            if constexpr (screen_rev == 1) {
                rates->rr = dual_mul(rates->rr, sc);
            }
        }

        if constexpr (eA == 3 && eB == 0 && eC == 0) {
            constexpr int Z1 = net_zion_i(sA);
            constexpr int A1 = net_aion_i(sA);

            const screen_factors_t scn_fac1 = screen_factor<Z1, A1, Z1, A1>();

            const Dual sc1 = actual_screen_d(pstate, scn_fac1, rp);

            const screen_factors_t scn_fac2 = screen_factor<Z1, A1, 2 * Z1, 2 * A1>();

            const Dual sc2 = actual_screen_d(pstate, scn_fac2, rp);

            const Dual sc = dual_mul(sc1, sc2);

            if constexpr (screen_fwd == 1) {
                rates->fr = dual_mul(rates->fr, sc);
            }
            if constexpr (screen_rev == 1) {
                rates->rr = dual_mul(rates->rr, sc);
            }
        }
    }
}

#endif // MPL_GEN_SCREENING

// from networks/rhs.H:548 apply_branching (number_t = Real)
template <int rate>
MPL_DEVICE void apply_branching_r(rate_t<Real>* rates)
{
    constexpr Real forward_branching_ratio = rd_forward_branching_ratio(rate);
    constexpr Real reverse_branching_ratio = rd_reverse_branching_ratio(rate);

    if constexpr (forward_branching_ratio != 1.0) {
        rates->fr *= forward_branching_ratio;
    }
    if constexpr (reverse_branching_ratio != 1.0) {
        rates->rr *= reverse_branching_ratio;
    }
}

// from networks/rhs.H:548 apply_branching (number_t = autodiff::dual)
template <int rate>
MPL_DEVICE void apply_branching_d(rate_t<Dual>* rates)
{
    constexpr Real forward_branching_ratio = rd_forward_branching_ratio(rate);
    constexpr Real reverse_branching_ratio = rd_reverse_branching_ratio(rate);

    if constexpr (forward_branching_ratio != 1.0) {
        rates->fr = dual_mul(rates->fr, dual_const(forward_branching_ratio));
    }
    if constexpr (reverse_branching_ratio != 1.0) {
        rates->rr = dual_mul(rates->rr, dual_const(reverse_branching_ratio));
    }
}

// from networks/rhs.H:1129 construct_rate (number_t = Real; analytic branch only)
template <int rate>
MPL_DEVICE rate_t<Real> construct_rate_r(const burn_t& burn_state, const NetYArray1D& y, const tf_t<Real>& tf,
                                         const plasma_state_t<Real>& pstate, const RatesParams& rp)
{
    rate_t<Real> rates{0.0, 0.0};

    // use_tables && rate_can_be_tabulated -> evaluate_tabulated_rate: not ported
    rates = net_analytical_rate_r<rate>(burn_state, y, tf, rp);

    // Set the density dependence
    apply_density_scaling_r<rate>(burn_state.rho, &rates);

#if MPL_GEN_SCREENING
    // Screen
    apply_screening_r<rate>(pstate, &rates, rp);
#endif

    // Branching ratios
    apply_branching_r<rate>(&rates);

    return rates;
}

// from networks/rhs.H:1129 construct_rate (number_t = autodiff::dual; analytic branch only)
template <int rate>
MPL_DEVICE rate_t<Dual> construct_rate_d(const burn_t& burn_state, const NetYArray1D& y, const tf_t<Dual>& tf,
                                         const plasma_state_t<Dual>& pstate, const RatesParams& rp)
{
    rate_t<Dual> rates{dual_const(0.0), dual_const(0.0)};

    rates = net_analytical_rate_d<rate>(burn_state, y, tf, rp);

    apply_density_scaling_d<rate>(burn_state.rho, &rates);

#if MPL_GEN_SCREENING
    apply_screening_d<rate>(pstate, &rates, rp);
#endif

    apply_branching_d<rate>(&rates);

    return rates;
}

// from networks/rhs.H:1163 fill_additional_rates (number_t = Real)
template <int rate>
MPL_DEVICE void fill_additional_rates_r(const net_inter_t<Real>& intermediate_rates, rate_t<Real>& rates1,
                                        rate_t<Real>& rates2, rate_t<Real>& rates3)
{
    constexpr int rate1 = rd_additional_reaction(rate, 1);
    constexpr int rate2 = rd_additional_reaction(rate, 2);
    constexpr int rate3 = rd_additional_reaction(rate, 3);

    if constexpr (rate1 >= 0) {
        constexpr int index = net_intermediate_index(rate1);
        static_assert(index >= 1 && index <= net_num_intermediate);
        rates1 = intermediate_rates(index);
    }
    if constexpr (rate2 >= 0) {
        constexpr int index = net_intermediate_index(rate2);
        static_assert(index >= 1 && index <= net_num_intermediate);
        rates2 = intermediate_rates(index);
    }
    if constexpr (rate3 >= 0) {
        constexpr int index = net_intermediate_index(rate3);
        static_assert(index >= 1 && index <= net_num_intermediate);
        rates3 = intermediate_rates(index);
    }
}

// from networks/rhs.H:1163 fill_additional_rates (number_t = autodiff::dual)
template <int rate>
MPL_DEVICE void fill_additional_rates_d(const net_inter_t<Dual>& intermediate_rates, rate_t<Dual>& rates1,
                                        rate_t<Dual>& rates2, rate_t<Dual>& rates3)
{
    constexpr int rate1 = rd_additional_reaction(rate, 1);
    constexpr int rate2 = rd_additional_reaction(rate, 2);
    constexpr int rate3 = rd_additional_reaction(rate, 3);

    if constexpr (rate1 >= 0) {
        constexpr int index = net_intermediate_index(rate1);
        static_assert(index >= 1 && index <= net_num_intermediate);
        rates1 = intermediate_rates(index);
    }
    if constexpr (rate2 >= 0) {
        constexpr int index = net_intermediate_index(rate2);
        static_assert(index >= 1 && index <= net_num_intermediate);
        rates2 = intermediate_rates(index);
    }
    if constexpr (rate3 >= 0) {
        constexpr int index = net_intermediate_index(rate3);
        static_assert(index >= 1 && index <= net_num_intermediate);
        rates3 = intermediate_rates(index);
    }
}

// from integration/utils/jacobian_utilities.H:18 temperature_to_energy_jacobian
MPL_DEVICE Real temperature_to_energy_jacobian(const burn_t& state, const Real jac_T)
{
    Real jac_e = jac_T / state.cv;
    return jac_e;
}

// ---------------------------------------------------------------------------
// networks/rhs.H:662 rhs_term and :875 jac_term. Both only read the rate
// values: rhs_term gets (fr, rr) for number_t = Real and
// (derivative(fr), derivative(rr)) for number_t = dual; jac_term gets
// (val(fr), val(rr)).
// ---------------------------------------------------------------------------

struct net_terms_t {
    Real forward_term;
    Real reverse_term;
};

// The per-species factor of rhs_term: Y**exponent / exponent!
// (networks/rhs.H:678-697 for A; B..F identical)
template <int rate, int slot>
MPL_DEVICE Real rhs_term_factor(const burn_t& state)
{
    constexpr int sp = rd_species(rate, slot);
    constexpr int exponent = rd_exponent(rate, slot);
    constexpr int ipf = rd_apply_identical_particle_factor(rate);
    constexpr Real ainv = net_aion_inv(sp);

    Real Y = state.xn[sp - 1] * ainv;

    static_assert(exponent <= 3);

    Real df = 1.0;

    if constexpr (exponent == 1) {
        df = Y;
    } else if constexpr (exponent == 2) {
        df = Y * Y;
    } else if constexpr (exponent == 3) {
        df = Y * Y * Y;
    }

    if constexpr (ipf != 0) {
        constexpr int identical_particle_factor = net_factorial(exponent);
        df *= 1.0 / static_cast<Real>(identical_particle_factor);
    }

    return df;
}

// from networks/rhs.H:662 rhs_term<species, rate, number_t>
template <int species, int rate>
MPL_DEVICE net_terms_t rhs_term(const burn_t& state, const Real fr, const Real rr)
{
    constexpr int sA = rd_species(rate, 0);
    constexpr int sB = rd_species(rate, 1);
    constexpr int sC = rd_species(rate, 2);
    constexpr int sD = rd_species(rate, 3);
    constexpr int sE = rd_species(rate, 4);
    constexpr int sF = rd_species(rate, 5);

    Real forward_term = fr;

    if constexpr (sA >= 0) {
        Real df = rhs_term_factor<rate, 0>(state);
        forward_term *= df;
    }
    if constexpr (sB >= 0) {
        Real df = rhs_term_factor<rate, 1>(state);
        forward_term *= df;
    }
    if constexpr (sC >= 0) {
        Real df = rhs_term_factor<rate, 2>(state);
        forward_term *= df;
    }

    Real reverse_term = rr;

    if constexpr (sD >= 0) {
        Real dr = rhs_term_factor<rate, 3>(state);
        reverse_term *= dr;
    }
    if constexpr (sE >= 0) {
        Real dr = rhs_term_factor<rate, 4>(state);
        reverse_term *= dr;
    }
    if constexpr (sF >= 0) {
        Real dr = rhs_term_factor<rate, 5>(state);
        reverse_term *= dr;
    }

    // Now compute the total contribution to this species.

    if constexpr (sA == species) {
        constexpr Real n = rd_number(rate, 0);
        forward_term *= -n;
        reverse_term *= n;
    } else if constexpr (sB == species) {
        constexpr Real n = rd_number(rate, 1);
        forward_term *= -n;
        reverse_term *= n;
    } else if constexpr (sC == species) {
        constexpr Real n = rd_number(rate, 2);
        forward_term *= -n;
        reverse_term *= n;
    } else if constexpr (sD == species) {
        constexpr Real n = rd_number(rate, 3);
        forward_term *= n;
        reverse_term *= -n;
    } else if constexpr (sE == species) {
        constexpr Real n = rd_number(rate, 4);
        forward_term *= n;
        reverse_term *= -n;
    } else if constexpr (sF == species) {
        constexpr Real n = rd_number(rate, 5);
        forward_term *= n;
        reverse_term *= -n;
    } else {
        forward_term = 0.0;
        reverse_term = 0.0;
    }

    net_terms_t t{forward_term, reverse_term};
    return t;
}

// The per-species factor of jac_term (networks/rhs.H:887-922 for A; B..F
// identical): d/dY(spec2) of Y**exponent / exponent! when spec2 is this
// species, Y**exponent / exponent! otherwise.
template <int spec2, int rate, int slot>
MPL_DEVICE Real jac_term_factor(const burn_t& state)
{
    constexpr int sp = rd_species(rate, slot);
    constexpr int exponent = rd_exponent(rate, slot);
    constexpr int ipf = rd_apply_identical_particle_factor(rate);
    constexpr Real ainv = net_aion_inv(sp);

    Real Y = state.xn[sp - 1] * ainv;

    constexpr int exp = (spec2 == sp) ? exponent - 1 : exponent;

    static_assert(exp <= 3);

    Real df = 1.0;

    if constexpr (exp == 1) {
        df = Y;
    } else if constexpr (exp == 2) {
        df = Y * Y;
    } else if constexpr (exp == 3) {
        df = Y * Y * Y;
    }

    if constexpr (spec2 == sp) {
        df *= exponent;
    }

    if constexpr (ipf != 0) {
        constexpr int identical_particle_factor = net_factorial(exponent);
        df *= 1.0 / static_cast<Real>(identical_particle_factor);
    }

    return df;
}

// from networks/rhs.H:875 jac_term<spec1, spec2, rate>
template <int spec1, int spec2, int rate>
MPL_DEVICE Real jac_term(const burn_t& state, const Real fr_val, const Real rr_val)
{
    constexpr int sA = rd_species(rate, 0);
    constexpr int sB = rd_species(rate, 1);
    constexpr int sC = rd_species(rate, 2);
    constexpr int sD = rd_species(rate, 3);
    constexpr int sE = rd_species(rate, 4);
    constexpr int sF = rd_species(rate, 5);
    constexpr int used = is_rate_used(spec1, rate);

    Real forward_term = 0.0;

    if constexpr (used != 0 && (spec2 == sA || spec2 == sB || spec2 == sC)) {

        forward_term = fr_val;

        if constexpr (sA >= 0) {
            Real df = jac_term_factor<spec2, rate, 0>(state);
            forward_term *= df;
        }
        if constexpr (sB >= 0) {
            Real df = jac_term_factor<spec2, rate, 1>(state);
            forward_term *= df;
        }
        if constexpr (sC >= 0) {
            Real df = jac_term_factor<spec2, rate, 2>(state);
            forward_term *= df;
        }
    }

    Real reverse_term = 0.0;

    if constexpr (used != 0 && (spec2 == sD || spec2 == sE || spec2 == sF)) {

        reverse_term = rr_val;

        if constexpr (sD >= 0) {
            Real dr = jac_term_factor<spec2, rate, 3>(state);
            reverse_term *= dr;
        }
        if constexpr (sE >= 0) {
            Real dr = jac_term_factor<spec2, rate, 4>(state);
            reverse_term *= dr;
        }
        if constexpr (sF >= 0) {
            Real dr = jac_term_factor<spec2, rate, 5>(state);
            reverse_term *= dr;
        }
    }

    // Now compute the total contribution to this species.

    Real term = 0.0;

    if constexpr (sA == spec1) {
        constexpr Real n = rd_number(rate, 0);
        term = n * (reverse_term - forward_term);
    }
    if constexpr (sB == spec1) {
        constexpr Real n = rd_number(rate, 1);
        term = n * (reverse_term - forward_term);
    }
    if constexpr (sC == spec1) {
        constexpr Real n = rd_number(rate, 2);
        term = n * (reverse_term - forward_term);
    }
    if constexpr (sD == spec1) {
        constexpr Real n = rd_number(rate, 3);
        term = n * (forward_term - reverse_term);
    }
    if constexpr (sE == spec1) {
        constexpr Real n = rd_number(rate, 4);
        term = n * (forward_term - reverse_term);
    }
    if constexpr (sF == spec1) {
        constexpr Real n = rd_number(rate, 5);
        term = n * (forward_term - reverse_term);
    }

    return term;
}

// ---------------------------------------------------------------------------
// networks/rhs.H:1216 rhs<neqs>
// ---------------------------------------------------------------------------

// rhs.H:1257 "Fill all intermediate rates first" (number_t = Real)
template <int rate>
struct NetInterR {
    static MPL_DEVICE void run(const burn_t& burn_state, const NetYArray1D& y, const tf_t<Real>& tf,
                               const plasma_state_t<Real>& pstate, const RatesParams& rp,
                               net_inter_t<Real>& intermediate_rates)
    {
        if constexpr (rate <= Rates::NumRates) {
            constexpr int index = net_intermediate_index(rate);
            if constexpr (index >= 1) {
                intermediate_rates(index) = construct_rate_r<rate>(burn_state, y, tf, pstate, rp);
            }
            NetInterR<rate + 1>::run(burn_state, y, tf, pstate, rp, intermediate_rates);
        }
    }
};

// rhs.H:1300 species loop of one rate (nrhs == neqs)
template <int species, int rate>
struct NetRhsSpecies {
    static MPL_DEVICE void run(const burn_t& burn_state, const rate_t<Real> rates, Array1D<Real, 1, neqs>& ydot)
    {
        if constexpr (species <= NumSpec) {
            if constexpr (is_rate_used(species, rate) != 0) {
                const net_terms_t t = rhs_term<species, rate>(burn_state, rates.fr, rates.rr);
                ydot(species) += t.forward_term + t.reverse_term;
            }
            NetRhsSpecies<species + 1, rate>::run(burn_state, rates, ydot);
        }
    }
};

// rhs.H:1270 loop over all rates
template <int rate>
struct NetRhsRate {
    static MPL_DEVICE void run(const burn_t& burn_state, const NetYArray1D& y, const tf_t<Real>& tf,
                               const plasma_state_t<Real>& pstate, const RatesParams& rp,
                               const net_inter_t<Real>& intermediate_rates, Array1D<Real, 1, neqs>& ydot)
    {
        if constexpr (rate <= Rates::NumRates) {
            rate_t<Real> rates{0.0, 0.0};

            constexpr int index = net_intermediate_index(rate);
            if constexpr (index < 0) {
                rates = construct_rate_r<rate>(burn_state, y, tf, pstate, rp);
            } else {
                rates = intermediate_rates(index);
            }

            rate_t<Real> rates1{0.0, 0.0};
            rate_t<Real> rates2{0.0, 0.0};
            rate_t<Real> rates3{0.0, 0.0};

            fill_additional_rates_r<rate>(intermediate_rates, rates1, rates2, rates3);

            postprocess_rate_r<rate>(burn_state, y, &rates, &rates1, &rates2, &rates3);

            NetRhsSpecies<1, rate>::run(burn_state, rates, ydot);

            NetRhsRate<rate + 1>::run(burn_state, y, tf, pstate, rp, intermediate_rates, ydot);
        }
    }
};

// rhs.H:1344 energy RHS term
template <int species>
struct NetRhsEner {
    static MPL_DEVICE void run(const burn_t& burn_state, Array1D<Real, 1, neqs>& ydot)
    {
        if constexpr (species <= NumSpec) {
            ydot(net_ienuc) += ener_gener_rate<species>(burn_state, ydot(species));
            NetRhsEner<species + 1>::run(burn_state, ydot);
        }
    }
};

// from networks/rhs.H:1216 rhs<nrhs = neqs> (RHS::rhs, actual_rhs)
// ydot(1..neqs)
MPL_DEVICE void actual_rhs(const burn_t& burn_state, const NetParams& params, Array1D<Real, 1, neqs>& ydot)
{
    const RatesParams rp = net_rates_params(params);

    // rhs_state_t<amrex::Real>: rho, eta, y_e are burn_state's; y below

    // network_properties.H aion_inv[0..NumSpec-1]
    const Real aion_inv[NumSpec] = {MPL_NET_AION_INV_LIST};

    // Convert X to Y.
    NetYArray1D y{};
    for (int n = 1; n <= NumSpec; ++n) {
        y(n) = burn_state.xn[n-1] * aion_inv[n-1];
    }

    // Set up the state data, which is the same for all screening factors.
    plasma_state_t<Real> pstate{};
#if MPL_GEN_SCREENING
    fill_plasma_state_r(&pstate, burn_state.T, burn_state.rho, y);
#endif

    // Initialize the rate temperature term.
    const tf_t<Real> tf = get_tfactors_r(burn_state.T);

    // Initialize the RHS terms.
    for (int n = 1; n <= neqs; ++n) {
        ydot(n) = 0.0;
    }

    net_inter_t<Real> intermediate_rates{};

    // Fill all intermediate rates first.
    NetInterR<1>::run(burn_state, y, tf, pstate, rp, intermediate_rates);

    // Loop over all rates, and then loop over all species, and if the
    // rate affects that given species, add its contribution to the RHS.
    NetRhsRate<1>::run(burn_state, y, tf, pstate, rp, intermediate_rates, ydot);

    // Evaluate the neutrino cooling (NEUTRINOS off).
    Real sneut = 0.0;

    // Compute the energy RHS term.
    ydot(net_ienuc) = -sneut;

    NetRhsEner<1>::run(burn_state, ydot);
}

// ---------------------------------------------------------------------------
// networks/rhs.H:1360 jac
// ---------------------------------------------------------------------------

// rhs.H:1407 "Fill all intermediate rates first" (number_t = autodiff::dual)
template <int rate>
struct NetInterD {
    static MPL_DEVICE void run(const burn_t& burn_state, const NetYArray1D& y, const tf_t<Dual>& tf,
                               const plasma_state_t<Dual>& pstate, const RatesParams& rp,
                               net_inter_t<Dual>& intermediate_rates)
    {
        if constexpr (rate <= Rates::NumRates) {
            constexpr int index = net_intermediate_index(rate);
            if constexpr (index >= 1) {
                intermediate_rates(index) = construct_rate_d<rate>(burn_state, y, tf, pstate, rp);
            }
            NetInterD<rate + 1>::run(burn_state, y, tf, pstate, rp, intermediate_rates);
        }
    }
};

// rhs.H:1453 columns spec2 of row spec1 (only instantiated when
// is_rate_used<spec1, rate>, which the stock code tests per element)
template <int spec1, int spec2, int rate>
struct NetJacCol {
    static MPL_DEVICE void run(const burn_t& burn_state, const rate_t<Dual> rates, JacNetArray2D& jac)
    {
        if constexpr (spec2 <= NumSpec) {
            jac(spec1, spec2) += jac_term<spec1, spec2, rate>(burn_state, rates.fr.val, rates.rr.val);
            NetJacCol<spec1, spec2 + 1, rate>::run(burn_state, rates, jac);
        }
    }
};

// rhs.H:1448 species Jacobian elements with respect to other species
template <int spec1, int rate>
struct NetJacRow {
    static MPL_DEVICE void run(const burn_t& burn_state, const rate_t<Dual> rates, JacNetArray2D& jac)
    {
        if constexpr (spec1 <= NumSpec) {
            if constexpr (is_rate_used(spec1, rate) != 0) {
                NetJacCol<spec1, 1, rate>::run(burn_state, rates, jac);
            }
            NetJacRow<spec1 + 1, rate>::run(burn_state, rates, jac);
        }
    }
};

// rhs.H:1464 Jacobian elements with respect to temperature
template <int species, int rate>
struct NetJacTemp {
    static MPL_DEVICE void run(const burn_t& burn_state, const rate_t<Dual> rates, JacNetArray2D& jac)
    {
        if constexpr (species <= NumSpec) {
            if constexpr (is_rate_used(species, rate) != 0) {
                const net_terms_t t = rhs_term<species, rate>(burn_state, rates.fr.grad, rates.rr.grad);
                jac(species, net_ienuc) += t.forward_term + t.reverse_term;
            }
            NetJacTemp<species + 1, rate>::run(burn_state, rates, jac);
        }
    }
};

// rhs.H:1419 loop over rates and compute Jacobian terms
// (rates1..3 live in jac(), outside the loop, as in the stock code)
template <int rate>
struct NetJacRate {
    static MPL_DEVICE void run(const burn_t& burn_state, const NetYArray1D& y, const tf_t<Dual>& tf,
                               const plasma_state_t<Dual>& pstate, const RatesParams& rp,
                               const net_inter_t<Dual>& intermediate_rates, rate_t<Dual>& rates1,
                               rate_t<Dual>& rates2, rate_t<Dual>& rates3, JacNetArray2D& jac)
    {
        if constexpr (rate <= Rates::NumRates) {
            rate_t<Dual> rates{dual_const(0.0), dual_const(0.0)};

            constexpr int index = net_intermediate_index(rate);
            if constexpr (index < 0) {
                rates = construct_rate_d<rate>(burn_state, y, tf, pstate, rp);
            } else {
                rates = intermediate_rates(index);
            }

            fill_additional_rates_d<rate>(intermediate_rates, rates1, rates2, rates3);

            postprocess_rate_d<rate>(burn_state, y, &rates, &rates1, &rates2, &rates3);

            NetJacRow<1, rate>::run(burn_state, rates, jac);

            NetJacTemp<1, rate>::run(burn_state, rates, jac);

            NetJacRate<rate + 1>::run(burn_state, y, tf, pstate, rp, intermediate_rates,
                                      rates1, rates2, rates3, jac);
        }
    }
};

// rhs.H:1500 jac(net_ienuc, species) += ener_gener_rate<s>(rhs_state, jac(s, species)), s = 1..NumSpec
template <int s, int species>
struct NetJacEnerRow {
    static MPL_DEVICE void run(const burn_t& burn_state, JacNetArray2D& jac)
    {
        if constexpr (s <= NumSpec) {
            jac(net_ienuc, species) += ener_gener_rate<s>(burn_state, jac(s, species));
            NetJacEnerRow<s + 1, species>::run(burn_state, jac);
        }
    }
};

// rhs.H:1488 energy generation rate Jacobian elements
template <int species>
struct NetJacEner {
    static MPL_DEVICE void run(const burn_t& burn_state, const Real dsnuda, const Real dsnudz, JacNetArray2D& jac)
    {
        if constexpr (species <= NumSpec) {
            constexpr Real zion = net_zion(species);

            // Energy generation rate Jacobian elements with respect to species.
            Real b1 = (-burn_state.abar * burn_state.abar * dsnuda +
                       (zion - burn_state.zbar) * burn_state.abar * dsnudz);
            jac(net_ienuc, species) = -b1;

            NetJacEnerRow<1, species>::run(burn_state, jac);

            // Convert previously computed terms from d/dT to d/de.
            jac(species, net_ienuc) =
                temperature_to_energy_jacobian(burn_state, jac(species, net_ienuc));

            // Compute df(e) / de term.
            jac(net_ienuc, net_ienuc) += ener_gener_rate<species>(burn_state, jac(species, net_ienuc));

            NetJacEner<species + 1>::run(burn_state, dsnuda, dsnudz, jac);
        }
    }
};

// from networks/rhs.H:1360 jac (RHS::jac, actual_jac)
// jac(i, j), i, j = 1..neqs
MPL_DEVICE void actual_jac(const burn_t& burn_state, const NetParams& params, JacNetArray2D& jac)
{
    const RatesParams rp = net_rates_params(params);

    // rhs_state_t<autodiff::dual>: rho, eta, y_e are burn_state's; y below

    // network_properties.H aion_inv[0..NumSpec-1]
    const Real aion_inv[NumSpec] = {MPL_NET_AION_INV_LIST};

    // Convert X to Y.
    NetYArray1D y{};
    for (int n = 1; n <= NumSpec; ++n) {
        y(n) = burn_state.xn[n-1] * aion_inv[n-1];
    }

    // seed the dual number for temperature before calculating anything with it
    const Dual dual_temp = dual_seed(burn_state.T);

    // Set up the state data, which is the same for all screening factors.
    plasma_state_t<Dual> pstate{};
#if MPL_GEN_SCREENING
    fill_plasma_state_d(&pstate, dual_temp, burn_state.rho, y);
#endif

    // Initialize the rate temperature term.
    const tf_t<Dual> tf = get_tfactors_d(dual_temp);

    // Initialize the Jacobian terms.
    for (int i = 1; i <= neqs; ++i) {
        for (int j = 1; j <= neqs; ++j) {
            jac(i, j) = 0.0;
        }
    }

    net_inter_t<Dual> intermediate_rates{};

    rate_t<Dual> rates1{};
    rate_t<Dual> rates2{};
    rate_t<Dual> rates3{};

    // Fill all intermediate rates first.
    NetInterD<1>::run(burn_state, y, tf, pstate, rp, intermediate_rates);

    // Loop over rates and compute Jacobian terms.
    NetJacRate<1>::run(burn_state, y, tf, pstate, rp, intermediate_rates,
                       rates1, rates2, rates3, jac);

    // Evaluate the neutrino cooling (NEUTRINOS off).
    Real dsneutdt = 0.0;
    Real dsnuda = 0.0;
    Real dsnudz = 0.0;

    jac(net_ienuc, net_ienuc) = -temperature_to_energy_jacobian(burn_state, dsneutdt);

    NetJacEner<1>::run(burn_state, dsnuda, dsnudz, jac);
}

// ---------------------------------------------------------------------------
// networks/rhs.H:276 dgesl, :310 dgefa (RHS::, called by VODE dvjac/dvnlsd
// for NEW_NETWORK_IMPLEMENTATION networks; no pivoting). Plain loops with the
// constexpr_for bounds written as strict upper bounds (`j < neqs + 1`, not
// `j <= neqs`): Loom proves those counted-loop indices in range, so no
// MPL_ASSUME is needed (with `<=`, dgesl's back substitution fails
// SUBRANGE/026); a(i, j) and b(i), i, j = 1..neqs.
// ---------------------------------------------------------------------------

// from networks/rhs.H:276 dgesl
MPL_DEVICE void rhs_dgesl(const JacNetArray2D& a, Array1D<Real, 1, neqs>& b)
{
    // solve a * x = b
    // first solve l * y = b
    for (int k = 1; k < neqs; ++k) {
        Real t = b(k);
        for (int j = k + 1; j < neqs + 1; ++j) {
            b(j) += t * a(j, k);
        }
    }

    // now solve u * x = y
    for (int kb = 1; kb < neqs + 1; ++kb) {
        const int k = neqs + 1 - kb;

        b(k) = b(k) / a(k, k);
        Real t = -b(k);

        for (int j = 1; j < k; ++j) {
            b(j) += t * a(j, k);
        }
    }
}

// from networks/rhs.H:310 dgefa: LU factorization in place without pivoting.
// Returns info (the last k with a zero pivot, 0 if none).
MPL_DEVICE int rhs_dgefa(JacNetArray2D& a)
{
    int info = 0;

    for (int k = 1; k < neqs; ++k) {
        // compute multipliers
        if (a(k, k) == 0.0) {
            info = k;
            // stock: `return;` from the constexpr_for body, i.e. continue
        } else {
            Real t = -1.0 / a(k, k);
            for (int j = k + 1; j < neqs + 1; ++j) {
                a(j, k) *= t;
            }

            // row elimination with column indexing
            for (int j = k + 1; j < neqs + 1; ++j) {
                t = a(k, j);
                for (int i = k + 1; i < neqs + 1; ++i) {
                    a(i, j) += t * a(i, k);
                }
            }
        }
    }

    return info;
}

} // namespace mpl

#endif
