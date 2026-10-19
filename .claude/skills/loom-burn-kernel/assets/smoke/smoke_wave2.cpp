// Loom import/compile smoke test for rates.H + screen.H (wave 2): every rate
// index in both instantiations (template-struct recursion, as the network
// module will use), plasma state and screening in both instantiations.
#include "../../burn.H"
#include "../../dual.H"
#include "../../rates.H"
#include "../../screen.H"

namespace mpl {

template <int rate>
struct SmokeRates {
    static MPL_DEVICE void run(const tf_t<Real>* tr, const tf_t<Dual>* td,
                               const RatesParams* p, Real* acc)
    {
        if constexpr (rate <= Rates::NumRates) {
            rate_t<Real> a = evaluate_analytical_rate_r<rate>(tr, p);
            rate_t<Dual> b = evaluate_analytical_rate_d<rate>(td, p);
            *acc = *acc + a.fr + a.rr + b.fr.val + b.fr.grad + b.rr.grad;
            SmokeRates<rate + 1>::run(tr, td, p, acc);
        }
    }
};

MPL_KERNEL void smoke_rates(const RatesParams* params, Real* out, const Real* in, int n)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n) {
        Real T = in[2 * idx];
        Real rho = in[2 * idx + 1];
        tf_t<Real> tr = get_tfactors_r(T);
        tf_t<Dual> td = get_tfactors_d(dual_seed(T));
        Real acc = 0.0;
        SmokeRates<1>::run(&tr, &td, params, &acc);

        Real y[NumSpec + 1] = {};
        for (int k = 1; k <= NumSpec; ++k) {
            y[k] = 1.0 / (4.0 * k);
        }
        plasma_state_t<Real> pr{};
        plasma_state_t<Dual> pd{};
        fill_plasma_state_r(&pr, T, rho, y);
        fill_plasma_state_d(&pd, dual_seed(T), rho, y);
        // network-neutral: the network's first screened pair (network_data.H)
        Real sr = 0.0;
        Dual sd{};
        if constexpr (net_num_screen_pairs > 0) {
            constexpr int pair = 0;
            screen_factors_t sf{net_screen_z1(pair), net_screen_z2(pair), net_screen_a1(pair),
                                net_screen_a2(pair), net_screen_zs53(pair), net_screen_z1_53(pair),
                                net_screen_z2_53(pair), net_screen_zs13(pair), net_screen_zs13inv(pair),
                                net_screen_zhat(pair), net_screen_zhat2(pair), net_screen_lzav(pair),
                                net_screen_aznut(pair)};
            sr = actual_screen_r(&pr, &sf, params);
            sd = actual_screen_d(&pd, &sf, params);
        }
        out[idx] = acc + sr + sd.val + sd.grad;
    }
}

} // namespace mpl
