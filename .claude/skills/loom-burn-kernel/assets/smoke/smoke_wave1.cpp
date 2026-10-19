// Loom import/compile smoke test for the wave-1 modules (eos.H, dual.H, the
// math recipes). Not a correctness test: it only checks that the translated
// code is accepted by loom-import-cxx and loom-compile for gfx90a.
#include "../../burn.H"
#include "../../eos.H"
#include "../../dual.H"

namespace mpl {

// Table structs must reach device code by value (a struct of pointers read
// through a pointer fails import: probes c03/c06), so the kernel takes each
// table pointer as its own argument.
MPL_KERNEL void smoke_eos(const EosParams* params,
                          const Real* f, const Real* dpdf, const Real* ef, const Real* xf,
                          const Real* d, const Real* t, const Real* dt_sav, const Real* dt2_sav,
                          const Real* dti_sav, const Real* dt2i_sav, const Real* dd_sav,
                          const Real* dd2_sav, const Real* ddi_sav, const Real* dd2i_sav,
                          Real* out, const Real* in, int n)
{
    EosTables tab{f, dpdf, ef, xf, d, t, dt_sav, dt2_sav, dti_sav, dt2i_sav,
                  dd_sav, dd2_sav, ddi_sav, dd2i_sav};
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n) {
        burn_t s{};
        s.rho = in[2 * idx];
        s.T = in[2 * idx + 1];
        for (int k = 0; k < NumSpec; ++k) {
            s.xn[k] = 1.0 / NumSpec;
        }
        eos(eos_input_rt, &s, params, tab);
        out[4 * idx] = s.e;
        out[4 * idx + 1] = s.cv;
        eos(eos_input_re, &s, params, tab);
        out[4 * idx + 2] = s.T;
    }
}

MPL_KERNEL void smoke_math(Real* out, const Real* in, int n)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n) {
        Real x = in[idx];
        Dual t = dual_seed(x);
        Dual r = dual_mul(dual_exp(t), dual_const(math::pow(x, 1.5)));
        out[idx] = r.val + r.grad + math::log10(x) + math::tanh(x) + math::cbrt(x);
    }
}

} // namespace mpl
