// Loom import/compile smoke test for the complete burn (burner.H): one
// zone per lane, the same call test_react's do_react makes.
#include "../../burner.H"

namespace mpl {

MPL_KERNEL void smoke_burn(const BurnParams* params,
                           const Real* f, const Real* dpdf, const Real* ef, const Real* xf,
                           const Real* d, const Real* t, const Real* dt_sav, const Real* dt2_sav,
                           const Real* dti_sav, const Real* dt2i_sav, const Real* dd_sav,
                           const Real* dd2_sav, const Real* ddi_sav, const Real* dd2i_sav,
                           Real* out, const Real* in, Real dt, int n)
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
        s.e = 0.0;
        s.T_fixed = -1.0;
        burner(&s, dt, params, tab);
        out[3 * idx] = s.e;
        out[3 * idx + 1] = s.T;
        out[3 * idx + 2] = static_cast<Real>(s.n_rhs) + (s.success ? 0.5 : 0.0);
    }
}

} // namespace mpl
