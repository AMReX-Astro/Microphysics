// Loom import/compile smoke test for network.H (wave 3): network RHS,
// Jacobian, and the network LU (dgefa/dgesl) on I - h J.
#include "../../burn.H"
#include "../../network.H"

namespace mpl {

MPL_KERNEL void smoke_net(const NetParams* params, Real* out, const Real* in, int n)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n) {
        burn_t s{};
        s.rho = in[4 * idx];
        s.T = in[4 * idx + 1];
        s.eta = in[4 * idx + 2];
        s.cv = in[4 * idx + 3];
        s.y_e = 0.5;
        s.abar = 12.0;
        s.zbar = 6.0;
        for (int k = 0; k < NumSpec; ++k) {
            s.xn[k] = 1.0 / NumSpec;
        }
        Real ydot[neqs + 1] = {};
        Real jac[MPL_JAC_SIZE] = {};
        actual_rhs(&s, params, ydot);
        actual_jac(&s, params, jac);
        Real a[MPL_JAC_SIZE] = {};
        for (int i = 1; i <= neqs; ++i) {
            for (int j = 1; j <= neqs; ++j) {
                Real ident = (i == j) ? 1.0 : 0.0;
                a[MPL_JAC(i, j)] = ident - 1.0e-6 * jac[MPL_JAC(i, j)];
            }
        }
        int info = rhs_dgefa(a);
        rhs_dgesl(a, ydot);
        Real acc = static_cast<Real>(info);
        for (int i = 1; i <= neqs; ++i) {
            acc = acc + ydot[i];
        }
        out[idx] = acc;
    }
}

} // namespace mpl
