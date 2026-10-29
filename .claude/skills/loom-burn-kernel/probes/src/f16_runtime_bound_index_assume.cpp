#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// same as f15 with MPL_ASSUME-style bounds on every derived index (current rewrite).
namespace probe {
struct Dvode { int NQ; int L; double el[14]; double yh[13 * 7]; };
__device__ static void adv(Dvode* vs) { const int NQ = vs->NQ; __builtin_assume(1 <= NQ && NQ <= 5); for (int j = 1; j <= NQ - 2; ++j) { __builtin_assume(1 <= j && j <= 3); for (int i = j + 1; i >= 2; --i) { __builtin_assume(2 <= i && i <= 6); vs->el[i] = vs->el[i] * 0.5 + vs->el[i - 1]; } } const int L = vs->L; __builtin_assume(1 <= L && L <= 6); for (int e = 0; e < 13; ++e) vs->yh[e * 7 + (L - 1)] += vs->el[L]; }
KERNEL void k(double* out, const double* in, const int* nq_in) { int i = blockIdx.x * blockDim.x + threadIdx.x; Dvode vs{}; for (int n = 0; n < 14; ++n) vs.el[n] = in[i] + n; vs.NQ = nq_in[i]; vs.L = vs.NQ + 1; adv(&vs); out[i] = vs.el[2] + vs.yh[6]; }
}
