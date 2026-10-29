#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// assumption over a struct field read through a pointer (`__builtin_assume(1 <= vs->L && vs->L <= 6)`) guarding `vs->el[vs->L]`.
namespace probe {
struct Dvode { int L; double el[7]; };
__device__ static double pick(const Dvode* vs) { __builtin_assume(1 <= vs->L && vs->L <= 6); return vs->el[vs->L]; }
KERNEL void k(double* out, const double* in, const int* l_in) { int i = blockIdx.x * blockDim.x + threadIdx.x; Dvode vs{}; for (int n = 0; n < 7; ++n) vs.el[n] = in[i] + n; vs.L = l_in[i]; out[i] = pick(&vs); }
}
