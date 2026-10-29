#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// BurnParams-style POD read through a const pointer (runtime params as kernel input).
namespace probe {
struct Params { double rtol_spec; double atol_spec; int ode_max_steps; int jacobian; bool use_burn_retry; long ncells; };
KERNEL void k(const Params* p, double* out, const double* in) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < p->ncells) { double r = p->use_burn_retry ? p->atol_spec : p->rtol_spec; out[i] = in[i] * r + p->ode_max_steps + p->jacobian; }
}
}
