#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// namespace-scope mutable global (runtime parameter) that the root never reads.
namespace probe {
namespace integrator_rp { double rtol_spec = 1.0e-12; }
KERNEL void k(double* out, const double* in) { int i = blockIdx.x * blockDim.x + threadIdx.x; out[i] = in[i] * 2.0; }
}
