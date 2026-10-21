#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// actual_eos: bool single_iter, converged; int var1, var2; (expected rejected).
namespace probe {
KERNEL void k(double* out, const double* in) { int i = blockIdx.x * blockDim.x + threadIdx.x; bool converged; int var1; if (in[i] > 0.0) { converged = true; var1 = 1; } else { converged = false; var1 = 2; } out[i] = converged ? var1 : -var1; }
}
