#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// struct local without initializer (`Burn s;`), every field assigned before use.
namespace probe {
struct Burn { double T; double rho; };
KERNEL void k(double* out, const double* in) { int i = blockIdx.x * blockDim.x + threadIdx.x; Burn s; s.T = in[i]; s.rho = 2.0 * in[i]; out[i] = s.T + s.rho; }
}
