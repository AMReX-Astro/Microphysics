#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// T x{} with default member initializers (eos_re_extra_t, dvode_t{}).
namespace probe {
struct E { double rho{0.0}; double T{0.0}; double dedT{0.0}; int n{0}; double xn[13]{}; };
KERNEL void k(double* out, const double* in) { int i = blockIdx.x * blockDim.x + threadIdx.x; E e{}; e.rho = in[i]; e.xn[2] = 1.0; out[i] = e.rho + e.xn[2] + e.n; }
}
