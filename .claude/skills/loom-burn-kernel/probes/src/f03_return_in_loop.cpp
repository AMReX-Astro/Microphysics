#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// VODE dvnlsd/dvstep/dvode: return from inside while(true) (expected rejected).
namespace probe {
__device__ int iterate(double* y, int maxit) { int m = 0; while (true) { y[0] *= 0.5; ++m; if (y[0] < 1.0) return 0; if (m >= maxit) return -1; } }
KERNEL void k(double* out, const double* in) { int i = blockIdx.x * blockDim.x + threadIdx.x; double y[1] = {in[i]}; int r = iterate(y, 50); out[i] = y[0] + r; }
}
