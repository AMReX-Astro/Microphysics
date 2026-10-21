#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// constant-bound loop with a data-dependent `break` indexing a private array (pivot search shape), no assumptions.
namespace probe {
KERNEL void k(double* out, const double* in) { int i = blockIdx.x * blockDim.x + threadIdx.x; double a[14] = {}; for (int n = 0; n < 14; ++n) a[n] = in[i] - n; int piv = 13; for (int j = 1; j <= 13; ++j) { a[j] += a[j - 1]; if (a[j] < 0.0) { piv = j; break; } } out[i] = a[piv]; }
}
