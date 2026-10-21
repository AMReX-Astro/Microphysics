#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// All tables in one buffer at compile-time offsets; one base pointer passed down.
namespace probe {
constexpr int OFF_A = 0;
constexpr int OFF_B = 1024;
__device__ double look(const double* tab, int i) { return tab[OFF_A + i] + tab[OFF_B + i + 1]; }
KERNEL void k(double* out, const double* tab) { int i = blockIdx.x * blockDim.x + threadIdx.x; out[i] = look(tab, i); }
}
