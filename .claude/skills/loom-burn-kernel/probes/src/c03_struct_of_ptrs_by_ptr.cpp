#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// Tables bundled as a struct of pointers, passed by pointer (EosTables pattern).
namespace probe {
struct Tabs { const double* a; const double* b; };
__device__ double look(const Tabs* t, int i) { return t->a[i] + t->b[i + 1]; }
KERNEL void k(double* out, const double* a, const double* b) { int i = blockIdx.x * blockDim.x + threadIdx.x; Tabs t{a, b}; out[i] = look(&t, i); }
}
