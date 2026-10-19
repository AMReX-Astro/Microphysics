#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// Struct of pointers living in device memory, kernel gets a pointer to it.
namespace probe {
struct Tabs { const double* a; const double* b; };
KERNEL void k(double* out, const Tabs* t) { int i = blockIdx.x * blockDim.x + threadIdx.x; out[i] = t->a[i] + t->b[i + 1]; }
}
