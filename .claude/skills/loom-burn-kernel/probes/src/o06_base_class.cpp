#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// GradArray : Array1D style inheritance (expected rejected).
namespace probe {
struct Base { double a[3]; };
struct Grad : Base { double scale; };
KERNEL void k(double* out, const double* in) { int i = blockIdx.x * blockDim.x + threadIdx.x; Grad g{}; g.a[0] = in[i]; g.scale = 2.0; out[i] = g.a[0] * g.scale; }
}
