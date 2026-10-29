#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// member function called through a pointer to a local struct (`vs->step(x)`).
namespace probe {
struct Acc { double s; int n; __device__ void add(double x) { s += x; ++n; } __device__ double mean() const { return s / n; } };
__device__ static void twice(Acc* a, double x) { a->add(x); a->add(2.0 * x); }
KERNEL void k(double* out, const double* in) { int i = blockIdx.x * blockDim.x + threadIdx.x; Acc a{}; twice(&a, in[i]); out[i] = a.mean(); }
}
