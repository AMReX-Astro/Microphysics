#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// Array1D::operator() style: free function returning T& (expected rejected).
namespace probe {
struct A { double d[14]; };
static double& at(A& a, int i) { return a.d[i - 1]; }
KERNEL void k(double* out, const double* in) { int i = blockIdx.x * blockDim.x + threadIdx.x; A a{}; at(a, 2) = in[i]; out[i] = a.d[1]; }
}
