#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// Array1D<Real,1,N>::operator()(int) returning T& (1-based) (expected rejected).
namespace probe {
struct A1 { double arr[14]; __device__ double& operator()(int i) { return arr[i - 1]; } };
KERNEL void k(double* out, const double* in) { int i = blockIdx.x * blockDim.x + threadIdx.x; A1 y{}; y(1) = in[i]; out[i] = y(1) * 2.0; }
}
