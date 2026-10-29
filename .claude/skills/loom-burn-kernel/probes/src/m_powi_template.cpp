#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// amrex::Math::powi<N> as recursive template, multiplication grouping preserved.
namespace probe {
template<int N> __device__ double powi(double x) {
  if constexpr (N == 0) return 1.0; else if constexpr (N == 1) return x; else if constexpr (N == 2) return x * x;
  else if constexpr (N % 2 == 0) return powi<2>(powi<N / 2>(x)); else return x * powi<N - 1>(x); }
KERNEL void k(double* out, const double* in) { int i = blockIdx.x * blockDim.x + threadIdx.x; out[i] = powi<13>(in[i]) + powi<6>(in[i]); }
}
