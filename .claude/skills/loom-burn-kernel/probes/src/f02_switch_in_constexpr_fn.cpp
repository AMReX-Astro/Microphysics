#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// aprox13 rhs_data(rate): switch inside a constexpr function returning a struct, used as constexpr.
namespace probe {
struct RhsT { int species_A{-1}; int species_B{-1}; int exponent_A{0}; double branching{1.0}; };
constexpr RhsT rhs_data(int rate) { RhsT d{}; switch (rate) { case 1: d.species_A = 1; d.exponent_A = 3; break; case 2: d.species_A = 2; d.species_B = 1; d.exponent_A = 1; break; default: break; } return d; }
template <int rate> __device__ double use(double x) { constexpr RhsT d = rhs_data(rate); if constexpr (d.exponent_A == 3) return x * x * x * d.branching; else return x * d.species_B; }
KERNEL void k(double* out, const double* in) { int i = blockIdx.x * blockDim.x + threadIdx.x; out[i] = use<1>(in[i]) + use<2>(in[i]); }
}
