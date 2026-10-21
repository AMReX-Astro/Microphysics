#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// enum NetworkRates : std::uint8_t, written without std headers.
namespace probe {
enum Rates : unsigned char { He4_He4_He4_to_C12 = 1, C12_He4_to_O16, NumRates = C12_He4_to_O16 };
template <int r> __device__ double f(double x) { if constexpr (r == C12_He4_to_O16) return 2.0 * x; else return x; }
KERNEL void k(double* out, const double* in) { int i = blockIdx.x * blockDim.x + threadIdx.x; out[i] = f<He4_He4_He4_to_C12>(in[i]) + f<C12_He4_to_O16>(in[i]) + NumRates; }
}
