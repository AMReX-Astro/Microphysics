#include <hip/hip_runtime.h>
#define KERNEL __global__ [[loom::workgroup_size(64, 1, 1)]]
// dvode_t-sized automatic struct (aprox13: neqs = 14, VODE_LMAX = 6) with runtime indices.
namespace probe {
struct Vode { double yh[14 * 6]; double jac[14 * 14]; double y[14]; double ewt[14]; double savf[14]; double acor[14]; double el[6]; double tau[6]; double tq[5]; short pivot[14]; int NQ; int L; double H; };
KERNEL void k(double* out, const double* in, const int* idx) { int i = blockIdx.x * blockDim.x + threadIdx.x; Vode v{}; v.NQ = idx[i] % 5 + 1;
  for (int a = 0; a < 14; ++a) { for (int b = 0; b < 14; ++b) { v.jac[a * 14 + b] = in[i] * (a + 1) + b; } v.pivot[a] = static_cast<short>(a); }
  int r = idx[i] % 14; v.yh[r * 6 + v.NQ] = v.jac[r * 14 + r]; out[i] = v.yh[r * 6 + v.NQ] + v.pivot[r]; }
}
