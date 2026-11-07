// Host stand-in for hip_runtime.h: run a probe kernel on the CPU, one thread at a time.
#pragma once
#define __global__
#define __device__
struct probe_dim3 { unsigned x, y, z; };
inline probe_dim3 blockIdx{0, 0, 0}, threadIdx{0, 0, 0}, blockDim{64, 1, 1};
