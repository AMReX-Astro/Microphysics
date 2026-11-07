#pragma once
#include <cstdio>
#include <cstring>
#include <cstdint>
constexpr int N = 64;
inline void fill_inputs(double* in, int* sp) {
    for (int i = 0; i < N; ++i) {
        in[i] = (i - 20) * 0.37 + 1.0 / (i + 3);   // negatives, values above and below 1
        sp[i] = 7 * i + 1;
    }
}
inline void dump(const char* path, const double* out) {
    FILE* f = std::fopen(path, "w");
    for (int i = 0; i < N; ++i) {
        std::uint64_t b; std::memcpy(&b, &out[i], 8);
        std::fprintf(f, "%2d %016llx %.17g\n", i, (unsigned long long)b, out[i]);
    }
    std::fclose(f);
}
