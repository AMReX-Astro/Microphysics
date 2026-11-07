// Build: c++ -std=c++20 -O0 -ffp-contract=off -I fakehip -include <probe.cpp> -DNARGS=2|3 host_driver.cpp
#include "common.H"
int main(int argc, char** argv) {
    double in[N], out[N]; int sp[N];
    fill_inputs(in, sp);
    for (int t = 0; t < N; ++t) {
        out[t] = -1.0;
    }
    for (int t = 0; t < N; ++t) {
        threadIdx.x = t;
#if NARGS == 3
        probe::k(out, in, sp);
#else
        probe::k(out, in);
#endif
    }
    dump(argv[1], out);
}
