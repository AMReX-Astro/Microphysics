#ifndef CONSTEXPR_FOR_H
#define CONSTEXPR_FOR_H

// Implementation of "constexpr for" based on
// https://artificial-mind.net/blog/2020/10/31/constexpr-for
//
// Approximates what one would get from a compile-time
// unrolling of the loop
// for (int i = 0; i < N; ++i) {
//    f(i);
// }
//
// The mechanism is recursive, we evaluate f(i) at the current
// i and then call the for loop at i+1. f() is a lambda function
// that provides the body of the loop and takes only an integer
// i as its argument. It is assumed that at the loop call site
// the lambda is declared [&] so that it can operate on whatever
// data it needs.

template<auto I, auto N, class F>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
constexpr void constexpr_for (F&& f)
{
    if constexpr (I < N) {
        f(std::integral_constant<decltype(I), I>());
        constexpr_for<I+1, N>(f);
    }
}

#endif
