#ifndef TEST_FAST_EXP_H
#define TEST_FAST_EXP_H

#include <cmath>
#include <iostream>
#include <limits>
#include <cassert>
#include <chrono>
#include <approx_math.H>


template <typename T>
void test_fast_exp_accuracy(T x) {
    //This tests fast_exp accuracy

    T fast_exp_result = fast_exp(x);
    T std_exp_result = std::exp(x);

    // Print results
    std::cout << "x: " << x
              << " fast_exp: " << fast_exp_result
              << " std::exp: " << std_exp_result << std::endl;

    T abs_err = std::abs(fast_exp_result - std_exp_result);
    T rel_err = abs_err / std_exp_result;
    auto rtol = static_cast<T>(0.05);

    std::cout << "absolute error: " << abs_err << std::endl;
    std::cout << "relative error: " << rel_err << std::endl;

    assert(rel_err < rtol);
}


template <typename T>
void test_fast_exp_speed(int loops, int iter, T x) {
    // This tests fast_exp performance

    std::cout << "Testing with loops: " << loops
              << " with iter: " << iter
              << " with initial x: " << x << std::endl;

    auto fac = static_cast<T>(1.e-4);

    auto start = std::chrono::high_resolution_clock::now();
    {
        T x_in;
        T result;
        for (int m = 0; m < loops; ++m) {
            x_in = x;
            result = 0.0;
            for (int i = 0; i < iter; ++i) {
                for (int j = 0; j < iter; ++j) {
                    for (int k = 0; k < iter; ++k) {
                        result += std::exp(x_in);
                        x_in -= fac * x_in;
                    }
                }
            }
        }
        // don't let the compiler elide this side-effect-free loop (at the cost of a memory write)
        volatile T volatile_result;
        volatile_result = result;
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<amrex::Real> std_exp_duration = end - start;


    start = std::chrono::high_resolution_clock::now();
    {
        T x_in;
        T result;
        for (int m = 0; m < loops; ++m) {
            x_in = x;
            result = 0.0;
            for (int i = 0; i < iter; ++i) {
                for (int j = 0; j < iter; ++j) {
                    for (int k = 0; k < iter; ++k) {
                        result += fast_exp(x_in);
                        x_in -= fac * x_in;
                    }
                }
            }
        }
        // don't let the compiler elide this side-effect-free loop (at the cost of a memory write)
        volatile T volatile_result;
        volatile_result = result;
    }
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<amrex::Real> fast_exp_duration = end - start;


    std::cout << "fast_exp duration: " << fast_exp_duration.count() << " seconds\n";
    std::cout << "std::exp duration: " << std_exp_duration.count() << " seconds\n";
}
#endif
