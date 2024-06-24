#ifndef TEST_FAST_LOG_H
#define TEST_FAST_LOG_H

#include <cmath>
#include <iostream>
#include <limits>
#include <cassert>
#include <chrono>
#include <approx_math.H>


template <typename T>
void test_fast_log_accuracy(T x) {
    //This tests fast_log accuracy

    T fast_log_result = fast_log(x);
    T std_log_result = std::log(x);

    // Print results
    std::cout << "x: " << x
              << " fast_log: " << fast_log_result
              << " std::log: " << std_log_result << std::endl;

    T abs_err = std::abs(fast_log_result - std_log_result);
    T tolerance = 2.e-3_rt;;

    std::cout << "absolute error: " << abs_err << std::endl;

    assert(abs_err < tolerance);
}


template <typename T>
void test_fast_log_speed(int loops, int iter, T x) {
    // This tests fast_log performance

    std::cout << "Testing with loops: " << loops
              << " with iter: " << iter
              << " with initial x: " << x << std::endl;

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
                        result += std::log(x_in);
                        x_in *= 0.99999_rt;
                    }
                }
            }
        }
        // don't let the compiler elide this side-effect-free loop (at the cost of a memory write)
        volatile T volatile_result;
        volatile_result = result;
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<amrex::Real> std_log_duration = end - start;


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
                        result += fast_log(x_in);
                        x_in *= 0.99999_rt;
                    }
                }
            }
        }
        // don't let the compiler elide this side-effect-free loop (at the cost of a memory write)
        volatile T volatile_result;
        volatile_result = result;
    }
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<amrex::Real> fast_log_duration = end - start;


    std::cout << "fast_log duration: " << fast_log_duration.count() << " seconds\n";
    std::cout << "std::log duration: " << std_log_duration.count() << " seconds\n";
}
#endif
