#ifndef TEST_FAST_POW_H
#define TEST_FAST_POW_H

#include <cmath>
#include <iostream>
#include <limits>
#include <cassert>
#include <chrono>
#include <approx_math.H>


template <typename T>
void test_fast_pow_accuracy(T x, T y) {
    //This tests fast_pow accuracy

    T fast_pow_result = fast_pow(x, y);
    T std_pow_result = std::pow(x, y);

    // Print results
    std::cout << "x: " << x
              << " fast_pow: " << fast_pow_result
              << " std::pow: " << std_pow_result << std::endl;

    T abs_err = std::abs(fast_pow_result - std_pow_result);
    T rel_err = std::abs(abs_err / std_pow_result);

    T tolerance = 4.0e-2_rt;;

    std::cout << "absolute error: " << abs_err << std::endl;
    std::cout << "relative error: " << rel_err << std::endl;

    assert(rel_err < tolerance);
}


template <typename T>
void test_fast_pow_speed(int loops, int iter, T x, T y) {
    // This tests fast_pow performance

    std::cout << "Testing with loops: " << loops
              << " with iter: " << iter
              << " with initial x: " << x
              << " and initial y: " << y << std::endl;

    auto start = std::chrono::high_resolution_clock::now();
    {
        T x_in;
        T y_in;
        T result;
        for (int m = 0; m < loops; ++m) {
            x_in = x;
            y_in = y;
            result = 0.0;
            for (int i = 0; i < iter; ++i) {
                for (int j = 0; j < iter; ++j) {
                    for (int k = 0; k < iter; ++k) {
                        result += std::pow(x_in, y_in);
                        y_in *= 0.99999_rt;
                    }
                }
            }
        }
        // don't let the compiler elide this side-effect-free loop (at the cost of a memory write)
        volatile T volatile_result;
        volatile_result = result;
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<amrex::Real> std_pow_duration = end - start;


    start = std::chrono::high_resolution_clock::now();
    {
        T x_in;
        T y_in;
        T result;
        for (int m = 0; m < loops; ++m) {
            x_in = x;
            y_in = y;
            result = 0.0;
            for (int i = 0; i < iter; ++i) {
                for (int j = 0; j < iter; ++j) {
                    for (int k = 0; k < iter; ++k) {
                        result += fast_pow(x_in, y_in);
                        y_in *= 0.99999_rt;
                    }
                }
            }
        }
        // don't let the compiler elide this side-effect-free loop (at the cost of a memory write)
        volatile T volatile_result;
        volatile_result = result;
    }
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<amrex::Real> fast_pow_duration = end - start;


    std::cout << "fast_pow duration: " << fast_pow_duration.count() << " seconds\n";
    std::cout << "std::pow duration: " << std_pow_duration.count() << " seconds\n";
}
#endif
