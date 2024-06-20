#ifndef TEST_FAST_ATAN_H
#define TEST_FAST_ATAN_H

#include <cmath>
#include <iostream>
#include <limits>
#include <cassert>
#include <chrono>
#include <approx_math.H>


template <typename T>
void test_fast_atan_accuracy(T x) {
    //This tests fast_atan accuracy

    T fast_atan_result = fast_atan(x);
    T std_atan_result = std::atan(x);

    // Print results
    std::cout << "x: " << x
              << " fast_atan: " << fast_atan_result
              << " std::atan: " << std_atan_result << std::endl;

    T abs_err = std::abs(fast_atan_result - std_atan_result);
    T tolerance = 2.e-3_rt;;

    std::cout << "absolute error: " << abs_err << std::endl;

    assert(abs_err < tolerance);
}


template <typename T>
void test_fast_atan_speed(int iter, T x) {
    // This tests fast_atan performance

    std::cout << "Testing with iter: " << iter
              << " with initial x: " << x << std::endl;

    volatile amrex::Real x_in = x;
    volatile T result = 0.0_rt;
    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < iter; ++i) {
        for (int j = 0; j < iter; ++j) {
            for (int k = 0; k < iter; ++k) {
                result += fast_atan(x_in);
                x_in -= 1e-8_rt * x_in;
            }
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<amrex::Real> fast_atan_duration = end - start;
    // std::cout << "Result from fast_atan is " << result << std::endl;

    x_in = x;
    result = 0.0_rt;
    start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < iter; ++i) {
        for (int j = 0; j < iter; ++j) {
            for (int k = 0; k < iter; ++k) {
                result += std::atan(x_in);
                x_in -= 1e-8_rt * x_in;
            }
        }
    }
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<amrex::Real> std_atan_duration = end - start;
    // std::cout << "Result from std::atan is " << result << std::endl;

    std::cout << "fast_atan duration: " << fast_atan_duration.count() << " seconds\n";
    std::cout << "std::atan duration: " << std_atan_duration.count() << " seconds\n";
}
#endif
