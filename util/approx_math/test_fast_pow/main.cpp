#include <test_fast_pow.H>

int main() {

    // Accuracy tests

    // Test values including edge cases and some typical values

    // test_fast_pow_accuracy(-1.0_rt);
    // test_fast_pow_accuracy(0.0_rt);
    test_fast_pow_accuracy(0.001_rt, 3.14_rt);
    test_fast_pow_accuracy(0.01_rt, 23.5_rt);
    test_fast_pow_accuracy(0.1_rt, 10.4_rt);
    test_fast_pow_accuracy(0.5_rt, 59.323_rt);
    test_fast_pow_accuracy(1.0_rt, 7.32_rt);
    test_fast_pow_accuracy(3.14_rt, 88.4_rt);
    test_fast_pow_accuracy(5.0_rt, 23.4_rt);
    test_fast_pow_accuracy(10.0_rt, -12.8_rt);
    test_fast_pow_accuracy(15.0_rt, 75.9_rt);
    test_fast_pow_accuracy(20.0_rt, -3.1_rt);
    test_fast_pow_accuracy(30.0_rt, -8.9_rt);
    test_fast_pow_accuracy(50.0_rt, -2.4_rt);
    test_fast_pow_accuracy(100.0_rt, 7.4_rt);
    test_fast_pow_accuracy(500.0_rt, 5.22_rt);
    // test_fast_pow_accuracy(std::numeric_limits<amrex::Real>::infinity());


    std::cout << "Accuracy tests passed!" << std::endl;

    // Now performance test

    int iters = 5;
    amrex::Real x_value = 3.14_rt;
    amrex::Real y_value = 122.22_rt;
    test_fast_pow_speed(100, iters, x_value, y_value);

    iters = 10;
    test_fast_pow_speed(100, iters, x_value, y_value);

    iters = 20;
    test_fast_pow_speed(100, iters, x_value, y_value);

    iters = 30;
    test_fast_pow_speed(100, iters, x_value, y_value);

    iters = 50;
    test_fast_pow_speed(100, iters, x_value, y_value);

    iters = 70;
    test_fast_pow_speed(100, iters, x_value, y_value);
}
