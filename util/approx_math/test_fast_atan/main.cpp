#include <test_fast_atan.H>

int main() {

    // Accuracy tests

    // Test values including edge cases and some typical values

    test_fast_atan_accuracy(0.0_rt);
    test_fast_atan_accuracy(0.01_rt);
    test_fast_atan_accuracy(0.05_rt);
    test_fast_atan_accuracy(0.1_rt);
    test_fast_atan_accuracy(0.5_rt);
    test_fast_atan_accuracy(1.0_rt);
    test_fast_atan_accuracy(5.0_rt);
    test_fast_atan_accuracy(8.0_rt);
    test_fast_atan_accuracy(10.0_rt);
    test_fast_atan_accuracy(23.0_rt);
    test_fast_atan_accuracy(100.0_rt);
    test_fast_atan_accuracy(500.0_rt);
    test_fast_atan_accuracy(3000.0_rt);

    test_fast_atan_accuracy(-0.01_rt);
    test_fast_atan_accuracy(-0.05_rt);
    test_fast_atan_accuracy(-0.1_rt);
    test_fast_atan_accuracy(-0.5_rt);
    test_fast_atan_accuracy(-1.0_rt);
    test_fast_atan_accuracy(-5.0_rt);
    test_fast_atan_accuracy(-8.0_rt);
    test_fast_atan_accuracy(-10.0_rt);
    test_fast_atan_accuracy(-23.0_rt);
    test_fast_atan_accuracy(-100.0_rt);
    test_fast_atan_accuracy(-500.0_rt);
    test_fast_atan_accuracy(-3000.0_rt);

    // Inf and NaN case
    test_fast_atan_accuracy(std::numeric_limits<amrex::Real>::infinity());

    std::cout << "Accuracy tests passed!" << std::endl;

    // Now performanc test
    int iters = 300;
    amrex::Real test_value = 0.43_rt;
    test_fast_atan_speed(iters, test_value);
}
