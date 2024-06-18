#include <test_fast_atan.H>

int main() {

    // Accuracy tests

    // Test values including edge cases and some typical values

    test_fast_atan_accuracy(0.0_rt);
    test_fast_atan_accuracy(0.01_rt);
    test_fast_atan_accuracy(-0.01_rt);
    test_fast_atan_accuracy(0.1_rt);
    test_fast_atan_accuracy(-0.1_rt);
    test_fast_atan_accuracy(0.5_rt);
    test_fast_atan_accuracy(-0.5_rt);
    test_fast_atan_accuracy(1.0_rt);
    test_fast_atan_accuracy(-1.0_rt);
    test_fast_atan_accuracy(10.0_rt);
    test_fast_atan_accuracy(-10.0_rt);
    test_fast_atan_accuracy(100.0_rt);
    test_fast_atan_accuracy(-100.0_rt);
    test_fast_atan_accuracy(-1000.0_rt);
    test_fast_atan_accuracy(-1000.0_rt);

    // Large positive and negative values
    test_fast_atan_accuracy(std::numeric_limits<amrex::Real>::infinity());
    test_fast_atan_accuracy(-std::numeric_limits<amrex::Real>::infinity());

    // NaN value
    test_fast_atan_accuracy(std::numeric_limits<amrex::Real>::quiet_NaN());

    std::cout << "Accuracy tests passed!" << std::endl;

    // Now performanc test
    int iters = 2000;
    amrex::Real test_value = 4146.3232_rt;
    test_fast_atan_speed(iters, test_value);
}
