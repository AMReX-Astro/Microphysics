#include <test_atanf.H>

int main() {

    // Accuracy tests

    // Test values including edge cases and some typical values

    test_atanf_accuracy(0.0_rt);
    test_atanf_accuracy(0.01_rt);
    test_atanf_accuracy(0.08_rt);
    test_atanf_accuracy(0.1_rt);
    test_atanf_accuracy(0.5_rt);
    test_atanf_accuracy(1.0_rt);
    test_atanf_accuracy(5.0_rt);
    test_atanf_accuracy(8.0_rt);
    test_atanf_accuracy(10.0_rt);
    test_atanf_accuracy(25.0_rt);
    test_atanf_accuracy(100.0_rt);
    test_atanf_accuracy(500.0_rt);
    test_atanf_accuracy(5.0e8_rt);

    test_atanf_accuracy(-0.01_rt);
    test_atanf_accuracy(-0.08_rt);
    test_atanf_accuracy(-0.1_rt);
    test_atanf_accuracy(-0.5_rt);
    test_atanf_accuracy(-1.0_rt);
    test_atanf_accuracy(-5.0_rt);
    test_atanf_accuracy(-8.0_rt);
    test_atanf_accuracy(-10.0_rt);
    test_atanf_accuracy(-25.0_rt);
    test_atanf_accuracy(-100.0_rt);
    test_atanf_accuracy(-500.0_rt);
    test_atanf_accuracy(-5.0e8_rt);

    // Inf Case

    test_atanf_accuracy(std::numeric_limits<amrex::Real>::infinity());
    test_atanf_accuracy(-std::numeric_limits<amrex::Real>::infinity());

    std::cout << "Accuracy tests passed!" << std::endl;

    // Now performanc test
    int iters = 5;
    amrex::Real test_value = 1.0e7_rt;
    test_atanf_speed(iters, test_value);

    iters = 10;
    test_atanf_speed(iters, test_value);

    iters = 50;
    test_atanf_speed(iters, test_value);

    iters = 75;
    test_atanf_speed(iters, test_value);

    iters = 100;
    test_atanf_speed(iters, test_value);

    iters = 300;
    test_atanf_speed(iters, test_value);

    iters = 500;
    test_atanf_speed(iters, test_value);
}
