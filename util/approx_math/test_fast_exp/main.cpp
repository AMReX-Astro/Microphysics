#include <test_fast_exp.H>

int main() {

    // Accuracy tests

    // Test values including edge cases and some typical values

    test_fast_exp_accuracy(0.0_rt);
    test_fast_exp_accuracy(0.01_rt);
    test_fast_exp_accuracy(0.08_rt);
    test_fast_exp_accuracy(0.1_rt);
    test_fast_exp_accuracy(0.5_rt);
    test_fast_exp_accuracy(1.0_rt);
    test_fast_exp_accuracy(5.0_rt);
    test_fast_exp_accuracy(8.0_rt);
    test_fast_exp_accuracy(10.0_rt);
    test_fast_exp_accuracy(25.0_rt);
    // test_fast_exp_accuracy(100.0_rt);
    // test_fast_exp_accuracy(500.0_rt);
    // test_fast_exp_accuracy(5.0e8_rt);

    test_fast_exp_accuracy(-0.01_rt);
    test_fast_exp_accuracy(-0.08_rt);
    test_fast_exp_accuracy(-0.1_rt);
    test_fast_exp_accuracy(-0.5_rt);
    test_fast_exp_accuracy(-1.0_rt);
    test_fast_exp_accuracy(-5.0_rt);
    test_fast_exp_accuracy(-8.0_rt);
    test_fast_exp_accuracy(-10.0_rt);
    test_fast_exp_accuracy(-25.0_rt);
    // test_fast_exp_accuracy(-100.0_rt);
    // test_fast_exp_accuracy(-500.0_rt);
    // test_fast_exp_accuracy(-5.0e8_rt);

    // Inf Case

    // test_fast_exp_accuracy(std::numeric_limits<amrex::Real>::infinity());
    // test_fast_exp_accuracy(-std::numeric_limits<amrex::Real>::infinity());

    std::cout << "Accuracy tests passed!" << std::endl;

    // Now performance test

    // int iters = 5;
    // amrex::Real test_value = 160.0_rt;
    // test_fast_exp_speed(100, iters, test_value);

    // iters = 10;
    // test_fast_exp_speed(100, iters, test_value);

    // iters = 20;
    // test_fast_exp_speed(100, iters, test_value);

    // iters = 30;
    // test_fast_exp_speed(100, iters, test_value);

    // iters = 50;
    // test_fast_exp_speed(100, iters, test_value);

    // iters = 70;
    // test_fast_exp_speed(100, iters, test_value);
}
