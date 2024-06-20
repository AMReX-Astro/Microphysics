#include <test_fast_exp.H>

int main() {

    // Accuracy tests

    // Test values including edge cases and some typical values

    test_fast_exp_accuracy(0.0_rt);
    test_fast_exp_accuracy(0.001_rt);
    test_fast_exp_accuracy(0.01_rt);
    test_fast_exp_accuracy(0.05_rt);
    test_fast_exp_accuracy(0.08_rt);
    test_fast_exp_accuracy(0.1_rt);
    test_fast_exp_accuracy(0.5_rt);
    test_fast_exp_accuracy(1.0_rt);
    test_fast_exp_accuracy(5.0_rt);
    test_fast_exp_accuracy(8.0_rt);
    test_fast_exp_accuracy(10.0_rt);
    test_fast_exp_accuracy(15.0_rt);
    test_fast_exp_accuracy(25.0_rt);
    test_fast_exp_accuracy(100.0_rt);
    test_fast_exp_accuracy(200.0_rt);
    test_fast_exp_accuracy(300.0_rt);

    test_fast_exp_accuracy(-0.001_rt);
    test_fast_exp_accuracy(-0.01_rt);
    test_fast_exp_accuracy(-0.05_rt);
    test_fast_exp_accuracy(-0.08_rt);;
    test_fast_exp_accuracy(-0.1_rt);
    test_fast_exp_accuracy(-0.5_rt);
    test_fast_exp_accuracy(-1.0_rt);
    test_fast_exp_accuracy(-5.0_rt);
    test_fast_exp_accuracy(-8.0_rt);
    test_fast_exp_accuracy(-9.0_rt);
    test_fast_exp_accuracy(-10.0_rt);
    test_fast_exp_accuracy(-11.0_rt);
    test_fast_exp_accuracy(-12.0_rt);
    test_fast_exp_accuracy(-13.0_rt);
    test_fast_exp_accuracy(-14.0_rt);
    test_fast_exp_accuracy(-15.0_rt);
    test_fast_exp_accuracy(-25.0_rt);
    test_fast_exp_accuracy(-100.0_rt);


    // Inf Case

    // test_fast_exp_accuracy(std::numeric_limits<amrex::Real>::infinity());
    // test_fast_exp_accuracy(-std::numeric_limits<amrex::Real>::infinity());

    std::cout << "Accuracy tests passed!" << std::endl;

    // Now performance test

    int iters = 5;
    amrex::Real test_value = 160.0_rt;
    test_fast_exp_speed(100, iters, test_value);

    iters = 10;
    test_fast_exp_speed(100, iters, test_value);

    iters = 20;
    test_fast_exp_speed(100, iters, test_value);

    iters = 30;
    test_fast_exp_speed(100, iters, test_value);

    iters = 50;
    test_fast_exp_speed(100, iters, test_value);

    // iters = 70;
    // test_fast_exp_speed(100, iters, test_value);
}
