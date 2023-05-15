#include "hybrj_type.H"
#include <test_func.H>
#include <hybrj.H>

int main() {

    constexpr int neqs = 9;

    hybrj_t<neqs> hj;

    // initial guess

    for (int j = 1; j <= neqs; ++j) {
        hj.x(j) = -1.0;
    }

    hj.xtol = std::sqrt(std::numeric_limits<Real>::epsilon());
    hj.mode = 2;

    for (int j = 1; j <= neqs; ++j) {
        hj.diag(j) = 1.0_rt;
    }

    // we'll pass a single Real as the extra data, but it can be any type

    Real x{-3.14159};

    // hybrj(hj, x, fcn<neqs, Real>, jcn<neqs, Real>);
    hybrj(hj, x);

    std::cout << "done! solution = " << std::endl;
    for (int j = 1; j <= neqs; ++j) {
        std::cout << hj.x(j) << std::endl;
    }

}
