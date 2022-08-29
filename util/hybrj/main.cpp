#include "hybrj_type.H"

#include <hybrj.H>

template<int neqs>
void fcn(Array1D<Real, 1, neqs>& x, Array1D<Real, 1, neqs>& fvec, int& iflag) {

    for (int k = 1; k <= neqs; ++k) {
        Real temp = (3.0_rt - 2.0_rt * x(k)) * x(k);
        Real temp1 = 0.0_rt;
        if (k != 1) {
            temp1 = x(k-1);
        }
        Real temp2 = 0.0_rt;
        if (k != neqs) {
            temp2 = x(k+1);
        }
        fvec(k) = temp - temp1 - 2.0_rt*temp2 + 1.0_rt;
    }
}

template<int neqs>
void jcn(Array1D<Real, 1, neqs>& x, Array2D<Real, 1, neqs, 1, neqs>& fjac, int& iflag) {

    for (int k = 1; k <= neqs; ++k) {
        for (int j = 1; j <= neqs; ++j) {
            fjac(k, j) = 0.0_rt;
        }
        fjac(k, k) = 3.0_rt - 4.0_rt * x(k);
        if (k != 1) {
            fjac(k, k-1) = -1.0_rt;
        }
        if (k != neqs) {
            fjac(k, k+1) = -2.0_rt;
        }
    }

}

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

    hybrj(hj, fcn<neqs>, jcn<neqs>);

    std::cout << "done! solution = " << std::endl;
    for (int j = 1; j <= neqs; ++j) {
        std::cout << hj.x(j) << std::endl;
    }

}
