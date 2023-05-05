#ifndef TEST_FUNC_H
#define TEST_FUNC_H

template<int neqs, typename T>
void fcn(Array1D<Real, 1, neqs>& x, Array1D<Real, 1, neqs>& fvec, const T& data, int& iflag) {

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

template<int neqs, typename T>
void jcn(Array1D<Real, 1, neqs>& x, Array2D<Real, 1, neqs, 1, neqs>& fjac, const T& data, int& iflag) {

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
#endif
