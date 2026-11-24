#ifndef TEST_FUNC_H
#define TEST_FUNC_H

template<int neqs, typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void fcn(amrex::Array1D<amrex::Real, 1, neqs>& x, amrex::Array1D<amrex::Real, 1, neqs>& fvec, const T& data, int& iflag) {

    for (int k = 1; k <= neqs; ++k) {
        amrex::Real temp = (3.0_rt - 2.0_rt * x(k)) * x(k);
        amrex::Real temp1 = 0.0_rt;
        if (k != 1) {
            temp1 = x(k-1);
        }
        amrex::Real temp2 = 0.0_rt;
        if (k != neqs) {
            temp2 = x(k+1);
        }
        fvec(k) = temp - temp1 - 2.0_rt*temp2 + 1.0_rt;
    }
}

template<int neqs, typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void jcn(amrex::Array1D<amrex::Real, 1, neqs>& x, amrex::Array2D<amrex::Real, 1, neqs, 1, neqs>& fjac, const T& data, int& iflag) {

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
