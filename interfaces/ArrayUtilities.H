#ifndef _array_utilities_H_
#define _array_utilities_H_

#include <AMReX_Array.H>
#include <AMReX_REAL.H>

namespace ArrayUtil
{
    template <int XLEN>
    struct RealArray1D : amrex::Array1D<amrex::Real, 0, XLEN-1>
    {};

    template <int XLEN, int YLEN>
    struct RealArray2D : amrex::Array2D<amrex::Real, 0, XLEN-1,
                                                     0, YLEN-1>
    {};

    template <int XLEN, int YLEN, int ZLEN>
    struct RealArray3D : amrex::Array3D<amrex::Real, 0, XLEN-1,
                                                     0, YLEN-1,
                                                     0, ZLEN-1>
    {};

    namespace Math
    {
        template <int XLO, int XHI>
        AMREX_GPU_HOST_DEVICE AMREX_INLINE
        Real sum(amrex::Array1D<Real, XLO, XHI> const& xarray,
                 const int xlo_s = XLO, const int xhi_s = XHI)
        {
            Real xsum = 0.0_rt;

            for (int i = xlo_s; i <= xhi_s; ++i) {
                xsum += xarray(i);
            }

            return xsum;
        }
    };
};

#endif
