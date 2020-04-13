#ifndef _array_utilities_H_
#define _array_utilities_H_

#include <AMReX_Array.H>
#include <AMReX_REAL.H>

namespace ArrayUtil
{
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
