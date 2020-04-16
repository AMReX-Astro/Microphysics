#ifndef _array_utilities_H_
#define _array_utilities_H_

#include <AMReX_Array.H>
#include <AMReX_REAL.H>

namespace ArrayUtil
{
    using namespace amrex;

    template <int XLO, int XHI>
    struct MathArray1D
    {
        AMREX_GPU_HOST_DEVICE AMREX_INLINE
        void zero()
        {
            for (int i = 0; i < (XHI-XLO+1); ++i) {
                arr[i] = 0.0_rt;
            }
        }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void operator= (const MathArray1D<XLO, XHI>& other) noexcept {
            for (int i = 0; i < (XHI-XLO+1); ++i) {
                arr[i] = other.arr[i];
            }
        }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        const Real& operator() (int i) const noexcept {
            return arr[i-XLO];
        }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        Real& operator() (int i) noexcept {
            return arr[i-XLO];
        }

        Real arr[(XHI-XLO+1)];
    };

    template <int XLO, int XHI, int YLO, int YHI>
    struct MathArray2D
    {
        AMREX_GPU_HOST_DEVICE AMREX_INLINE
        void zero()
        {
            for (int i = 0; i < (YHI-YLO+1)*(XHI-XLO+1); ++i) {
                arr[i] = 0.0_rt;
            }
        }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void mul (const Real x) noexcept {
            for (int i = 0; i < (YHI-YLO+1)*(XHI-XLO+1); ++i) {
                arr[i] *= x;
            }
        }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void set (const int i, const int j, const Real x) noexcept {
            arr[i+j*(XHI-XLO+1)-(YLO*(XHI-XLO+1)+XLO)] = x;
        }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void add (const int i, const int j, const Real x) noexcept {
            arr[i+j*(XHI-XLO+1)-(YLO*(XHI-XLO+1)+XLO)] += x;
        }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void mul (const int i, const int j, const Real x) noexcept {
            arr[i+j*(XHI-XLO+1)-(YLO*(XHI-XLO+1)+XLO)] *= x;
        }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        Real get (const int i, const int j) const noexcept {
            return arr[i+j*(XHI-XLO+1)-(YLO*(XHI-XLO+1)+XLO)];
        }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void add_identity () noexcept {
            for (int i = XLO; i <= XHI; ++i) {
                arr[i+i*(XHI-XLO+1)-(YLO*(XHI-XLO+1)+XLO)] += 1.0_rt;
            }
        }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void operator= (const MathArray2D<XLO, XHI, YLO, YHI>& other) noexcept {
            for (int i = 0; i < (YHI-YLO+1)*(XHI-XLO+1); ++i) {
                arr[i] = other.arr[i];
            }
        }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        const Real& operator() (int i, int j) const noexcept {
            return arr[i+j*(XHI-XLO+1)-(YLO*(XHI-XLO+1)+XLO)];
        }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        Real& operator() (int i, int j) noexcept {
            return arr[i+j*(XHI-XLO+1)-(YLO*(XHI-XLO+1)+XLO)];
        }

        Real arr[(XHI-XLO+1)*(YHI-YLO+1)];
    };

    namespace Math
    {
        template <int XLO, int XHI>
        AMREX_GPU_HOST_DEVICE AMREX_INLINE
        amrex::Real sum(amrex::Array1D<amrex::Real, XLO, XHI> const& xarray,
                 const int xlo_s = XLO, const int xhi_s = XHI)
        {
            amrex::Real xsum = 0.0_rt;

            for (int i = xlo_s; i <= xhi_s; ++i) {
                xsum += xarray(i);
            }

            return xsum;
        }
    };
};

#endif
