#ifndef LINPACK_H
#define LINPACK_H

#include <AMReX_REAL.H>
#include <AMReX_Array.H>

#include <ArrayUtilities.H>

namespace linpack_impl {
    // Implementation of "constexpr for" based on
    // https://artificial-mind.net/blog/2020/10/31/constexpr-for
    template<auto I, auto N, class F>
    constexpr void linalg_for (F&& f)
    {
        if constexpr (I < N) {
            f(std::integral_constant<decltype(I), I>());
            linalg_for<I+1, N>(f);
        }
    }
}

template <int num_eqs>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void dgesl (RArray2D& a1, IArray1D& pivot1, RArray1D& b1)
{
    auto const& a = reinterpret_cast<Array2D<Real, 0, num_eqs-1, 0, num_eqs-1>&>(a1);
    auto const& pivot = reinterpret_cast<Array1D<int, 0, num_eqs-1>&>(pivot1);
    auto& b = reinterpret_cast<Array1D<Real, 0, num_eqs-1>&>(b1);

    constexpr int nm1 = num_eqs - 1;

    // solve a * x = b
    // first solve l * y = b
    if constexpr (nm1 >= 1) {
        linpack_impl::linalg_for<0, nm1>([&] (auto k)
        {
            int l = pivot(k) - 1;
            Real t = b(l);
            if (l != k) {
                b(l) = b(k);
                b(k) = t;
            }

            linpack_impl::linalg_for<k+1, num_eqs>([&] (auto j)
            {
                b(j) += t * a(j,k);
            });
        });
    }

    // now solve u * x = y
    linpack_impl::linalg_for<0, num_eqs>([&] (auto kb)
    {
        constexpr int k = num_eqs - kb - 1;
        b(k) = b(k) / a(k,k);
        Real t = -b(k);
        linpack_impl::linalg_for<0, k>([&] (auto j)
        {
            b(j) += t * a(j,k);
        });
    });
}



template <int num_eqs>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void dgefa (RArray2D& a1, IArray1D& pivot1, int& info)
{
    auto& a = reinterpret_cast<Array2D<Real, 0, num_eqs-1, 0, num_eqs-1>&>(a1);
    auto& pivot = reinterpret_cast<Array1D<int, 0, num_eqs-1>&>(pivot1);

    // dgefa factors a matrix by gaussian elimination.
    // a is returned in the form a = l * u where
    // l is a product of permutation and unit lower
    // triangular matrices and u is upper triangular.

    // gaussian elimination with partial pivoting

    info = 0;
    constexpr int nm1 = num_eqs - 1;

    Real t;

    if constexpr (nm1 >= 1) {

        linpack_impl::linalg_for<0, nm1>([&] (auto k)
        {
            // find l = pivot index
            int l = k;
            Real dmax = std::abs(a(k,k));
            linpack_impl::linalg_for<k+1, num_eqs>([&] (auto i)
            {
                if (std::abs(a(i,k)) > dmax) {
                    l = i;
                    dmax = std::abs(a(i,k));
                }
            });

            pivot(k) = l + 1;

            // zero pivot implies this column already triangularized
            if (a(l,k) != 0.0e0_rt) {

                // interchange if necessary
                if (l != k) {
                    t = a(l,k);
                    a(l,k) = a(k,k);
                    a(k,k) = t;
                }

                // compute multipliers
                t = -1.0e0_rt / a(k,k);
                linpack_impl::linalg_for<k+1, num_eqs>([&] (auto j)
                {
                    a(j,k) *= t;
                });

                constexpr int start = k+1;

                // row elimination with column indexing
                linpack_impl::linalg_for<start, num_eqs>([&] (auto j)
                {
                    t = a(l,j);
                    if (l != k) {
                        a(l,j) = a(k,j);
                        a(k,j) = t;
                    }
                    linpack_impl::linalg_for<start, num_eqs>([&] (auto i)
                    {
                        a(i,j) += t * a(i,k);
                    });
                });
            }
            else {

                info = k+1;

            }

        });

    }

    pivot(num_eqs-1) = num_eqs;

    if (a(num_eqs-1,num_eqs-1) == 0.0e0_rt) {
        info = num_eqs;
    }
}

#endif
