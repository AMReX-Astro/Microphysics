#ifndef TEST_LINEAR_ALGEBRA_H
#define TEST_LINEAR_ALGEBRA_H

#include <iostream>
#include <iomanip>

#include <integrator_data.H>
#include <rhs.H>
#include <cmath>

#include <linpack.H>

using namespace amrex::literals;

constexpr Real eps = 1.e-2_rt;

void create_A(RArray2D& A) {

    // create a diagonally-dominant matrix with non-zero
    // elements only where the Jacobian terms are non-zero

    constexpr_for<1, INT_NEQS+1>([&] (auto irow)
    {
        constexpr_for<1, INT_NEQS+1>([&] (auto jcol)
        {

            if constexpr (RHS::is_jacobian_term_used<irow, jcol>()) {
                if constexpr (irow == jcol) {
                    A(irow, jcol) = 1.0 - eps * static_cast<double>(irow);
                } else {
                    A(irow, jcol) = 1.0_rt / (1.0_rt + static_cast<double>(irow + jcol));
                }
            } else {
                A(irow, jcol) = 0.0;
            }
        });
    });


}

RArray1D Ax(RArray2D& A, RArray1D& x) {

    RArray1D b;

    for (int irow = 1; irow <= INT_NEQS; ++irow) {
        b(irow) = 0.0_rt;
        for (int jcol = 1; jcol <= INT_NEQS; ++jcol) {
            b(irow) += A(irow, jcol) * x(jcol);
        }
    }

    return b;
}

void linear_algebra() {

    RArray2D A;

    RArray1D x;
    RArray1D b;

    create_A(A);

    std::cout << "matrix A: " << std::endl;
    std::cout << std::setprecision(2);

    for (int irow = 1; irow <= INT_NEQS; ++irow) {
        for (int jcol = 1; jcol <= INT_NEQS; ++jcol) {
            std::cout << std::setw(5) << A(irow, jcol) << " ";
        }
        std::cout << std::endl;
    }

    std::cout << std::endl;


    for (int jcol = 1; jcol <= INT_NEQS; ++jcol) {
        if (jcol % 2 == 1) {
            x(jcol) = static_cast<double>(jcol);
        } else {
            x(jcol) = 10.0 * static_cast<double>(jcol);
        }
    }

    // now create b = A x

    b = Ax(A, x);

    std::cout << "vectors x and b (before the solve): " << std::endl;

    std::cout << std::setprecision(10);

    for (int jcol = 1; jcol <= INT_NEQS; ++jcol) {
        std::cout << std::setw(15) << x(jcol) << " " << std::setw(15) << b(jcol) << std::endl;
    }

    std::cout << std::endl;

    // solve the linear system with the templated solver

    if (jacobian == 1) {
        constexpr bool allow_skipping_zero_entries{true};
        RHS::dgefa<allow_skipping_zero_entries>(A);
        RHS::dgesl<allow_skipping_zero_entries>(A, b);
    }
    else {
        constexpr bool allow_skipping_zero_entries{false};
        RHS::dgefa<allow_skipping_zero_entries>(A);
        RHS::dgesl<allow_skipping_zero_entries>(A, b);
    }

    std::cout << std::setprecision(16);

    std::cout << "original x and x from the solve (RHS:: solve): " << std::endl;

    for (int jcol = 1; jcol <= INT_NEQS; ++jcol) {
        std::cout << std::setw(20) << x(jcol) << " " << std::setw(20) << b(jcol) << std::endl;
    }

    std::cout << std::endl;


    // now use the linpack.H solver

    // recreate A and b
    create_A(A);
    b = Ax(A, x);

    IArray1D pivot;
    int info;

    constexpr bool allow_pivot{true};

    dgefa<INT_NEQS, allow_pivot>(A, pivot, info);
    dgesl<INT_NEQS, allow_pivot>(A, pivot, b);


    std::cout << "original x and x from the solve (linpack.H solve): " << std::endl;

    for (int jcol = 1; jcol <= INT_NEQS; ++jcol) {
        std::cout << std::setw(20) << x(jcol) << " " << std::setw(20) << b(jcol) << std::endl;
    }

}

#endif
