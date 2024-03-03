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

    RHS::dgefa(A);
    RHS::dgesl(A, b);

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

    std::cout << std::endl;

    std::cout << "the Jacobian mask seen by RHS::is_jacobian_term_used()" << std::endl;

    // now output the Jacobian mask as seen by `is_jacobian_term_used<>()`

    std::cout << std::setw(4) << "*" << " ";

    constexpr_for<1, INT_NEQS+1>([&] (auto i)
    {
        if constexpr (i < INT_NEQS) {
            std::cout << std::setw(4) << short_spec_names_cxx[i-1] << " ";
        } else {
            std::cout << std::setw(4) << "enuc" << " ";
        }
    });

    std::cout << std::endl;

    constexpr_for<1, INT_NEQS+1>([&] (auto irow)
    {

        if constexpr (irow < INT_NEQS) {
            std::cout << std::setw(4) << short_spec_names_cxx[irow-1] << " ";
        } else {
            std::cout << std::setw(4) << "enuc" << " ";
        }

        constexpr_for<1, INT_NEQS+1>([&] (auto jcol)
        {
            std::cout << std::setw(4) << RHS::is_jacobian_term_used<irow, jcol>() << " ";
        });
        std::cout << std::endl;
    });

    std::cout << std::endl;

    // now try to output a Jacobian mask based on the actual Jacobian

    ArrayUtil::MathArray2D<1, INT_NEQS, 1, INT_NEQS> jac;

    burn_t burn_state;
    burn_state.rho = 1.e8;
    burn_state.T = 1.e9;
    for (int i = 0; i < NumSpec; ++i) {
        burn_state.xn[i] = 1.0 / static_cast<amrex::Real>(NumSpec);
    }

    actual_jac(burn_state, jac);

    std::cout << "the Jacobian mask from calling actual_jac" << std::endl;

    std::cout << std::setw(4) << "*" << " ";

    for (int i = 1; i <= INT_NEQS; ++i) {
        if (i < INT_NEQS) {
            std::cout << std::setw(4) << short_spec_names_cxx[i-1] << " ";
        } else {
            std::cout << std::setw(4) << "enuc" << " ";
        }
    }

    std::cout << std::endl;

    for (int irow = 1; irow <= INT_NEQS; ++irow) {
        if (irow < INT_NEQS) {
            std::cout << std::setw(4) << short_spec_names_cxx[irow-1] << " ";
        } else {
            std::cout << std::setw(4) << "enuc" << " ";
        }

        for (int jcol = 1; jcol <= INT_NEQS; ++jcol) {

            if (jac(irow, jcol) != 0.0) {
                std::cout << std::setw(4) << 1 << " ";
            } else {
                std::cout << std::setw(4) << 0 << " ";
            }
        }
        std::cout << std::endl;
    }


}

#endif
