#ifndef TEST_AUTODIFF_ARRAYS_H
#define TEST_AUTODIFF_ARRAYS_H

#include <algorithm>
#include <iostream>
#include <sstream>
#include <type_traits>

#include <AMReX_Extension.H>
#include <AMReX_REAL.H>
#include <microphysics_autodiff.H>

using namespace amrex::literals;
using microphysics_autodiff::GradArray;
using amrex::Real;

template<typename T, int XLO, int XHI>
AMREX_GPU_HOST_DEVICE
bool operator==(const GradArray<T, XLO, XHI>& s,
                const GradArray<T, XLO, XHI>& p) noexcept
{
    return std::equal(s.begin(), s.end(), p.begin(), p.end());
}

template <int XLO, int XHI>
std::ostream& operator<<(std::ostream& os, const GradArray<Real, XLO, XHI>& arr)
{
    os << "GradArray<Real, " << XLO << ", " << XHI << ">{";
    for (int i = XLO; i <= XHI; ++i) {
        if (i != XLO) {
            os << ", ";
        }
        os << arr(i);
    }
    os << "}";
    return os;
}

template <typename T, typename U>
void assert_equal_impl(T actual, U expected,
        const char* EX, const char* file, int line) {
    if (!(actual == expected)) {
        std::ostringstream oss{};
        oss << "(lhs=" << actual << ", rhs=" << expected << ")";
        amrex::Assert(EX, file, line, oss.str().c_str());
    }
}
#define ASSERT_EQUAL(actual, ...) assert_equal_impl(actual, __VA_ARGS__, #actual " == " #__VA_ARGS__, __FILE__, __LINE__)

template <int XLO, int XHI, int YLO, int YHI>
void check_result_types() {
    using result_t = GradArray<Real, amrex::min(XLO, YLO), amrex::max(XHI, YHI)>;
    GradArray<Real, XLO, XHI> x{1.0_rt};
    GradArray<Real, YLO, YHI> y{2.0_rt};

    static_assert(std::is_same_v<decltype(x + y), result_t>, "incorrect result type for addition");
    static_assert(std::is_same_v<decltype(x - y), result_t>, "incorrect result type for subtraction");
}

void test_constructors() {
    std::cout << "testing constructors..." << std::endl;
    // this macro just removes the parentheses around an expression, which are
    // required if an argument contains commas
#define SINGLE_ARG(...) __VA_ARGS__
    // NOLINTBEGIN(cppcoreguidelines-avoid-do-while)
#define TEST_CONSTRUCTOR(decl, expected) do { \
        std::string decl_str = #decl; \
        std::string expected_str = #expected; \
        /* remove everything not inside parentheses */ \
        decl_str.erase(0, decl_str.find('(') + 1); \
        decl_str.erase(decl_str.rfind(')')); \
        expected_str.erase(0, expected_str.find('(') + 1); \
        expected_str.erase(expected_str.rfind(')')); \
        std::cout << "testing arr == " << expected_str << " for `" << decl_str << ";`" << std::endl; \
        \
        decl; \
        std::initializer_list<Real> l expected; \
        ASSERT_EQUAL(l.size(), arr.size()); \
        auto it = l.begin(); \
        for (int i = arr.lo(); it != l.end(); ++i, ++it) { \
            ASSERT_EQUAL(arr(i), *it); \
        } \
    } while(0)

    TEST_CONSTRUCTOR(SINGLE_ARG(GradArray<Real, 1, 3> arr{}),
                     SINGLE_ARG({0.0_rt, 0.0_rt, 0.0_rt}));
    TEST_CONSTRUCTOR(SINGLE_ARG(GradArray<Real, 1, 3> arr{1.0_rt}),
                     SINGLE_ARG({1.0_rt, 0.0_rt, 0.0_rt}));
    TEST_CONSTRUCTOR(SINGLE_ARG(GradArray<Real, 1, 3> arr{1.0_rt, 2.0_rt}),
                     SINGLE_ARG({1.0_rt, 2.0_rt, 0.0_rt}));
    std::cout << std::endl;
    TEST_CONSTRUCTOR(SINGLE_ARG(GradArray<Real, 1, 3> arr),
                     SINGLE_ARG({0.0_rt, 0.0_rt, 0.0_rt}));
    TEST_CONSTRUCTOR(SINGLE_ARG(GradArray<Real, 1, 3> arr(1.0_rt)),
                     SINGLE_ARG({1.0_rt, 1.0_rt, 1.0_rt}));
    std::cout << std::endl;
    TEST_CONSTRUCTOR(SINGLE_ARG(GradArray<Real, 1, 3> arr = {}),
                     SINGLE_ARG({0.0_rt, 0.0_rt, 0.0_rt}));
    TEST_CONSTRUCTOR(SINGLE_ARG(GradArray<Real, 1, 3> arr = {1.0_rt}),
                     SINGLE_ARG({1.0_rt, 0.0_rt, 0.0_rt}));
    TEST_CONSTRUCTOR(SINGLE_ARG(GradArray<Real, 1, 3> arr = {1.0_rt, 2.0_rt}),
                     SINGLE_ARG({1.0_rt, 2.0_rt, 0.0_rt}));
    std::cout << std::endl;
    TEST_CONSTRUCTOR(SINGLE_ARG(auto arr = autodiff::detail::Zero<GradArray<Real, 1, 3>>()),
                     SINGLE_ARG({0.0_rt, 0.0_rt, 0.0_rt}));
    TEST_CONSTRUCTOR(SINGLE_ARG(auto arr = autodiff::detail::One<GradArray<Real, 1, 3>>()),
                     SINGLE_ARG({1.0_rt, 1.0_rt, 1.0_rt}));
#undef TEST_CONSTRUCTOR
    // NOLINTEND(cppcoreguidelines-avoid-do-while)
#undef SINGLE_ARG
    std::cout << std::endl;
}

void test_autodiff_arrays() {

    check_result_types<1, 1, 1, 3>();
    check_result_types<1, 3, 1, 1>();
    check_result_types<1, 3, 2, 2>();
    check_result_types<1, 1, 3, 3>();

    test_constructors();

    std::cout << "testing GradArray addition & subtraction..." << std::endl;
    {
        ASSERT_EQUAL((GradArray<Real, 1, 1>{1.0_rt} + GradArray<Real, 2, 2>{2.0_rt}), GradArray<Real, 1, 2>{1.0_rt, 2.0_rt});
        ASSERT_EQUAL((GradArray<Real, 1, 1>{1.0_rt} - GradArray<Real, 2, 2>{2.0_rt}), GradArray<Real, 1, 2>{1.0_rt, -2.0_rt});
    }

    std::cout << "testing dual_array narrowing..." << std::endl;
    {
        autodiff::dual_array<1, 3> x = 1.0_rt;
        x.grad = {3.0_rt, 4.0_rt, 5.0_rt};
        amrex::constexpr_for<1, 3>([&] (auto n) {
                auto tmp = microphysics_autodiff::narrow_array<autodiff::dual_array<n, n>>(x);
                ASSERT_EQUAL(tmp.val, 1.0_rt);
                ASSERT_EQUAL(tmp.grad.size(), 1U);
                ASSERT_EQUAL(tmp.grad(n), x.grad(n));
                });
        amrex::constexpr_for<1, 3>([&] (auto n) {
                auto tmp = microphysics_autodiff::narrow_array<n, n>(x);
                ASSERT_EQUAL(tmp.val, 1.0_rt);
                ASSERT_EQUAL(tmp.grad.size(), 1U);
                ASSERT_EQUAL(tmp.grad(n), x.grad(n));
                });
        {
            auto tmp = microphysics_autodiff::narrow_array<autodiff::dual_array<2, 3>>(x);
            ASSERT_EQUAL(tmp.val, 1.0_rt);
            ASSERT_EQUAL(tmp.grad.size(), 2U);
            ASSERT_EQUAL(tmp.grad(2), 4.0_rt);
            ASSERT_EQUAL(tmp.grad(3), 5.0_rt);
        }
    }

    std::cout << "testing make_partial_arrays()..." << std::endl;
    {
        auto [x, y, z] = autodiff::make_partial_arrays<1, 3>(2.0, 3.0, 4.0);
        ASSERT_EQUAL(autodiff::derivative(x), GradArray<Real, 1, 1>{1.0_rt});
        ASSERT_EQUAL(autodiff::derivative(y), GradArray<Real, 2, 2>{1.0_rt});
        ASSERT_EQUAL(autodiff::derivative(z), GradArray<Real, 3, 3>{1.0_rt});

        std::cout << "testing operations with partial dual_arrays..." << std::endl;
        using autodiff::eval;

        const auto x_cubed = admath::powi<3>(x);
        ASSERT_EQUAL(autodiff::val(x_cubed), 8.0_rt);
        ASSERT_EQUAL(autodiff::derivative(x_cubed), GradArray<Real, 1, 1>{12.0_rt});

        const auto x_cubed_y = eval(x_cubed * y);
        ASSERT_EQUAL(autodiff::val(x_cubed_y), 24.0_rt);
        ASSERT_EQUAL(autodiff::derivative(x_cubed_y), GradArray<Real, 1, 2>{36.0_rt, 8.0_rt});

        const auto x_plus_z = eval(x + z);
        ASSERT_EQUAL(autodiff::val(x_plus_z), 6.0_rt);
        ASSERT_EQUAL(autodiff::derivative(x_plus_z), GradArray<Real, 1, 3>{1.0_rt, 0.0_rt, 1.0_rt});

        const auto result = eval(x_cubed_y * admath::sqrt(z));
        ASSERT_EQUAL(autodiff::val(result), 48.0_rt);
        ASSERT_EQUAL(autodiff::derivative(result), GradArray<Real, 1, 3>{72.0_rt, 16.0_rt, 6.0_rt});
    }

    {
        autodiff::dual_array<1, 2> x = 1.0_rt;
        x.grad(1) = 1.0_rt;

        autodiff::dual_array<2, 2> y = 2.0_rt;
        y.grad(2) = 1.0_rt;

        x += 1.0_rt;
        x -= 2.0_rt;
        x *= 3.0_rt;
        x /= 4.0_rt;
        x += y;
        x *= y;
        x -= y;
        x /= y;
    }

    std::cout << "\nall tests passed!" << std::endl;
}

#endif
