// One violation per lint rule; tests/selftest.sh checks that every rule in
// rules/ fires here at least once and that no rule file fails to parse.
#include <cmath>
namespace t {
struct B {};
struct A : B { static constexpr double tab[2] = {1, 2}; };
struct D { double v; };
D operator*(D x, D y) { return D{x.v * y.v}; }
double g;
constexpr double arr[2] = {1, 2};
void g(A& a) { a.v = 1.0; }
void f(double& r, int k)
{
    D u;
    double ok_scalar; Real ok_real; unsigned int ok_u; bool ok_b;
    auto lam = [](int x) { return x; };
    switch (k) { default: break; }
    for (;;) { if (k) return; }
    for (int i : {1, 2}) { (void)i; }
    A a{};
    double w = __builtin_exp(1.0) + std::exp(1.0) + amrex::Math::powi<2>(r) + integrator_rp::rtol_spec;
    throw 1;
}
}
