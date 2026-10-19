// One violation per lint rule; tests/selftest.sh checks that every rule in
// rules/ fires here at least once and that no rule file fails to parse.
#include <cmath>
namespace t {
struct B {};
struct A : B { static constexpr double tab[2] = {1, 2}; };
double g;
constexpr double arr[2] = {1, 2};
void f(double& r, int k)
{
    double u;
    auto lam = [](int x) { return x; };
    switch (k) { default: break; }
    for (;;) { if (k) return; }
    for (int i : {1, 2}) { (void)i; }
    A a{};
    a.m(1);
    double v = a.operator()(1);
    double w = __builtin_exp(1.0) + std::exp(1.0) + amrex::Math::powi<2>(r) + integrator_rp::rtol_spec;
    throw 1;
}
}
