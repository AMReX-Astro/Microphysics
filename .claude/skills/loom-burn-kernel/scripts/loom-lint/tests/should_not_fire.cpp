// Constructs that import at the pin; tests/selftest.sh checks that no rule
// fires here. Add a case whenever a rule is relaxed.
namespace t {
template <int I> struct IC { static constexpr int value = I; constexpr operator int() const { return I; } };
template <int B, int E, class F> void cfor(F&& f) { if constexpr (B < E) { f(IC<B>{}); cfor<B + 1, E>(f); } }
struct S { double x; double e; double T; };
struct Dual {
    double v; double d;
    Dual operator*(Dual b) const { return Dual{v * b.v, d * b.v + v * b.d}; }
    Dual& operator+=(Dual b) { v += b.v; d += b.d; return *this; }
};
struct Pair { double first; double second; };
template <class T> constexpr const T& mx(const T& a, const T& b) { return (a < b) ? b : a; }
constexpr double arr[2] = {1, 2};
static void upd(S& s, double& t) { t = 3.0 * s.e; }
static void set(double x, double& o) { o = 2.0 * x; }
double lookup(int n) { static constexpr double lt[3] = {4.0, 12.0, 16.0}; return lt[n]; }
double f(double in)
{
    double s = 0.0;
    cfor<1, 4>([&](auto n) { constexpr int m = n; s += m * in; });   // l01
    S st;                                                            // d14
    st.x = in; st.T = 0.0;
    double& r = st.x;                                                // r10
    r += 1.0;
    upd(st, st.T);                                                   // r08
    double o;                                                        // r25
    set(in, o);
    Dual a{in, 1.0};
    a += a * a;                                                      // o09
    auto [p, q] = Pair{in, 2.0};                                     // o08
    const double& m = mx(p, q);                                      // r12
    return s + r + o + a.v + m + lookup(1) + arr[0];                 // d22
}
}
