#ifndef LOOM_VALIDATE_H
#define LOOM_VALIDATE_H

// Native bit-exact validation of the Loom burn-kernel translation
// (see the port's README.md). Each module's validate_<module>.H defines
// `int validate_<module>()` returning its mismatch count, using these helpers.

#include <cstdint>
#include <cstring>
#include <iomanip>
#include <string>
#include <AMReX_Print.H>

struct LoomCheck {
    std::string module;
    long checks{0};
    long mismatches{0};
    int reported{0};
};

// Compare two doubles bit for bit; any NaN matches any NaN (payload and sign ignored).
inline bool loom_check_bits(LoomCheck& c, const char* what, double stock, double loom)
{
    std::uint64_t a = 0;
    std::uint64_t b = 0;
    std::memcpy(&a, &stock, sizeof a);
    std::memcpy(&b, &loom, sizeof b);
    ++c.checks;
    if (a == b) return true;
    // any NaN matches any NaN: sign and payload bits are not compared
    if (stock != stock && loom != loom) return true;
    ++c.mismatches;
    if (c.reported < 20) {
        ++c.reported;
        amrex::Print() << std::setprecision(17) << "LOOM_MISMATCH " << c.module << " " << what
                       << " stock=" << stock << " loom=" << loom << "\n";
    }
    return false;
}

inline bool loom_check_int(LoomCheck& c, const char* what, long stock, long loom)
{
    ++c.checks;
    if (stock == loom) return true;
    ++c.mismatches;
    if (c.reported < 20) {
        ++c.reported;
        amrex::Print() << "LOOM_MISMATCH " << c.module << " " << what
                       << " stock=" << stock << " loom=" << loom << "\n";
    }
    return false;
}

inline long loom_report(const LoomCheck& c)
{
    amrex::Print() << "LOOM_VALIDATE " << c.module << " checks=" << c.checks
                   << " mismatches=" << c.mismatches << "\n";
    return c.mismatches;
}

long loom_validate();

#endif
