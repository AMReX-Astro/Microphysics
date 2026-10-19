// Runs the native validators of the Loom burn-kernel translation selected at
// build time (-DLOOM_VALIDATE_<MODULE>). See the port's README.md.
#include <iomanip>

#include <extern_parameters.H>
#include <eos.H>
#include <network.H>
#include <burner.H>

#include <loom_validate.H>

#ifdef LOOM_VALIDATE_DUAL
#include <validate_dual.H>
#endif
#ifdef LOOM_VALIDATE_EOS
#include <validate_eos.H>
#endif
#ifdef LOOM_VALIDATE_RATES
#include <validate_rates.H>
#endif
#ifdef LOOM_VALIDATE_NETWORK
#include <validate_network.H>
#endif
#ifdef LOOM_VALIDATE_INTEGRATOR
#include <validate_integrator.H>
#endif
#ifdef LOOM_VALIDATE_KERNEL
#include <validate_kernel.H>
#endif

long loom_validate()
{
    long mismatches = 0;
#ifdef LOOM_VALIDATE_DUAL
    mismatches += validate_dual();
#endif
#ifdef LOOM_VALIDATE_EOS
    mismatches += validate_eos();
#endif
#ifdef LOOM_VALIDATE_RATES
    mismatches += validate_rates();
#endif
#ifdef LOOM_VALIDATE_NETWORK
    mismatches += validate_network();
#endif
#ifdef LOOM_VALIDATE_INTEGRATOR
    mismatches += validate_integrator();
#endif
#ifdef LOOM_VALIDATE_KERNEL
    mismatches += validate_kernel();
#endif
    amrex::Print() << "LOOM_VALIDATE total mismatches=" << mismatches << "\n";
    return mismatches;
}
