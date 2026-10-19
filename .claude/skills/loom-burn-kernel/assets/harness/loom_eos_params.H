// Host-side EOS parameter and table setup for the Loom burn kernel:
// fill_eos_params and fill_eos_tables, for the translated helmholtz EOS
// (eos.H). Network-independent: copy as is to
// unit_test/test_react/loom_eos_params.H. Included by <loom_params.H> (the
// launcher and the burn validators) and by validate_eos.H, which therefore
// builds as soon as the kernel's eos.H and network_data.H exist, before
// rates.H, network.H or loom_params.H.
#ifndef LOOM_EOS_PARAMS_H
#define LOOM_EOS_PARAMS_H

#include <extern_parameters.H>
#include <eos.H>          // stock helmholtz:: / EOSData:: globals

// The kernel directory is LAST on the include path (its headers share names
// with stock ones), so its eos.H is reached by path from this directory.
// LOOM_KERNEL_DIR_PATH is the port's kernel directory relative to test_react,
// unquoted (-DLOOM_KERNEL_DIR_PATH=loom/<network>_<integrator>, set by the
// GNUmakefile from LOOM_KERNEL_DIR; test_react.patch). LOOM_KERNEL_PATH(f)
// makes the include string "loom/<network>_<integrator>/f".
#ifndef LOOM_KERNEL_PATH
#ifndef LOOM_KERNEL_DIR_PATH
#error "LOOM_KERNEL_DIR_PATH undefined: build with LOOM_KERNEL_DIR=loom/<network>_<integrator> (assets/harness/test_react.patch)"
#endif
#define LOOM_KERNEL_STR_(x) #x
#define LOOM_KERNEL_STR(x) LOOM_KERNEL_STR_(x)
#define LOOM_KERNEL_PATH(f) LOOM_KERNEL_STR(LOOM_KERNEL_DIR_PATH/f)
#endif
#include LOOM_KERNEL_PATH(eos.H)

// Read the stock EOS runtime parameters (after eos_init()).
inline void fill_eos_params(mpl::EosParams* p)
{
    p->eos_prad_limiter_rho_c = eos_rp::prad_limiter_rho_c;
    p->eos_prad_limiter_delta_rho = eos_rp::prad_limiter_delta_rho;
    p->helmholtz_do_coulomb = helmholtz::do_coulomb;
    p->helmholtz_input_is_constant = helmholtz::input_is_constant;
    p->helmholtz_ttol = helmholtz::ttol;
    p->helmholtz_dtol = helmholtz::dtol;
    p->EOSData_mintemp = EOSData::mintemp;
    p->EOSData_maxtemp = EOSData::maxtemp;
    p->EOSData_mindens = EOSData::mindens;
    p->EOSData_maxdens = EOSData::maxdens;
    p->EOSData_mine = EOSData::mine;
    p->EOSData_maxe = EOSData::maxe;
}

// Point at the stock helmholtz tables. On the host these are ordinary
// globals; in a HIP build (AMREX_USE_GPU) they are AMREX_GPU_MANAGED, so the
// same addresses are valid on the device (managed memory), and the launcher
// passes them straight to the kernel. Without managed memory a launcher
// must copy each array (same layout) to device memory instead.
inline void fill_eos_tables(mpl::EosTables* t)
{
    t->f = &helmholtz::f[0][0][0];
    t->dpdf = &helmholtz::dpdf[0][0][0];
    t->ef = &helmholtz::ef[0][0][0];
    t->xf = &helmholtz::xf[0][0][0];
    t->d = helmholtz::d;
    t->t = helmholtz::t;
    t->dt_sav = helmholtz::dt_sav;
    t->dt2_sav = helmholtz::dt2_sav;
    t->dti_sav = helmholtz::dti_sav;
    t->dt2i_sav = helmholtz::dt2i_sav;
    t->dd_sav = helmholtz::dd_sav;
    t->dd2_sav = helmholtz::dd2_sav;
    t->ddi_sav = helmholtz::ddi_sav;
    t->dd2i_sav = helmholtz::dd2i_sav;
}

#endif
