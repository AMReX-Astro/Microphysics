#ifndef eos_query_H
#define eos_query_H

AMREX_INLINE
bool is_filled_eos(const std::string& field) {

    if (field == "p") {
        return true;

    } else if (field == "e") {
        return true;

    } else if (field == "h") {
        return true;

    } else if (field == "s") {
        return true;

    } else if (field == "dpdT") {
        return true;

    } else if (field == "dpdr") {
        return true;

    } else if (field == "dedT") {
        return true;

    } else if (field == "dedr") {
        return true;

    } else if (field == "dhdT") {
        return true;

    } else if (field == "dhdr") {
        return true;

    } else if (field == "dsdT") {
        return true;

    } else if (field == "dsdr") {
        return true;

    } else if (field == "dpde") {
        return true;

    } else if (field == "dpdr_e") {
        return true;

    } else if (field == "cv") {
        return true;

    } else if (field == "cp") {
        return true;

    } else if (field == "xne") {
        return true;

    } else if (field == "xnp") {
        return true;

    } else if (field == "eta") {
        return true;

    } else if (field == "pele") {
        return true;

    } else if (field == "ppos") {
        return true;

    } else if (field == "mu") {
        return true;

    } else if (field == "mu_e") {
        return true;

    } else if (field == "y_e") {
        return true;

    } else if (field == "gam1") {
        return true;

    } else if (field == "cs") {
        return true;

    } else if (field == "abar") {
        return true;

    } else if (field == "zbar") {
        return true;


    } else if (field == "dpdA") {
#ifdef EXTRA_THERMO
        return true;
#else
        return false;
#endif

    } else if (field == "dpdZ") {
#ifdef EXTRA_THERMO
        return true;
#else
        return false;
#endif

    } else if (field == "dedA") {
#ifdef EXTRA_THERMO
        return true;
#else
        return false;
#endif

    } else if (field == "dedZ") {
#ifdef EXTRA_THERMO
        return true;
#else
        return false;
#endif

    }

    return false;
}
#endif
