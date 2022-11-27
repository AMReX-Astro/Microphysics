#ifndef rhs_type_H
#define rhs_type_H

#include <AMReX_REAL.H>
#include <network_properties.H>
#ifdef SCREENING
#include <screen.H>
#endif
#include <tfactors.H>

namespace RHS
{

struct rhs_t {
    int species_A;
    int species_B;
    int species_C;
    int species_D;
    int species_E;
    int species_F;

    int number_A;
    int number_B;
    int number_C;
    int number_D;
    int number_E;
    int number_F;

    int exponent_A;
    int exponent_B;
    int exponent_C;
    int exponent_D;
    int exponent_E;
    int exponent_F;

    amrex::Real forward_branching_ratio;
    amrex::Real reverse_branching_ratio;

    int apply_identical_particle_factor;

    int rate_can_be_tabulated;

    int screen_forward_reaction;
    int screen_reverse_reaction;

    int additional_reaction_1;
    int additional_reaction_2;
    int additional_reaction_3;
};

constexpr amrex::Real tab_tlo = 6.0e0_rt;
constexpr amrex::Real tab_thi = 10.0e0_rt;
constexpr int tab_per_decade = 2000;
constexpr int nrattab = static_cast<int>(tab_thi - tab_tlo) * tab_per_decade + 1;
constexpr int tab_imax = static_cast<int>(tab_thi - tab_tlo) * tab_per_decade + 1;
constexpr amrex::Real tab_tstp = (tab_thi - tab_tlo) / static_cast<Real>(tab_imax - 1);

extern AMREX_GPU_MANAGED Array1D<amrex::Real, 1, nrattab> ttab;

// A struct that contains the terms needed to evaluate a tabulated rate.
struct rate_tab_t
{
    amrex::Real alfa, beta, gama, delt;
    int iat;

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    void initialize (amrex::Real temp) {
        // hash locate
        constexpr int mp = 4;
        iat = static_cast<int>((std::log10(temp) - tab_tlo) / tab_tstp) + 1;
        iat = amrex::max(1, amrex::min(iat - 1, tab_imax - mp + 1));

        // setup the lagrange interpolation coefficients for a cubic
        amrex::Real x  = temp;
        amrex::Real x1 = ttab(iat);
        amrex::Real x2 = ttab(iat+1);
        amrex::Real x3 = ttab(iat+2);
        amrex::Real x4 = ttab(iat+3);
        amrex::Real a  = x - x1;
        amrex::Real b  = x - x2;
        amrex::Real c  = x - x3;
        amrex::Real d  = x - x4;
        amrex::Real e  = x1 - x2;
        amrex::Real f  = x1 - x3;
        amrex::Real g  = x1 - x4;
        amrex::Real h  = x2 - x3;
        amrex::Real p  = x2 - x4;
        amrex::Real q  = x3 - x4;
        alfa =  b * c * d / (e * f * g);
        beta = -a * c * d / (e * h * p);
        gama =  a * b * d / (f * h * q);
        delt = -a * b * c / (g * p * q);
    }
};

struct rhs_state_t
{
    amrex::Real rho;
    tf_t tf;
    rate_tab_t tab;
#ifdef SCREENING
    plasma_state_t pstate;
#endif
    amrex::Real y_e;
    amrex::Real eta;
    Array1D<amrex::Real, 1, NumSpec> y;
};

struct rate_t
{
    amrex::Real fr;
    amrex::Real rr;
    amrex::Real frdt;
    amrex::Real rrdt;
};

} // namespace RHS

#endif
