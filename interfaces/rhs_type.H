#ifndef rhs_type_H
#define rhs_type_H

#include <AMReX_REAL.H>
#include <AMReX_Loop.H>
#include <network_properties.H>
#ifdef SCREENING
#include <screen.H>
#endif
#include <tfactors.H>


namespace RHS
{

struct rhs_t {
    int species_A{-1};
    int species_B{-1};
    int species_C{-1};
    int species_D{-1};
    int species_E{-1};
    int species_F{-1};

    amrex::Real number_A{0.0_rt};
    amrex::Real number_B{0.0_rt};
    amrex::Real number_C{0.0_rt};
    amrex::Real number_D{0.0_rt};
    amrex::Real number_E{0.0_rt};
    amrex::Real number_F{0.0_rt};

    int exponent_A{0};
    int exponent_B{0};
    int exponent_C{0};
    int exponent_D{0};
    int exponent_E{0};
    int exponent_F{0};

    amrex::Real forward_branching_ratio{1.0_rt};
    amrex::Real reverse_branching_ratio{1.0_rt};

    int apply_identical_particle_factor{1};

    int rate_can_be_tabulated{1};

    int screen_forward_reaction{1};
    int screen_reverse_reaction{1};

    int additional_reaction_1{-1};
    int additional_reaction_2{-1};
    int additional_reaction_3{-1};
};

constexpr amrex::Real tab_tlo = 6.0e0_rt;
constexpr amrex::Real tab_thi = 10.0e0_rt;
constexpr int tab_per_decade = 2000;
constexpr int nrattab = static_cast<int>(tab_thi - tab_tlo) * tab_per_decade + 1;
constexpr int tab_imax = static_cast<int>(tab_thi - tab_tlo) * tab_per_decade + 1;
constexpr amrex::Real tab_tstp = (tab_thi - tab_tlo) / static_cast<amrex::Real>(tab_imax - 1);

extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, nrattab> ttab;

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
    amrex::Array1D<amrex::Real, 1, NumSpec> y;
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
