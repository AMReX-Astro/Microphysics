#ifndef rhs_type_H
#define rhs_type_H

#include <AMReX_REAL.H>

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

} // namespace RHS

#endif
