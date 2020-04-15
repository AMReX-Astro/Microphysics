#ifndef _actual_linear_solver_H_
#define _actual_linear_solver_H_

#include <AMReX_Array.H>
#include <AMReX_REAL.H>

using namespace amrex;

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
actual_solve(Array2D<Real, 1, 9, 1, 9> const& A,
             Array1D<Real, 1, 9>& b)
{
    Real scratch_0 = 1.0/A(1,1);
    Real scratch_1 = b(1)*scratch_0;
    Real scratch_2 = A(1,3)*scratch_0;
    Real scratch_3 = A(1,2)*scratch_0;
    Real scratch_4 = 1.0/(-A(2,1)*scratch_3 + A(2,2));
    Real scratch_5 = scratch_4*(-A(2,1)*scratch_2 + A(2,3));
    Real scratch_6 = scratch_3*scratch_5;
    Real scratch_7 = scratch_2 - scratch_6;
    Real scratch_8 = -A(3,1)*scratch_3 + A(3,2);
    Real scratch_9 = 1.0/(-A(3,1)*scratch_2 + A(3,3) - scratch_5*scratch_8);
    Real scratch_10 = scratch_4*(-A(2,1)*scratch_1 + b(2));
    Real scratch_11 = scratch_9*(-A(3,1)*scratch_1 + b(3) - scratch_10*scratch_8);
    Real scratch_12 = A(1,4)*scratch_0;
    Real scratch_13 = A(2,1)*scratch_4;
    Real scratch_14 = scratch_12*scratch_13;
    Real scratch_15 = scratch_9*(-A(3,1)*scratch_12 + A(3,4) + scratch_14*scratch_8);
    Real scratch_16 = A(1,2)*scratch_13/pow(A(1,1), 2);
    Real scratch_17 = A(1,4)*scratch_16;
    Real scratch_18 = scratch_12 - scratch_15*scratch_7 + scratch_17;
    Real scratch_19 = -A(4,1)*scratch_3 + A(4,2);
    Real scratch_20 = -A(4,1)*scratch_2 + A(4,3) - scratch_19*scratch_5;
    Real scratch_21 = 1.0/(-A(4,1)*scratch_12 + A(4,4) + scratch_14*scratch_19 - scratch_15*scratch_20);
    Real scratch_22 = scratch_21*(-A(4,1)*scratch_1 + b(4) - scratch_10*scratch_19 - scratch_11*scratch_20);
    Real scratch_23 = A(1,5)*scratch_0;
    Real scratch_24 = scratch_13*scratch_23;
    Real scratch_25 = scratch_9*(-A(3,1)*scratch_23 + scratch_24*scratch_8);
    Real scratch_26 = scratch_21*(-A(4,1)*scratch_23 + A(4,5) + scratch_19*scratch_24 - scratch_20*scratch_25);
    Real scratch_27 = A(1,5)*scratch_16;
    Real scratch_28 = -scratch_18*scratch_26 + scratch_23 - scratch_25*scratch_7 + scratch_27;
    Real scratch_29 = -A(5,1)*scratch_3 + A(5,2);
    Real scratch_30 = -A(5,1)*scratch_2 + A(5,3) - scratch_29*scratch_5;
    Real scratch_31 = -A(5,1)*scratch_12 + A(5,4) + scratch_14*scratch_29 - scratch_15*scratch_30;
    Real scratch_32 = 1.0/(-A(5,1)*scratch_23 + A(5,5) + scratch_24*scratch_29 - scratch_25*scratch_30 - scratch_26*scratch_31);
    Real scratch_33 = scratch_32*(-A(5,1)*scratch_1 + b(5) - scratch_10*scratch_29 - scratch_11*scratch_30 - scratch_22*scratch_31);
    Real scratch_34 = A(1,6)*scratch_0;
    Real scratch_35 = scratch_13*scratch_34;
    Real scratch_36 = scratch_9*(-A(3,1)*scratch_34 + scratch_35*scratch_8);
    Real scratch_37 = scratch_21*(-A(4,1)*scratch_34 + scratch_19*scratch_35 - scratch_20*scratch_36);
    Real scratch_38 = scratch_32*(-A(5,1)*scratch_34 + A(5,6) + scratch_29*scratch_35 - scratch_30*scratch_36 - scratch_31*scratch_37);
    Real scratch_39 = A(1,6)*scratch_16;
    Real scratch_40 = -scratch_18*scratch_37 - scratch_28*scratch_38 + scratch_34 - scratch_36*scratch_7 + scratch_39;
    Real scratch_41 = -A(6,1)*scratch_3 + A(6,2);
    Real scratch_42 = -A(6,1)*scratch_2 + A(6,3) - scratch_41*scratch_5;
    Real scratch_43 = -A(6,1)*scratch_12 + scratch_14*scratch_41 - scratch_15*scratch_42;
    Real scratch_44 = -A(6,1)*scratch_23 + A(6,5) + scratch_24*scratch_41 - scratch_25*scratch_42 - scratch_26*scratch_43;
    Real scratch_45 = 1.0/(-A(6,1)*scratch_34 + A(6,6) + scratch_35*scratch_41 - scratch_36*scratch_42 - scratch_37*scratch_43 - scratch_38*scratch_44);
    Real scratch_46 = scratch_45*(-A(6,1)*scratch_1 + b(6) - scratch_10*scratch_41 - scratch_11*scratch_42 - scratch_22*scratch_43 - scratch_33*scratch_44);
    Real scratch_47 = A(1,7)*scratch_0;
    Real scratch_48 = scratch_13*scratch_47;
    Real scratch_49 = scratch_9*(-A(3,1)*scratch_47 + scratch_48*scratch_8);
    Real scratch_50 = scratch_21*(-A(4,1)*scratch_47 + scratch_19*scratch_48 - scratch_20*scratch_49);
    Real scratch_51 = scratch_32*(-A(5,1)*scratch_47 + scratch_29*scratch_48 - scratch_30*scratch_49 - scratch_31*scratch_50);
    Real scratch_52 = scratch_45*(-A(6,1)*scratch_47 + A(6,7) + scratch_41*scratch_48 - scratch_42*scratch_49 - scratch_43*scratch_50 - scratch_44*scratch_51);
    Real scratch_53 = A(1,7)*scratch_16;
    Real scratch_54 = -scratch_18*scratch_50 - scratch_28*scratch_51 - scratch_40*scratch_52 + scratch_47 - scratch_49*scratch_7 + scratch_53;
    Real scratch_55 = -A(7,1)*scratch_2 + A(7,1)*scratch_6;
    Real scratch_56 = -A(7,1)*scratch_12 - A(7,1)*scratch_17 - scratch_15*scratch_55;
    Real scratch_57 = -A(7,1)*scratch_23 - A(7,1)*scratch_27 - scratch_25*scratch_55 - scratch_26*scratch_56;
    Real scratch_58 = -A(7,1)*scratch_34 - A(7,1)*scratch_39 + A(7,6) - scratch_36*scratch_55 - scratch_37*scratch_56 - scratch_38*scratch_57;
    Real scratch_59 = 1.0/(-A(7,1)*scratch_47 - A(7,1)*scratch_53 + A(7,7) - scratch_49*scratch_55 - scratch_50*scratch_56 - scratch_51*scratch_57 - scratch_52*scratch_58);
    Real scratch_60 = scratch_10*scratch_3;
    Real scratch_61 = scratch_59*(-A(7,1)*scratch_1 + A(7,1)*scratch_60 + b(7) - scratch_11*scratch_55 - scratch_22*scratch_56 - scratch_33*scratch_57 - scratch_46*scratch_58);
    Real scratch_62 = A(1,8)*scratch_0;
    Real scratch_63 = scratch_4*(-A(2,1)*scratch_62 + A(2,8));
    Real scratch_64 = scratch_9*(-A(3,1)*scratch_62 + A(3,8) - scratch_63*scratch_8);
    Real scratch_65 = scratch_21*(-A(4,1)*scratch_62 + A(4,8) - scratch_19*scratch_63 - scratch_20*scratch_64);
    Real scratch_66 = scratch_32*(-A(5,1)*scratch_62 + A(5,8) - scratch_29*scratch_63 - scratch_30*scratch_64 - scratch_31*scratch_65);
    Real scratch_67 = scratch_45*(-A(6,1)*scratch_62 + A(6,8) - scratch_41*scratch_63 - scratch_42*scratch_64 - scratch_43*scratch_65 - scratch_44*scratch_66);
    Real scratch_68 = scratch_3*scratch_63;
    Real scratch_69 = scratch_59*(-A(7,1)*scratch_62 + A(7,1)*scratch_68 + A(7,8) - scratch_55*scratch_64 - scratch_56*scratch_65 - scratch_57*scratch_66 - scratch_58*scratch_67);
    Real scratch_70 = -A(8,1)*scratch_3 + A(8,2);
    Real scratch_71 = -A(8,1)*scratch_2 + A(8,3) - scratch_5*scratch_70;
    Real scratch_72 = -A(8,1)*scratch_12 + A(8,4) + scratch_14*scratch_70 - scratch_15*scratch_71;
    Real scratch_73 = -A(8,1)*scratch_23 + A(8,5) + scratch_24*scratch_70 - scratch_25*scratch_71 - scratch_26*scratch_72;
    Real scratch_74 = -A(8,1)*scratch_34 + A(8,6) + scratch_35*scratch_70 - scratch_36*scratch_71 - scratch_37*scratch_72 - scratch_38*scratch_73;
    Real scratch_75 = -A(8,1)*scratch_47 + A(8,7) + scratch_48*scratch_70 - scratch_49*scratch_71 - scratch_50*scratch_72 - scratch_51*scratch_73 - scratch_52*scratch_74;
    Real scratch_76 = (-A(8,1)*scratch_1 + b(8) - scratch_10*scratch_70 - scratch_11*scratch_71 - scratch_22*scratch_72 - scratch_33*scratch_73 - scratch_46*scratch_74 - scratch_61*scratch_75)/(-A(8,1)*scratch_62 + A(8,8) - scratch_63*scratch_70 - scratch_64*scratch_71 - scratch_65*scratch_72 - scratch_66*scratch_73 - scratch_67*scratch_74 - scratch_69*scratch_75);
    Real scratch_77 = -scratch_14 - scratch_15*scratch_5;
    Real scratch_78 = -scratch_24 - scratch_25*scratch_5 - scratch_26*scratch_77;
    Real scratch_79 = -scratch_35 - scratch_36*scratch_5 - scratch_37*scratch_77 - scratch_38*scratch_78;
    Real scratch_80 = -scratch_48 - scratch_49*scratch_5 - scratch_50*scratch_77 - scratch_51*scratch_78 - scratch_52*scratch_79;
    Real scratch_81 = -scratch_15*scratch_26 + scratch_25;
    Real scratch_82 = -scratch_15*scratch_37 + scratch_36 - scratch_38*scratch_81;
    Real scratch_83 = -scratch_15*scratch_50 + scratch_49 - scratch_51*scratch_81 - scratch_52*scratch_82;
    Real scratch_84 = -scratch_26*scratch_38 + scratch_37;
    Real scratch_85 = -scratch_26*scratch_51 + scratch_50 - scratch_52*scratch_84;
    Real scratch_86 = -scratch_38*scratch_52 + scratch_51;
    Real scratch_87 = -A(9,1)*scratch_3 + A(9,2);
    Real scratch_88 = -A(9,1)*scratch_2 + A(9,3) - scratch_5*scratch_87;
    Real scratch_89 = -A(9,1)*scratch_12 + A(9,4) + scratch_14*scratch_87 - scratch_15*scratch_88;
    Real scratch_90 = -A(9,1)*scratch_23 + A(9,5) + scratch_24*scratch_87 - scratch_25*scratch_88 - scratch_26*scratch_89;
    Real scratch_91 = -A(9,1)*scratch_34 + A(9,6) + scratch_35*scratch_87 - scratch_36*scratch_88 - scratch_37*scratch_89 - scratch_38*scratch_90;
    Real scratch_92 = -A(9,1)*scratch_47 + A(9,7) + scratch_48*scratch_87 - scratch_49*scratch_88 - scratch_50*scratch_89 - scratch_51*scratch_90 - scratch_52*scratch_91;

    Array1D<Real, 1, 9> x;

    x(1) = scratch_1 - scratch_11*scratch_7 - scratch_18*scratch_22 - scratch_28*scratch_33 - scratch_40*scratch_46 - scratch_54*scratch_61 - scratch_60 - scratch_76*(-scratch_18*scratch_65 - scratch_28*scratch_66 - scratch_40*scratch_67 - scratch_54*scratch_69 + scratch_62 - scratch_64*scratch_7 - scratch_68);
    x(2) = scratch_10 - scratch_11*scratch_5 - scratch_22*scratch_77 - scratch_33*scratch_78 - scratch_46*scratch_79 - scratch_61*scratch_80 - scratch_76*(-scratch_5*scratch_64 + scratch_63 - scratch_65*scratch_77 - scratch_66*scratch_78 - scratch_67*scratch_79 - scratch_69*scratch_80);
    x(3) = scratch_11 - scratch_15*scratch_22 - scratch_33*scratch_81 - scratch_46*scratch_82 - scratch_61*scratch_83 - scratch_76*(-scratch_15*scratch_65 + scratch_64 - scratch_66*scratch_81 - scratch_67*scratch_82 - scratch_69*scratch_83);
    x(4) = scratch_22 - scratch_26*scratch_33 - scratch_46*scratch_84 - scratch_61*scratch_85 - scratch_76*(-scratch_26*scratch_66 + scratch_65 - scratch_67*scratch_84 - scratch_69*scratch_85);
    x(5) = scratch_33 - scratch_38*scratch_46 - scratch_61*scratch_86 - scratch_76*(-scratch_38*scratch_67 + scratch_66 - scratch_69*scratch_86);
    x(6) = scratch_46 - scratch_52*scratch_61 - scratch_76*(-scratch_52*scratch_69 + scratch_67);
    x(7) = scratch_61 - scratch_69*scratch_76;
    x(8) = scratch_76;
    x(9) = (-A(9,1)*scratch_1 + b(9) - scratch_10*scratch_87 - scratch_11*scratch_88 - scratch_22*scratch_89 - scratch_33*scratch_90 - scratch_46*scratch_91 - scratch_61*scratch_92 - scratch_76*(-A(9,1)*scratch_62 + A(9,8) - scratch_63*scratch_87 - scratch_64*scratch_88 - scratch_65*scratch_89 - scratch_66*scratch_90 - scratch_67*scratch_91 - scratch_69*scratch_92))/A(9,9);

    for (int i = 1; i <= 9; ++i) {
        b(i) = x(i);
    }
}

#endif
