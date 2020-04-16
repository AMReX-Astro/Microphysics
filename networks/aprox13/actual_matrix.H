#ifndef _actual_matrix_H_
#define _actual_matrix_H_

#include <cmath>

#include <AMReX_Array.H>
#include <AMReX_REAL.H>

#include <network_properties.H>

using namespace amrex;

struct SparseMatrix
{
    static constexpr int flatten(const int i, const int j) {
        using namespace Species;

        if (i==He4) {
            if (j==He4)
                return 0;
            else if (j==C12)
                return 1;
            else if (j==O16)
                return 2;
            else if (j==Ne20)
                return 3;
            else if (j==Mg24)
                return 4;
            else if (j==Si28)
                return 5;
            else if (j==S32)
                return 6;
            else if (j==Ar36)
                return 7;
            else if (j==Ca40)
                return 8;
            else if (j==Ti44)
                return 9;
            else if (j==Cr48)
                return 10;
            else if (j==Fe52)
                return 11;
            else if (j==Ni56)
                return 12;
            else if (j==net_itemp)
                return 13;
        } else if (i==C12) {
            if (j==He4)
                return 14;
            else if (j==C12)
                return 15;
            else if (j==O16)
                return 16;
            else if (j==net_itemp)
                return 17;
        } else if (i==O16) {
            if (j==He4)
                return 18;
            else if (j==C12)
                return 19;
            else if (j==O16)
                return 20;
            else if (j==Ne20)
                return 21;
            else if (j==net_itemp)
                return 22;
        } else if (i==Ne20) {
            if (j==He4)
                return 23;
            else if (j==C12)
                return 24;
            else if (j==O16)
                return 25;
            else if (j==Ne20)
                return 26;
            else if (j==Mg24)
                return 27;
            else if (j==net_itemp)
                return 28;
        } else if (i==Mg24) {
            if (j==He4)
                return 29;
            else if (j==C12)
                return 30;
            else if (j==O16)
                return 31;
            else if (j==Ne20)
                return 32;
            else if (j==Mg24)
                return 33;
            else if (j==Si28)
                return 34;
            else if (j==net_itemp)
                return 35;
        } else if (i==Si28) {
            if (j==He4)
                return 36;
            else if (j==C12)
                return 37;
            else if (j==O16)
                return 38;
            else if (j==Mg24)
                return 39;
            else if (j==Si28)
                return 40;
            else if (j==S32)
                return 41;
            else if (j==net_itemp)
                return 42;
        } else if (i==S32) {
            if (j==He4)
                return 43;
            else if (j==O16)
                return 44;
            else if (j==Si28)
                return 45;
            else if (j==S32)
                return 46;
            else if (j==Ar36)
                return 47;
            else if (j==net_itemp)
                return 48;
        } else if (i==Ar36) {
            if (j==He4)
                return 49;
            else if (j==S32)
                return 50;
            else if (j==Ar36)
                return 51;
            else if (j==Ca40)
                return 52;
            else if (j==net_itemp)
                return 53;
        } else if (i==Ca40) {
            if (j==He4)
                return 54;
            else if (j==Ar36)
                return 55;
            else if (j==Ca40)
                return 56;
            else if (j==Ti44)
                return 57;
            else if (j==net_itemp)
                return 58;
        } else if (i==Ti44) {
            if (j==He4)
                return 59;
            else if (j==Ca40)
                return 60;
            else if (j==Ti44)
                return 61;
            else if (j==Cr48)
                return 62;
            else if (j==net_itemp)
                return 63;
        } else if (i==Cr48) {
            if (j==He4)
                return 64;
            else if (j==Ti44)
                return 65;
            else if (j==Cr48)
                return 66;
            else if (j==Fe52)
                return 67;
            else if (j==net_itemp)
                return 68;
        } else if (i==Fe52) {
            if (j==He4)
                return 69;
            else if (j==Cr48)
                return 70;
            else if (j==Fe52)
                return 71;
            else if (j==Ni56)
                return 72;
            else if (j==net_itemp)
                return 73;
        } else if (i==Ni56) {
            if (j==He4)
                return 74;
            else if (j==Fe52)
                return 75;
            else if (j==Ni56)
                return 76;
            else if (j==net_itemp)
                return 77;
        } else if (i==net_itemp) {
            if (j==He4)
                return 78;
            else if (j==C12)
                return 79;
            else if (j==O16)
                return 80;
            else if (j==Ne20)
                return 81;
            else if (j==Mg24)
                return 82;
            else if (j==Si28)
                return 83;
            else if (j==S32)
                return 84;
            else if (j==Ar36)
                return 85;
            else if (j==Ca40)
                return 86;
            else if (j==Ti44)
                return 87;
            else if (j==Cr48)
                return 88;
            else if (j==Fe52)
                return 89;
            else if (j==Ni56)
                return 90;
            else if (j==net_itemp)
                return 91;
        } else if (i==net_ienuc) {
            if (j==He4)
                return 92;
            else if (j==C12)
                return 93;
            else if (j==O16)
                return 94;
            else if (j==Ne20)
                return 95;
            else if (j==Mg24)
                return 96;
            else if (j==Si28)
                return 97;
            else if (j==S32)
                return 98;
            else if (j==Ar36)
                return 99;
            else if (j==Ca40)
                return 100;
            else if (j==Ti44)
                return 101;
            else if (j==Cr48)
                return 102;
            else if (j==Fe52)
                return 103;
            else if (j==Ni56)
                return 104;
            else if (j==net_itemp)
                return 105;
            else if (j==net_ienuc)
                return 106;
        }

        // If we got here, then out-of-bounds
        // indices i,j were passed so return an invalid index.
        return -1;
    }


    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void zero () noexcept {
        for (int i = 0; i < 107; ++i)
            arr[i] = 0.0_rt;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void mul (const Real x) noexcept {
        for (int i = 0; i < 107; ++i)
            arr[i] *= x;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void mul (const int i, const int j, const Real x) noexcept {
        const int k = SparseMatrix::flatten(i,j);
        // Don't do asserts here because we just skip invalid indices
        if (k >= 0 && k < 107) {
            arr[k] *= x;
        }
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void set (const int i, const int j, const Real x) noexcept {
        const int k = SparseMatrix::flatten(i,j);
        // Don't do asserts here because we just skip invalid indices
        if (k >= 0 && k < 107) {
            arr[k] = x;
        }
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real get (const int i, const int j) noexcept {
        const int k = SparseMatrix::flatten(i,j);
        // Don't do asserts here because we return 0 for invalid indices
        if (k >= 0 && k < 107) {
            return arr[k];
        } else {
            return 0.0_rt;
        }
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void add (const int i, const int j, const Real x) noexcept {
        const int k = SparseMatrix::flatten(i,j);
        // Don't do asserts here because we just skip invalid indices
        if (k >= 0 && k < 107) {
            arr[k] += x;
        }
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void add_identity () noexcept {
        for (int i = 1; i <= 15; ++i) {
            const int k = SparseMatrix::flatten(i,i);
            arr[k] += 1.0_rt;
        }
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    const Real& operator() (const int i, const int j) const noexcept {
        const int k = SparseMatrix::flatten(i,j);
        AMREX_ASSERT(k >= 0 && k < 107);
        return arr[k];
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real& operator() (const int i, const int j) noexcept {
        const int k = SparseMatrix::flatten(i,j);
        AMREX_ASSERT(k >= 0 && k < 107);
        return arr[k];
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void operator= (const SparseMatrix& other) noexcept {
        for (int i = 0; i < 107; ++i)
            arr[i] = other.arr[i];
    }

    Real arr[107];
};

template<class MatrixType>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void actual_solve(MatrixType const& A,
                  Array1D<Real, 1, 15>& b)
{
    const Real scratch_0 = 1.0/A(1,1);
    const Real scratch_1 = b(1)*scratch_0;
    const Real scratch_2 = A(1,3)*scratch_0;
    const Real scratch_3 = A(1,2)*scratch_0;
    const Real scratch_4 = 1.0/(-A(2,1)*scratch_3 + A(2,2));
    const Real scratch_5 = scratch_4*(-A(2,1)*scratch_2 + A(2,3));
    const Real scratch_6 = scratch_3*scratch_5;
    const Real scratch_7 = scratch_2 - scratch_6;
    const Real scratch_8 = -A(3,1)*scratch_3 + A(3,2);
    const Real scratch_9 = 1.0/(-A(3,1)*scratch_2 + A(3,3) - scratch_5*scratch_8);
    const Real scratch_10 = scratch_4*(-A(2,1)*scratch_1 + b(2));
    const Real scratch_11 = scratch_9*(-A(3,1)*scratch_1 + b(3) - scratch_10*scratch_8);
    const Real scratch_12 = A(1,4)*scratch_0;
    const Real scratch_13 = A(2,1)*scratch_4;
    const Real scratch_14 = scratch_12*scratch_13;
    const Real scratch_15 = scratch_9*(-A(3,1)*scratch_12 + A(3,4) + scratch_14*scratch_8);
    const Real scratch_16 = A(1,2)*scratch_13/std::pow(A(1,1), 2);
    const Real scratch_17 = A(1,4)*scratch_16;
    const Real scratch_18 = scratch_12 - scratch_15*scratch_7 + scratch_17;
    const Real scratch_19 = -A(4,1)*scratch_3 + A(4,2);
    const Real scratch_20 = -A(4,1)*scratch_2 + A(4,3) - scratch_19*scratch_5;
    const Real scratch_21 = 1.0/(-A(4,1)*scratch_12 + A(4,4) + scratch_14*scratch_19 - scratch_15*scratch_20);
    const Real scratch_22 = scratch_21*(-A(4,1)*scratch_1 + b(4) - scratch_10*scratch_19 - scratch_11*scratch_20);
    const Real scratch_23 = A(1,5)*scratch_0;
    const Real scratch_24 = scratch_13*scratch_23;
    const Real scratch_25 = scratch_9*(-A(3,1)*scratch_23 + scratch_24*scratch_8);
    const Real scratch_26 = scratch_21*(-A(4,1)*scratch_23 + A(4,5) + scratch_19*scratch_24 - scratch_20*scratch_25);
    const Real scratch_27 = A(1,5)*scratch_16;
    const Real scratch_28 = -scratch_18*scratch_26 + scratch_23 - scratch_25*scratch_7 + scratch_27;
    const Real scratch_29 = -A(5,1)*scratch_3 + A(5,2);
    const Real scratch_30 = -A(5,1)*scratch_2 + A(5,3) - scratch_29*scratch_5;
    const Real scratch_31 = -A(5,1)*scratch_12 + A(5,4) + scratch_14*scratch_29 - scratch_15*scratch_30;
    const Real scratch_32 = 1.0/(-A(5,1)*scratch_23 + A(5,5) + scratch_24*scratch_29 - scratch_25*scratch_30 - scratch_26*scratch_31);
    const Real scratch_33 = scratch_32*(-A(5,1)*scratch_1 + b(5) - scratch_10*scratch_29 - scratch_11*scratch_30 - scratch_22*scratch_31);
    const Real scratch_34 = A(1,6)*scratch_0;
    const Real scratch_35 = scratch_13*scratch_34;
    const Real scratch_36 = scratch_9*(-A(3,1)*scratch_34 + scratch_35*scratch_8);
    const Real scratch_37 = scratch_21*(-A(4,1)*scratch_34 + scratch_19*scratch_35 - scratch_20*scratch_36);
    const Real scratch_38 = scratch_32*(-A(5,1)*scratch_34 + A(5,6) + scratch_29*scratch_35 - scratch_30*scratch_36 - scratch_31*scratch_37);
    const Real scratch_39 = A(1,6)*scratch_16;
    const Real scratch_40 = -scratch_18*scratch_37 - scratch_28*scratch_38 + scratch_34 - scratch_36*scratch_7 + scratch_39;
    const Real scratch_41 = -A(6,1)*scratch_3 + A(6,2);
    const Real scratch_42 = -A(6,1)*scratch_2 + A(6,3) - scratch_41*scratch_5;
    const Real scratch_43 = -A(6,1)*scratch_12 + scratch_14*scratch_41 - scratch_15*scratch_42;
    const Real scratch_44 = -A(6,1)*scratch_23 + A(6,5) + scratch_24*scratch_41 - scratch_25*scratch_42 - scratch_26*scratch_43;
    const Real scratch_45 = 1.0/(-A(6,1)*scratch_34 + A(6,6) + scratch_35*scratch_41 - scratch_36*scratch_42 - scratch_37*scratch_43 - scratch_38*scratch_44);
    const Real scratch_46 = scratch_45*(-A(6,1)*scratch_1 + b(6) - scratch_10*scratch_41 - scratch_11*scratch_42 - scratch_22*scratch_43 - scratch_33*scratch_44);
    const Real scratch_47 = A(1,7)*scratch_0;
    const Real scratch_48 = scratch_13*scratch_47;
    const Real scratch_49 = scratch_9*(-A(3,1)*scratch_47 + scratch_48*scratch_8);
    const Real scratch_50 = scratch_21*(-A(4,1)*scratch_47 + scratch_19*scratch_48 - scratch_20*scratch_49);
    const Real scratch_51 = scratch_32*(-A(5,1)*scratch_47 + scratch_29*scratch_48 - scratch_30*scratch_49 - scratch_31*scratch_50);
    const Real scratch_52 = scratch_45*(-A(6,1)*scratch_47 + A(6,7) + scratch_41*scratch_48 - scratch_42*scratch_49 - scratch_43*scratch_50 - scratch_44*scratch_51);
    const Real scratch_53 = A(1,7)*scratch_16;
    const Real scratch_54 = -scratch_18*scratch_50 - scratch_28*scratch_51 - scratch_40*scratch_52 + scratch_47 - scratch_49*scratch_7 + scratch_53;
    const Real scratch_55 = -A(7,1)*scratch_2 + A(7,1)*scratch_6 + A(7,3);
    const Real scratch_56 = -A(7,1)*scratch_12 - A(7,1)*scratch_17 - scratch_15*scratch_55;
    const Real scratch_57 = -A(7,1)*scratch_23 - A(7,1)*scratch_27 - scratch_25*scratch_55 - scratch_26*scratch_56;
    const Real scratch_58 = -A(7,1)*scratch_34 - A(7,1)*scratch_39 + A(7,6) - scratch_36*scratch_55 - scratch_37*scratch_56 - scratch_38*scratch_57;
    const Real scratch_59 = 1.0/(-A(7,1)*scratch_47 - A(7,1)*scratch_53 + A(7,7) - scratch_49*scratch_55 - scratch_50*scratch_56 - scratch_51*scratch_57 - scratch_52*scratch_58);
    const Real scratch_60 = scratch_10*scratch_3;
    const Real scratch_61 = scratch_59*(-A(7,1)*scratch_1 + A(7,1)*scratch_60 + b(7) - scratch_11*scratch_55 - scratch_22*scratch_56 - scratch_33*scratch_57 - scratch_46*scratch_58);
    const Real scratch_62 = A(1,8)*scratch_0;
    const Real scratch_63 = scratch_13*scratch_62;
    const Real scratch_64 = scratch_9*(-A(3,1)*scratch_62 + scratch_63*scratch_8);
    const Real scratch_65 = scratch_21*(-A(4,1)*scratch_62 + scratch_19*scratch_63 - scratch_20*scratch_64);
    const Real scratch_66 = scratch_32*(-A(5,1)*scratch_62 + scratch_29*scratch_63 - scratch_30*scratch_64 - scratch_31*scratch_65);
    const Real scratch_67 = scratch_45*(-A(6,1)*scratch_62 + scratch_41*scratch_63 - scratch_42*scratch_64 - scratch_43*scratch_65 - scratch_44*scratch_66);
    const Real scratch_68 = A(1,8)*scratch_16;
    const Real scratch_69 = scratch_59*(-A(7,1)*scratch_62 - A(7,1)*scratch_68 + A(7,8) - scratch_55*scratch_64 - scratch_56*scratch_65 - scratch_57*scratch_66 - scratch_58*scratch_67);
    const Real scratch_70 = -scratch_18*scratch_65 - scratch_28*scratch_66 - scratch_40*scratch_67 - scratch_54*scratch_69 + scratch_62 - scratch_64*scratch_7 + scratch_68;
    const Real scratch_71 = -A(8,1)*scratch_2 + A(8,1)*scratch_6;
    const Real scratch_72 = -A(8,1)*scratch_12 - A(8,1)*scratch_17 - scratch_15*scratch_71;
    const Real scratch_73 = -A(8,1)*scratch_23 - A(8,1)*scratch_27 - scratch_25*scratch_71 - scratch_26*scratch_72;
    const Real scratch_74 = -A(8,1)*scratch_34 - A(8,1)*scratch_39 - scratch_36*scratch_71 - scratch_37*scratch_72 - scratch_38*scratch_73;
    const Real scratch_75 = -A(8,1)*scratch_47 - A(8,1)*scratch_53 + A(8,7) - scratch_49*scratch_71 - scratch_50*scratch_72 - scratch_51*scratch_73 - scratch_52*scratch_74;
    const Real scratch_76 = 1.0/(-A(8,1)*scratch_62 - A(8,1)*scratch_68 + A(8,8) - scratch_64*scratch_71 - scratch_65*scratch_72 - scratch_66*scratch_73 - scratch_67*scratch_74 - scratch_69*scratch_75);
    const Real scratch_77 = scratch_76*(-A(8,1)*scratch_1 + A(8,1)*scratch_60 + b(8) - scratch_11*scratch_71 - scratch_22*scratch_72 - scratch_33*scratch_73 - scratch_46*scratch_74 - scratch_61*scratch_75);
    const Real scratch_78 = A(1,9)*scratch_0;
    const Real scratch_79 = scratch_13*scratch_78;
    const Real scratch_80 = scratch_9*(-A(3,1)*scratch_78 + scratch_79*scratch_8);
    const Real scratch_81 = scratch_21*(-A(4,1)*scratch_78 + scratch_19*scratch_79 - scratch_20*scratch_80);
    const Real scratch_82 = scratch_32*(-A(5,1)*scratch_78 + scratch_29*scratch_79 - scratch_30*scratch_80 - scratch_31*scratch_81);
    const Real scratch_83 = scratch_45*(-A(6,1)*scratch_78 + scratch_41*scratch_79 - scratch_42*scratch_80 - scratch_43*scratch_81 - scratch_44*scratch_82);
    const Real scratch_84 = A(1,9)*scratch_16;
    const Real scratch_85 = scratch_59*(-A(7,1)*scratch_78 - A(7,1)*scratch_84 - scratch_55*scratch_80 - scratch_56*scratch_81 - scratch_57*scratch_82 - scratch_58*scratch_83);
    const Real scratch_86 = scratch_76*(-A(8,1)*scratch_78 - A(8,1)*scratch_84 + A(8,9) - scratch_71*scratch_80 - scratch_72*scratch_81 - scratch_73*scratch_82 - scratch_74*scratch_83 - scratch_75*scratch_85);
    const Real scratch_87 = -scratch_18*scratch_81 - scratch_28*scratch_82 - scratch_40*scratch_83 - scratch_54*scratch_85 - scratch_7*scratch_80 - scratch_70*scratch_86 + scratch_78 + scratch_84;
    const Real scratch_88 = -A(9,1)*scratch_2 + A(9,1)*scratch_6;
    const Real scratch_89 = -A(9,1)*scratch_12 - A(9,1)*scratch_17 - scratch_15*scratch_88;
    const Real scratch_90 = -A(9,1)*scratch_23 - A(9,1)*scratch_27 - scratch_25*scratch_88 - scratch_26*scratch_89;
    const Real scratch_91 = -A(9,1)*scratch_34 - A(9,1)*scratch_39 - scratch_36*scratch_88 - scratch_37*scratch_89 - scratch_38*scratch_90;
    const Real scratch_92 = -A(9,1)*scratch_47 - A(9,1)*scratch_53 - scratch_49*scratch_88 - scratch_50*scratch_89 - scratch_51*scratch_90 - scratch_52*scratch_91;
    const Real scratch_93 = -A(9,1)*scratch_62 - A(9,1)*scratch_68 + A(9,8) - scratch_64*scratch_88 - scratch_65*scratch_89 - scratch_66*scratch_90 - scratch_67*scratch_91 - scratch_69*scratch_92;
    const Real scratch_94 = 1.0/(-A(9,1)*scratch_78 - A(9,1)*scratch_84 + A(9,9) - scratch_80*scratch_88 - scratch_81*scratch_89 - scratch_82*scratch_90 - scratch_83*scratch_91 - scratch_85*scratch_92 - scratch_86*scratch_93);
    const Real scratch_95 = scratch_94*(-A(9,1)*scratch_1 + A(9,1)*scratch_60 + b(9) - scratch_11*scratch_88 - scratch_22*scratch_89 - scratch_33*scratch_90 - scratch_46*scratch_91 - scratch_61*scratch_92 - scratch_77*scratch_93);
    const Real scratch_96 = A(1,10)*scratch_0;
    const Real scratch_97 = scratch_13*scratch_96;
    const Real scratch_98 = scratch_9*(-A(3,1)*scratch_96 + scratch_8*scratch_97);
    const Real scratch_99 = scratch_21*(-A(4,1)*scratch_96 + scratch_19*scratch_97 - scratch_20*scratch_98);
    const Real scratch_100 = scratch_32*(-A(5,1)*scratch_96 + scratch_29*scratch_97 - scratch_30*scratch_98 - scratch_31*scratch_99);
    const Real scratch_101 = scratch_45*(-A(6,1)*scratch_96 - scratch_100*scratch_44 + scratch_41*scratch_97 - scratch_42*scratch_98 - scratch_43*scratch_99);
    const Real scratch_102 = A(1,10)*scratch_16;
    const Real scratch_103 = scratch_59*(-A(7,1)*scratch_102 - A(7,1)*scratch_96 - scratch_100*scratch_57 - scratch_101*scratch_58 - scratch_55*scratch_98 - scratch_56*scratch_99);
    const Real scratch_104 = scratch_76*(-A(8,1)*scratch_102 - A(8,1)*scratch_96 - scratch_100*scratch_73 - scratch_101*scratch_74 - scratch_103*scratch_75 - scratch_71*scratch_98 - scratch_72*scratch_99);
    const Real scratch_105 = scratch_94*(A(9,10) - A(9,1)*scratch_102 - A(9,1)*scratch_96 - scratch_100*scratch_90 - scratch_101*scratch_91 - scratch_103*scratch_92 - scratch_104*scratch_93 - scratch_88*scratch_98 - scratch_89*scratch_99);
    const Real scratch_106 = -scratch_100*scratch_28 - scratch_101*scratch_40 + scratch_102 - scratch_103*scratch_54 - scratch_104*scratch_70 - scratch_105*scratch_87 - scratch_18*scratch_99 - scratch_7*scratch_98 + scratch_96;
    const Real scratch_107 = -A(10,1)*scratch_2 + A(10,1)*scratch_6;
    const Real scratch_108 = -A(10,1)*scratch_12 - A(10,1)*scratch_17 - scratch_107*scratch_15;
    const Real scratch_109 = -A(10,1)*scratch_23 - A(10,1)*scratch_27 - scratch_107*scratch_25 - scratch_108*scratch_26;
    const Real scratch_110 = -A(10,1)*scratch_34 - A(10,1)*scratch_39 - scratch_107*scratch_36 - scratch_108*scratch_37 - scratch_109*scratch_38;
    const Real scratch_111 = -A(10,1)*scratch_47 - A(10,1)*scratch_53 - scratch_107*scratch_49 - scratch_108*scratch_50 - scratch_109*scratch_51 - scratch_110*scratch_52;
    const Real scratch_112 = -A(10,1)*scratch_62 - A(10,1)*scratch_68 - scratch_107*scratch_64 - scratch_108*scratch_65 - scratch_109*scratch_66 - scratch_110*scratch_67 - scratch_111*scratch_69;
    const Real scratch_113 = -A(10,1)*scratch_78 - A(10,1)*scratch_84 + A(10,9) - scratch_107*scratch_80 - scratch_108*scratch_81 - scratch_109*scratch_82 - scratch_110*scratch_83 - scratch_111*scratch_85 - scratch_112*scratch_86;
    const Real scratch_114 = 1.0/(A(10,10) - A(10,1)*scratch_102 - A(10,1)*scratch_96 - scratch_100*scratch_109 - scratch_101*scratch_110 - scratch_103*scratch_111 - scratch_104*scratch_112 - scratch_105*scratch_113 - scratch_107*scratch_98 - scratch_108*scratch_99);
    const Real scratch_115 = scratch_114*(-A(10,1)*scratch_1 + A(10,1)*scratch_60 + b(10) - scratch_107*scratch_11 - scratch_108*scratch_22 - scratch_109*scratch_33 - scratch_110*scratch_46 - scratch_111*scratch_61 - scratch_112*scratch_77 - scratch_113*scratch_95);
    const Real scratch_116 = A(1,11)*scratch_0;
    const Real scratch_117 = scratch_116*scratch_13;
    const Real scratch_118 = scratch_9*(-A(3,1)*scratch_116 + scratch_117*scratch_8);
    const Real scratch_119 = scratch_21*(-A(4,1)*scratch_116 + scratch_117*scratch_19 - scratch_118*scratch_20);
    const Real scratch_120 = scratch_32*(-A(5,1)*scratch_116 + scratch_117*scratch_29 - scratch_118*scratch_30 - scratch_119*scratch_31);
    const Real scratch_121 = scratch_45*(-A(6,1)*scratch_116 + scratch_117*scratch_41 - scratch_118*scratch_42 - scratch_119*scratch_43 - scratch_120*scratch_44);
    const Real scratch_122 = A(1,11)*scratch_16;
    const Real scratch_123 = scratch_59*(-A(7,1)*scratch_116 - A(7,1)*scratch_122 - scratch_118*scratch_55 - scratch_119*scratch_56 - scratch_120*scratch_57 - scratch_121*scratch_58);
    const Real scratch_124 = scratch_76*(-A(8,1)*scratch_116 - A(8,1)*scratch_122 - scratch_118*scratch_71 - scratch_119*scratch_72 - scratch_120*scratch_73 - scratch_121*scratch_74 - scratch_123*scratch_75);
    const Real scratch_125 = scratch_94*(-A(9,1)*scratch_116 - A(9,1)*scratch_122 - scratch_118*scratch_88 - scratch_119*scratch_89 - scratch_120*scratch_90 - scratch_121*scratch_91 - scratch_123*scratch_92 - scratch_124*scratch_93);
    const Real scratch_126 = scratch_114*(A(10,11) - A(10,1)*scratch_116 - A(10,1)*scratch_122 - scratch_107*scratch_118 - scratch_108*scratch_119 - scratch_109*scratch_120 - scratch_110*scratch_121 - scratch_111*scratch_123 - scratch_112*scratch_124 - scratch_113*scratch_125);
    const Real scratch_127 = -scratch_106*scratch_126 + scratch_116 - scratch_118*scratch_7 - scratch_119*scratch_18 - scratch_120*scratch_28 - scratch_121*scratch_40 + scratch_122 - scratch_123*scratch_54 - scratch_124*scratch_70 - scratch_125*scratch_87;
    const Real scratch_128 = -A(11,1)*scratch_2 + A(11,1)*scratch_6;
    const Real scratch_129 = -A(11,1)*scratch_12 - A(11,1)*scratch_17 - scratch_128*scratch_15;
    const Real scratch_130 = -A(11,1)*scratch_23 - A(11,1)*scratch_27 - scratch_128*scratch_25 - scratch_129*scratch_26;
    const Real scratch_131 = -A(11,1)*scratch_34 - A(11,1)*scratch_39 - scratch_128*scratch_36 - scratch_129*scratch_37 - scratch_130*scratch_38;
    const Real scratch_132 = -A(11,1)*scratch_47 - A(11,1)*scratch_53 - scratch_128*scratch_49 - scratch_129*scratch_50 - scratch_130*scratch_51 - scratch_131*scratch_52;
    const Real scratch_133 = -A(11,1)*scratch_62 - A(11,1)*scratch_68 - scratch_128*scratch_64 - scratch_129*scratch_65 - scratch_130*scratch_66 - scratch_131*scratch_67 - scratch_132*scratch_69;
    const Real scratch_134 = -A(11,1)*scratch_78 - A(11,1)*scratch_84 - scratch_128*scratch_80 - scratch_129*scratch_81 - scratch_130*scratch_82 - scratch_131*scratch_83 - scratch_132*scratch_85 - scratch_133*scratch_86;
    const Real scratch_135 = A(11,10) - A(11,1)*scratch_102 - A(11,1)*scratch_96 - scratch_100*scratch_130 - scratch_101*scratch_131 - scratch_103*scratch_132 - scratch_104*scratch_133 - scratch_105*scratch_134 - scratch_128*scratch_98 - scratch_129*scratch_99;
    const Real scratch_136 = 1.0/(A(11,11) - A(11,1)*scratch_116 - A(11,1)*scratch_122 - scratch_118*scratch_128 - scratch_119*scratch_129 - scratch_120*scratch_130 - scratch_121*scratch_131 - scratch_123*scratch_132 - scratch_124*scratch_133 - scratch_125*scratch_134 - scratch_126*scratch_135);
    const Real scratch_137 = scratch_136*(-A(11,1)*scratch_1 + A(11,1)*scratch_60 + b(11) - scratch_11*scratch_128 - scratch_115*scratch_135 - scratch_129*scratch_22 - scratch_130*scratch_33 - scratch_131*scratch_46 - scratch_132*scratch_61 - scratch_133*scratch_77 - scratch_134*scratch_95);
    const Real scratch_138 = A(1,12)*scratch_0;
    const Real scratch_139 = scratch_13*scratch_138;
    const Real scratch_140 = scratch_9*(-A(3,1)*scratch_138 + scratch_139*scratch_8);
    const Real scratch_141 = scratch_21*(-A(4,1)*scratch_138 + scratch_139*scratch_19 - scratch_140*scratch_20);
    const Real scratch_142 = scratch_32*(-A(5,1)*scratch_138 + scratch_139*scratch_29 - scratch_140*scratch_30 - scratch_141*scratch_31);
    const Real scratch_143 = scratch_45*(-A(6,1)*scratch_138 + scratch_139*scratch_41 - scratch_140*scratch_42 - scratch_141*scratch_43 - scratch_142*scratch_44);
    const Real scratch_144 = A(1,12)*scratch_16;
    const Real scratch_145 = scratch_59*(-A(7,1)*scratch_138 - A(7,1)*scratch_144 - scratch_140*scratch_55 - scratch_141*scratch_56 - scratch_142*scratch_57 - scratch_143*scratch_58);
    const Real scratch_146 = scratch_76*(-A(8,1)*scratch_138 - A(8,1)*scratch_144 - scratch_140*scratch_71 - scratch_141*scratch_72 - scratch_142*scratch_73 - scratch_143*scratch_74 - scratch_145*scratch_75);
    const Real scratch_147 = scratch_94*(-A(9,1)*scratch_138 - A(9,1)*scratch_144 - scratch_140*scratch_88 - scratch_141*scratch_89 - scratch_142*scratch_90 - scratch_143*scratch_91 - scratch_145*scratch_92 - scratch_146*scratch_93);
    const Real scratch_148 = scratch_114*(-A(10,1)*scratch_138 - A(10,1)*scratch_144 - scratch_107*scratch_140 - scratch_108*scratch_141 - scratch_109*scratch_142 - scratch_110*scratch_143 - scratch_111*scratch_145 - scratch_112*scratch_146 - scratch_113*scratch_147);
    const Real scratch_149 = scratch_136*(A(11,12) - A(11,1)*scratch_138 - A(11,1)*scratch_144 - scratch_128*scratch_140 - scratch_129*scratch_141 - scratch_130*scratch_142 - scratch_131*scratch_143 - scratch_132*scratch_145 - scratch_133*scratch_146 - scratch_134*scratch_147 - scratch_135*scratch_148);
    const Real scratch_150 = -scratch_106*scratch_148 - scratch_127*scratch_149 + scratch_138 - scratch_140*scratch_7 - scratch_141*scratch_18 - scratch_142*scratch_28 - scratch_143*scratch_40 + scratch_144 - scratch_145*scratch_54 - scratch_146*scratch_70 - scratch_147*scratch_87;
    const Real scratch_151 = -A(12,1)*scratch_2 + A(12,1)*scratch_6;
    const Real scratch_152 = -A(12,1)*scratch_12 - A(12,1)*scratch_17 - scratch_15*scratch_151;
    const Real scratch_153 = -A(12,1)*scratch_23 - A(12,1)*scratch_27 - scratch_151*scratch_25 - scratch_152*scratch_26;
    const Real scratch_154 = -A(12,1)*scratch_34 - A(12,1)*scratch_39 - scratch_151*scratch_36 - scratch_152*scratch_37 - scratch_153*scratch_38;
    const Real scratch_155 = -A(12,1)*scratch_47 - A(12,1)*scratch_53 - scratch_151*scratch_49 - scratch_152*scratch_50 - scratch_153*scratch_51 - scratch_154*scratch_52;
    const Real scratch_156 = -A(12,1)*scratch_62 - A(12,1)*scratch_68 - scratch_151*scratch_64 - scratch_152*scratch_65 - scratch_153*scratch_66 - scratch_154*scratch_67 - scratch_155*scratch_69;
    const Real scratch_157 = -A(12,1)*scratch_78 - A(12,1)*scratch_84 - scratch_151*scratch_80 - scratch_152*scratch_81 - scratch_153*scratch_82 - scratch_154*scratch_83 - scratch_155*scratch_85 - scratch_156*scratch_86;
    const Real scratch_158 = -A(12,1)*scratch_102 - A(12,1)*scratch_96 - scratch_100*scratch_153 - scratch_101*scratch_154 - scratch_103*scratch_155 - scratch_104*scratch_156 - scratch_105*scratch_157 - scratch_151*scratch_98 - scratch_152*scratch_99;
    const Real scratch_159 = A(12,11) - A(12,1)*scratch_116 - A(12,1)*scratch_122 - scratch_118*scratch_151 - scratch_119*scratch_152 - scratch_120*scratch_153 - scratch_121*scratch_154 - scratch_123*scratch_155 - scratch_124*scratch_156 - scratch_125*scratch_157 - scratch_126*scratch_158;
    const Real scratch_160 = 1.0/(A(12,12) - A(12,1)*scratch_138 - A(12,1)*scratch_144 - scratch_140*scratch_151 - scratch_141*scratch_152 - scratch_142*scratch_153 - scratch_143*scratch_154 - scratch_145*scratch_155 - scratch_146*scratch_156 - scratch_147*scratch_157 - scratch_148*scratch_158 - scratch_149*scratch_159);
    const Real scratch_161 = scratch_160*(-A(12,1)*scratch_1 + A(12,1)*scratch_60 + b(12) - scratch_11*scratch_151 - scratch_115*scratch_158 - scratch_137*scratch_159 - scratch_152*scratch_22 - scratch_153*scratch_33 - scratch_154*scratch_46 - scratch_155*scratch_61 - scratch_156*scratch_77 - scratch_157*scratch_95);
    const Real scratch_162 = A(1,13)*scratch_0;
    const Real scratch_163 = scratch_13*scratch_162;
    const Real scratch_164 = scratch_9*(-A(3,1)*scratch_162 + scratch_163*scratch_8);
    const Real scratch_165 = scratch_21*(-A(4,1)*scratch_162 + scratch_163*scratch_19 - scratch_164*scratch_20);
    const Real scratch_166 = scratch_32*(-A(5,1)*scratch_162 + scratch_163*scratch_29 - scratch_164*scratch_30 - scratch_165*scratch_31);
    const Real scratch_167 = scratch_45*(-A(6,1)*scratch_162 + scratch_163*scratch_41 - scratch_164*scratch_42 - scratch_165*scratch_43 - scratch_166*scratch_44);
    const Real scratch_168 = A(1,13)*scratch_16;
    const Real scratch_169 = scratch_59*(-A(7,1)*scratch_162 - A(7,1)*scratch_168 - scratch_164*scratch_55 - scratch_165*scratch_56 - scratch_166*scratch_57 - scratch_167*scratch_58);
    const Real scratch_170 = scratch_76*(-A(8,1)*scratch_162 - A(8,1)*scratch_168 - scratch_164*scratch_71 - scratch_165*scratch_72 - scratch_166*scratch_73 - scratch_167*scratch_74 - scratch_169*scratch_75);
    const Real scratch_171 = scratch_94*(-A(9,1)*scratch_162 - A(9,1)*scratch_168 - scratch_164*scratch_88 - scratch_165*scratch_89 - scratch_166*scratch_90 - scratch_167*scratch_91 - scratch_169*scratch_92 - scratch_170*scratch_93);
    const Real scratch_172 = scratch_114*(-A(10,1)*scratch_162 - A(10,1)*scratch_168 - scratch_107*scratch_164 - scratch_108*scratch_165 - scratch_109*scratch_166 - scratch_110*scratch_167 - scratch_111*scratch_169 - scratch_112*scratch_170 - scratch_113*scratch_171);
    const Real scratch_173 = scratch_136*(-A(11,1)*scratch_162 - A(11,1)*scratch_168 - scratch_128*scratch_164 - scratch_129*scratch_165 - scratch_130*scratch_166 - scratch_131*scratch_167 - scratch_132*scratch_169 - scratch_133*scratch_170 - scratch_134*scratch_171 - scratch_135*scratch_172);
    const Real scratch_174 = scratch_160*(A(12,13) - A(12,1)*scratch_162 - A(12,1)*scratch_168 - scratch_151*scratch_164 - scratch_152*scratch_165 - scratch_153*scratch_166 - scratch_154*scratch_167 - scratch_155*scratch_169 - scratch_156*scratch_170 - scratch_157*scratch_171 - scratch_158*scratch_172 - scratch_159*scratch_173);
    const Real scratch_175 = -scratch_106*scratch_172 - scratch_127*scratch_173 - scratch_150*scratch_174 + scratch_162 - scratch_164*scratch_7 - scratch_165*scratch_18 - scratch_166*scratch_28 - scratch_167*scratch_40 + scratch_168 - scratch_169*scratch_54 - scratch_170*scratch_70 - scratch_171*scratch_87;
    const Real scratch_176 = -A(13,1)*scratch_2 + A(13,1)*scratch_6;
    const Real scratch_177 = -A(13,1)*scratch_12 - A(13,1)*scratch_17 - scratch_15*scratch_176;
    const Real scratch_178 = -A(13,1)*scratch_23 - A(13,1)*scratch_27 - scratch_176*scratch_25 - scratch_177*scratch_26;
    const Real scratch_179 = -A(13,1)*scratch_34 - A(13,1)*scratch_39 - scratch_176*scratch_36 - scratch_177*scratch_37 - scratch_178*scratch_38;
    const Real scratch_180 = -A(13,1)*scratch_47 - A(13,1)*scratch_53 - scratch_176*scratch_49 - scratch_177*scratch_50 - scratch_178*scratch_51 - scratch_179*scratch_52;
    const Real scratch_181 = -A(13,1)*scratch_62 - A(13,1)*scratch_68 - scratch_176*scratch_64 - scratch_177*scratch_65 - scratch_178*scratch_66 - scratch_179*scratch_67 - scratch_180*scratch_69;
    const Real scratch_182 = -A(13,1)*scratch_78 - A(13,1)*scratch_84 - scratch_176*scratch_80 - scratch_177*scratch_81 - scratch_178*scratch_82 - scratch_179*scratch_83 - scratch_180*scratch_85 - scratch_181*scratch_86;
    const Real scratch_183 = -A(13,1)*scratch_102 - A(13,1)*scratch_96 - scratch_100*scratch_178 - scratch_101*scratch_179 - scratch_103*scratch_180 - scratch_104*scratch_181 - scratch_105*scratch_182 - scratch_176*scratch_98 - scratch_177*scratch_99;
    const Real scratch_184 = -A(13,1)*scratch_116 - A(13,1)*scratch_122 - scratch_118*scratch_176 - scratch_119*scratch_177 - scratch_120*scratch_178 - scratch_121*scratch_179 - scratch_123*scratch_180 - scratch_124*scratch_181 - scratch_125*scratch_182 - scratch_126*scratch_183;
    const Real scratch_185 = A(13,12) - A(13,1)*scratch_138 - A(13,1)*scratch_144 - scratch_140*scratch_176 - scratch_141*scratch_177 - scratch_142*scratch_178 - scratch_143*scratch_179 - scratch_145*scratch_180 - scratch_146*scratch_181 - scratch_147*scratch_182 - scratch_148*scratch_183 - scratch_149*scratch_184;
    const Real scratch_186 = 1.0/(A(13,13) - A(13,1)*scratch_162 - A(13,1)*scratch_168 - scratch_164*scratch_176 - scratch_165*scratch_177 - scratch_166*scratch_178 - scratch_167*scratch_179 - scratch_169*scratch_180 - scratch_170*scratch_181 - scratch_171*scratch_182 - scratch_172*scratch_183 - scratch_173*scratch_184 - scratch_174*scratch_185);
    const Real scratch_187 = scratch_186*(-A(13,1)*scratch_1 + A(13,1)*scratch_60 + b(13) - scratch_11*scratch_176 - scratch_115*scratch_183 - scratch_137*scratch_184 - scratch_161*scratch_185 - scratch_177*scratch_22 - scratch_178*scratch_33 - scratch_179*scratch_46 - scratch_180*scratch_61 - scratch_181*scratch_77 - scratch_182*scratch_95);
    const Real scratch_188 = A(1,14)*scratch_0;
    const Real scratch_189 = scratch_4*(A(2,14) - A(2,1)*scratch_188);
    const Real scratch_190 = scratch_9*(A(3,14) - A(3,1)*scratch_188 - scratch_189*scratch_8);
    const Real scratch_191 = scratch_21*(A(4,14) - A(4,1)*scratch_188 - scratch_189*scratch_19 - scratch_190*scratch_20);
    const Real scratch_192 = scratch_32*(A(5,14) - A(5,1)*scratch_188 - scratch_189*scratch_29 - scratch_190*scratch_30 - scratch_191*scratch_31);
    const Real scratch_193 = scratch_45*(A(6,14) - A(6,1)*scratch_188 - scratch_189*scratch_41 - scratch_190*scratch_42 - scratch_191*scratch_43 - scratch_192*scratch_44);
    const Real scratch_194 = scratch_189*scratch_3;
    const Real scratch_195 = scratch_59*(A(7,14) - A(7,1)*scratch_188 + A(7,1)*scratch_194 - scratch_190*scratch_55 - scratch_191*scratch_56 - scratch_192*scratch_57 - scratch_193*scratch_58);
    const Real scratch_196 = scratch_76*(A(8,14) - A(8,1)*scratch_188 + A(8,1)*scratch_194 - scratch_190*scratch_71 - scratch_191*scratch_72 - scratch_192*scratch_73 - scratch_193*scratch_74 - scratch_195*scratch_75);
    const Real scratch_197 = scratch_94*(A(9,14) - A(9,1)*scratch_188 + A(9,1)*scratch_194 - scratch_190*scratch_88 - scratch_191*scratch_89 - scratch_192*scratch_90 - scratch_193*scratch_91 - scratch_195*scratch_92 - scratch_196*scratch_93);
    const Real scratch_198 = scratch_114*(A(10,14) - A(10,1)*scratch_188 + A(10,1)*scratch_194 - scratch_107*scratch_190 - scratch_108*scratch_191 - scratch_109*scratch_192 - scratch_110*scratch_193 - scratch_111*scratch_195 - scratch_112*scratch_196 - scratch_113*scratch_197);
    const Real scratch_199 = scratch_136*(A(11,14) - A(11,1)*scratch_188 + A(11,1)*scratch_194 - scratch_128*scratch_190 - scratch_129*scratch_191 - scratch_130*scratch_192 - scratch_131*scratch_193 - scratch_132*scratch_195 - scratch_133*scratch_196 - scratch_134*scratch_197 - scratch_135*scratch_198);
    const Real scratch_200 = scratch_160*(A(12,14) - A(12,1)*scratch_188 + A(12,1)*scratch_194 - scratch_151*scratch_190 - scratch_152*scratch_191 - scratch_153*scratch_192 - scratch_154*scratch_193 - scratch_155*scratch_195 - scratch_156*scratch_196 - scratch_157*scratch_197 - scratch_158*scratch_198 - scratch_159*scratch_199);
    const Real scratch_201 = scratch_186*(A(13,14) - A(13,1)*scratch_188 + A(13,1)*scratch_194 - scratch_176*scratch_190 - scratch_177*scratch_191 - scratch_178*scratch_192 - scratch_179*scratch_193 - scratch_180*scratch_195 - scratch_181*scratch_196 - scratch_182*scratch_197 - scratch_183*scratch_198 - scratch_184*scratch_199 - scratch_185*scratch_200);
    const Real scratch_202 = -A(14,1)*scratch_3 + A(14,2);
    const Real scratch_203 = -A(14,1)*scratch_2 + A(14,3) - scratch_202*scratch_5;
    const Real scratch_204 = -A(14,1)*scratch_12 + A(14,4) + scratch_14*scratch_202 - scratch_15*scratch_203;
    const Real scratch_205 = -A(14,1)*scratch_23 + A(14,5) + scratch_202*scratch_24 - scratch_203*scratch_25 - scratch_204*scratch_26;
    const Real scratch_206 = -A(14,1)*scratch_34 + A(14,6) + scratch_202*scratch_35 - scratch_203*scratch_36 - scratch_204*scratch_37 - scratch_205*scratch_38;
    const Real scratch_207 = -A(14,1)*scratch_47 + A(14,7) + scratch_202*scratch_48 - scratch_203*scratch_49 - scratch_204*scratch_50 - scratch_205*scratch_51 - scratch_206*scratch_52;
    const Real scratch_208 = -A(14,1)*scratch_62 + A(14,8) + scratch_202*scratch_63 - scratch_203*scratch_64 - scratch_204*scratch_65 - scratch_205*scratch_66 - scratch_206*scratch_67 - scratch_207*scratch_69;
    const Real scratch_209 = -A(14,1)*scratch_78 + A(14,9) + scratch_202*scratch_79 - scratch_203*scratch_80 - scratch_204*scratch_81 - scratch_205*scratch_82 - scratch_206*scratch_83 - scratch_207*scratch_85 - scratch_208*scratch_86;
    const Real scratch_210 = A(14,10) - A(14,1)*scratch_96 - scratch_100*scratch_205 - scratch_101*scratch_206 - scratch_103*scratch_207 - scratch_104*scratch_208 - scratch_105*scratch_209 + scratch_202*scratch_97 - scratch_203*scratch_98 - scratch_204*scratch_99;
    const Real scratch_211 = A(14,11) - A(14,1)*scratch_116 + scratch_117*scratch_202 - scratch_118*scratch_203 - scratch_119*scratch_204 - scratch_120*scratch_205 - scratch_121*scratch_206 - scratch_123*scratch_207 - scratch_124*scratch_208 - scratch_125*scratch_209 - scratch_126*scratch_210;
    const Real scratch_212 = A(14,12) - A(14,1)*scratch_138 + scratch_139*scratch_202 - scratch_140*scratch_203 - scratch_141*scratch_204 - scratch_142*scratch_205 - scratch_143*scratch_206 - scratch_145*scratch_207 - scratch_146*scratch_208 - scratch_147*scratch_209 - scratch_148*scratch_210 - scratch_149*scratch_211;
    const Real scratch_213 = A(14,13) - A(14,1)*scratch_162 + scratch_163*scratch_202 - scratch_164*scratch_203 - scratch_165*scratch_204 - scratch_166*scratch_205 - scratch_167*scratch_206 - scratch_169*scratch_207 - scratch_170*scratch_208 - scratch_171*scratch_209 - scratch_172*scratch_210 - scratch_173*scratch_211 - scratch_174*scratch_212;
    const Real scratch_214 = (-A(14,1)*scratch_1 + b(14) - scratch_10*scratch_202 - scratch_11*scratch_203 - scratch_115*scratch_210 - scratch_137*scratch_211 - scratch_161*scratch_212 - scratch_187*scratch_213 - scratch_204*scratch_22 - scratch_205*scratch_33 - scratch_206*scratch_46 - scratch_207*scratch_61 - scratch_208*scratch_77 - scratch_209*scratch_95)/(A(14,14) - A(14,1)*scratch_188 - scratch_189*scratch_202 - scratch_190*scratch_203 - scratch_191*scratch_204 - scratch_192*scratch_205 - scratch_193*scratch_206 - scratch_195*scratch_207 - scratch_196*scratch_208 - scratch_197*scratch_209 - scratch_198*scratch_210 - scratch_199*scratch_211 - scratch_200*scratch_212 - scratch_201*scratch_213);
    const Real scratch_215 = -scratch_14 - scratch_15*scratch_5;
    const Real scratch_216 = -scratch_215*scratch_26 - scratch_24 - scratch_25*scratch_5;
    const Real scratch_217 = -scratch_215*scratch_37 - scratch_216*scratch_38 - scratch_35 - scratch_36*scratch_5;
    const Real scratch_218 = -scratch_215*scratch_50 - scratch_216*scratch_51 - scratch_217*scratch_52 - scratch_48 - scratch_49*scratch_5;
    const Real scratch_219 = -scratch_215*scratch_65 - scratch_216*scratch_66 - scratch_217*scratch_67 - scratch_218*scratch_69 - scratch_5*scratch_64 - scratch_63;
    const Real scratch_220 = -scratch_215*scratch_81 - scratch_216*scratch_82 - scratch_217*scratch_83 - scratch_218*scratch_85 - scratch_219*scratch_86 - scratch_5*scratch_80 - scratch_79;
    const Real scratch_221 = -scratch_100*scratch_216 - scratch_101*scratch_217 - scratch_103*scratch_218 - scratch_104*scratch_219 - scratch_105*scratch_220 - scratch_215*scratch_99 - scratch_5*scratch_98 - scratch_97;
    const Real scratch_222 = -scratch_117 - scratch_118*scratch_5 - scratch_119*scratch_215 - scratch_120*scratch_216 - scratch_121*scratch_217 - scratch_123*scratch_218 - scratch_124*scratch_219 - scratch_125*scratch_220 - scratch_126*scratch_221;
    const Real scratch_223 = -scratch_139 - scratch_140*scratch_5 - scratch_141*scratch_215 - scratch_142*scratch_216 - scratch_143*scratch_217 - scratch_145*scratch_218 - scratch_146*scratch_219 - scratch_147*scratch_220 - scratch_148*scratch_221 - scratch_149*scratch_222;
    const Real scratch_224 = -scratch_163 - scratch_164*scratch_5 - scratch_165*scratch_215 - scratch_166*scratch_216 - scratch_167*scratch_217 - scratch_169*scratch_218 - scratch_170*scratch_219 - scratch_171*scratch_220 - scratch_172*scratch_221 - scratch_173*scratch_222 - scratch_174*scratch_223;
    const Real scratch_225 = -scratch_15*scratch_26 + scratch_25;
    const Real scratch_226 = -scratch_15*scratch_37 - scratch_225*scratch_38 + scratch_36;
    const Real scratch_227 = -scratch_15*scratch_50 - scratch_225*scratch_51 - scratch_226*scratch_52 + scratch_49;
    const Real scratch_228 = -scratch_15*scratch_65 - scratch_225*scratch_66 - scratch_226*scratch_67 - scratch_227*scratch_69 + scratch_64;
    const Real scratch_229 = -scratch_15*scratch_81 - scratch_225*scratch_82 - scratch_226*scratch_83 - scratch_227*scratch_85 - scratch_228*scratch_86 + scratch_80;
    const Real scratch_230 = -scratch_100*scratch_225 - scratch_101*scratch_226 - scratch_103*scratch_227 - scratch_104*scratch_228 - scratch_105*scratch_229 - scratch_15*scratch_99 + scratch_98;
    const Real scratch_231 = scratch_118 - scratch_119*scratch_15 - scratch_120*scratch_225 - scratch_121*scratch_226 - scratch_123*scratch_227 - scratch_124*scratch_228 - scratch_125*scratch_229 - scratch_126*scratch_230;
    const Real scratch_232 = scratch_140 - scratch_141*scratch_15 - scratch_142*scratch_225 - scratch_143*scratch_226 - scratch_145*scratch_227 - scratch_146*scratch_228 - scratch_147*scratch_229 - scratch_148*scratch_230 - scratch_149*scratch_231;
    const Real scratch_233 = -scratch_15*scratch_165 + scratch_164 - scratch_166*scratch_225 - scratch_167*scratch_226 - scratch_169*scratch_227 - scratch_170*scratch_228 - scratch_171*scratch_229 - scratch_172*scratch_230 - scratch_173*scratch_231 - scratch_174*scratch_232;
    const Real scratch_234 = -scratch_26*scratch_38 + scratch_37;
    const Real scratch_235 = -scratch_234*scratch_52 - scratch_26*scratch_51 + scratch_50;
    const Real scratch_236 = -scratch_234*scratch_67 - scratch_235*scratch_69 - scratch_26*scratch_66 + scratch_65;
    const Real scratch_237 = -scratch_234*scratch_83 - scratch_235*scratch_85 - scratch_236*scratch_86 - scratch_26*scratch_82 + scratch_81;
    const Real scratch_238 = -scratch_100*scratch_26 - scratch_101*scratch_234 - scratch_103*scratch_235 - scratch_104*scratch_236 - scratch_105*scratch_237 + scratch_99;
    const Real scratch_239 = scratch_119 - scratch_120*scratch_26 - scratch_121*scratch_234 - scratch_123*scratch_235 - scratch_124*scratch_236 - scratch_125*scratch_237 - scratch_126*scratch_238;
    const Real scratch_240 = scratch_141 - scratch_142*scratch_26 - scratch_143*scratch_234 - scratch_145*scratch_235 - scratch_146*scratch_236 - scratch_147*scratch_237 - scratch_148*scratch_238 - scratch_149*scratch_239;
    const Real scratch_241 = scratch_165 - scratch_166*scratch_26 - scratch_167*scratch_234 - scratch_169*scratch_235 - scratch_170*scratch_236 - scratch_171*scratch_237 - scratch_172*scratch_238 - scratch_173*scratch_239 - scratch_174*scratch_240;
    const Real scratch_242 = -scratch_38*scratch_52 + scratch_51;
    const Real scratch_243 = -scratch_242*scratch_69 - scratch_38*scratch_67 + scratch_66;
    const Real scratch_244 = -scratch_242*scratch_85 - scratch_243*scratch_86 - scratch_38*scratch_83 + scratch_82;
    const Real scratch_245 = scratch_100 - scratch_101*scratch_38 - scratch_103*scratch_242 - scratch_104*scratch_243 - scratch_105*scratch_244;
    const Real scratch_246 = scratch_120 - scratch_121*scratch_38 - scratch_123*scratch_242 - scratch_124*scratch_243 - scratch_125*scratch_244 - scratch_126*scratch_245;
    const Real scratch_247 = scratch_142 - scratch_143*scratch_38 - scratch_145*scratch_242 - scratch_146*scratch_243 - scratch_147*scratch_244 - scratch_148*scratch_245 - scratch_149*scratch_246;
    const Real scratch_248 = scratch_166 - scratch_167*scratch_38 - scratch_169*scratch_242 - scratch_170*scratch_243 - scratch_171*scratch_244 - scratch_172*scratch_245 - scratch_173*scratch_246 - scratch_174*scratch_247;
    const Real scratch_249 = -scratch_52*scratch_69 + scratch_67;
    const Real scratch_250 = -scratch_249*scratch_86 - scratch_52*scratch_85 + scratch_83;
    const Real scratch_251 = scratch_101 - scratch_103*scratch_52 - scratch_104*scratch_249 - scratch_105*scratch_250;
    const Real scratch_252 = scratch_121 - scratch_123*scratch_52 - scratch_124*scratch_249 - scratch_125*scratch_250 - scratch_126*scratch_251;
    const Real scratch_253 = scratch_143 - scratch_145*scratch_52 - scratch_146*scratch_249 - scratch_147*scratch_250 - scratch_148*scratch_251 - scratch_149*scratch_252;
    const Real scratch_254 = scratch_167 - scratch_169*scratch_52 - scratch_170*scratch_249 - scratch_171*scratch_250 - scratch_172*scratch_251 - scratch_173*scratch_252 - scratch_174*scratch_253;
    const Real scratch_255 = -scratch_69*scratch_86 + scratch_85;
    const Real scratch_256 = scratch_103 - scratch_104*scratch_69 - scratch_105*scratch_255;
    const Real scratch_257 = scratch_123 - scratch_124*scratch_69 - scratch_125*scratch_255 - scratch_126*scratch_256;
    const Real scratch_258 = scratch_145 - scratch_146*scratch_69 - scratch_147*scratch_255 - scratch_148*scratch_256 - scratch_149*scratch_257;
    const Real scratch_259 = scratch_169 - scratch_170*scratch_69 - scratch_171*scratch_255 - scratch_172*scratch_256 - scratch_173*scratch_257 - scratch_174*scratch_258;
    const Real scratch_260 = scratch_104 - scratch_105*scratch_86;
    const Real scratch_261 = scratch_124 - scratch_125*scratch_86 - scratch_126*scratch_260;
    const Real scratch_262 = scratch_146 - scratch_147*scratch_86 - scratch_148*scratch_260 - scratch_149*scratch_261;
    const Real scratch_263 = scratch_170 - scratch_171*scratch_86 - scratch_172*scratch_260 - scratch_173*scratch_261 - scratch_174*scratch_262;
    const Real scratch_264 = -scratch_105*scratch_126 + scratch_125;
    const Real scratch_265 = -scratch_105*scratch_148 + scratch_147 - scratch_149*scratch_264;
    const Real scratch_266 = -scratch_105*scratch_172 + scratch_171 - scratch_173*scratch_264 - scratch_174*scratch_265;
    const Real scratch_267 = -scratch_126*scratch_149 + scratch_148;
    const Real scratch_268 = -scratch_126*scratch_173 + scratch_172 - scratch_174*scratch_267;
    const Real scratch_269 = -scratch_149*scratch_174 + scratch_173;
    const Real scratch_270 = -A(15,1)*scratch_3 + A(15,2);
    const Real scratch_271 = -A(15,1)*scratch_2 + A(15,3) - scratch_270*scratch_5;
    const Real scratch_272 = -A(15,1)*scratch_12 + A(15,4) + scratch_14*scratch_270 - scratch_15*scratch_271;
    const Real scratch_273 = -A(15,1)*scratch_23 + A(15,5) + scratch_24*scratch_270 - scratch_25*scratch_271 - scratch_26*scratch_272;
    const Real scratch_274 = -A(15,1)*scratch_34 + A(15,6) + scratch_270*scratch_35 - scratch_271*scratch_36 - scratch_272*scratch_37 - scratch_273*scratch_38;
    const Real scratch_275 = -A(15,1)*scratch_47 + A(15,7) + scratch_270*scratch_48 - scratch_271*scratch_49 - scratch_272*scratch_50 - scratch_273*scratch_51 - scratch_274*scratch_52;
    const Real scratch_276 = -A(15,1)*scratch_62 + A(15,8) + scratch_270*scratch_63 - scratch_271*scratch_64 - scratch_272*scratch_65 - scratch_273*scratch_66 - scratch_274*scratch_67 - scratch_275*scratch_69;
    const Real scratch_277 = -A(15,1)*scratch_78 + A(15,9) + scratch_270*scratch_79 - scratch_271*scratch_80 - scratch_272*scratch_81 - scratch_273*scratch_82 - scratch_274*scratch_83 - scratch_275*scratch_85 - scratch_276*scratch_86;
    const Real scratch_278 = A(15,10) - A(15,1)*scratch_96 - scratch_100*scratch_273 - scratch_101*scratch_274 - scratch_103*scratch_275 - scratch_104*scratch_276 - scratch_105*scratch_277 + scratch_270*scratch_97 - scratch_271*scratch_98 - scratch_272*scratch_99;
    const Real scratch_279 = A(15,11) - A(15,1)*scratch_116 + scratch_117*scratch_270 - scratch_118*scratch_271 - scratch_119*scratch_272 - scratch_120*scratch_273 - scratch_121*scratch_274 - scratch_123*scratch_275 - scratch_124*scratch_276 - scratch_125*scratch_277 - scratch_126*scratch_278;
    const Real scratch_280 = A(15,12) - A(15,1)*scratch_138 + scratch_139*scratch_270 - scratch_140*scratch_271 - scratch_141*scratch_272 - scratch_142*scratch_273 - scratch_143*scratch_274 - scratch_145*scratch_275 - scratch_146*scratch_276 - scratch_147*scratch_277 - scratch_148*scratch_278 - scratch_149*scratch_279;
    const Real scratch_281 = A(15,13) - A(15,1)*scratch_162 + scratch_163*scratch_270 - scratch_164*scratch_271 - scratch_165*scratch_272 - scratch_166*scratch_273 - scratch_167*scratch_274 - scratch_169*scratch_275 - scratch_170*scratch_276 - scratch_171*scratch_277 - scratch_172*scratch_278 - scratch_173*scratch_279 - scratch_174*scratch_280;

    Array1D<Real, 1, 15> x;

    x(1) = scratch_1 - scratch_106*scratch_115 - scratch_11*scratch_7 - scratch_127*scratch_137 - scratch_150*scratch_161 - scratch_175*scratch_187 - scratch_18*scratch_22 - scratch_214*(-scratch_106*scratch_198 - scratch_127*scratch_199 - scratch_150*scratch_200 - scratch_175*scratch_201 - scratch_18*scratch_191 + scratch_188 - scratch_190*scratch_7 - scratch_192*scratch_28 - scratch_193*scratch_40 - scratch_194 - scratch_195*scratch_54 - scratch_196*scratch_70 - scratch_197*scratch_87) - scratch_28*scratch_33 - scratch_40*scratch_46 - scratch_54*scratch_61 - scratch_60 - scratch_70*scratch_77 - scratch_87*scratch_95;
    x(2) = scratch_10 - scratch_11*scratch_5 - scratch_115*scratch_221 - scratch_137*scratch_222 - scratch_161*scratch_223 - scratch_187*scratch_224 - scratch_214*(scratch_189 - scratch_190*scratch_5 - scratch_191*scratch_215 - scratch_192*scratch_216 - scratch_193*scratch_217 - scratch_195*scratch_218 - scratch_196*scratch_219 - scratch_197*scratch_220 - scratch_198*scratch_221 - scratch_199*scratch_222 - scratch_200*scratch_223 - scratch_201*scratch_224) - scratch_215*scratch_22 - scratch_216*scratch_33 - scratch_217*scratch_46 - scratch_218*scratch_61 - scratch_219*scratch_77 - scratch_220*scratch_95;
    x(3) = scratch_11 - scratch_115*scratch_230 - scratch_137*scratch_231 - scratch_15*scratch_22 - scratch_161*scratch_232 - scratch_187*scratch_233 - scratch_214*(-scratch_15*scratch_191 + scratch_190 - scratch_192*scratch_225 - scratch_193*scratch_226 - scratch_195*scratch_227 - scratch_196*scratch_228 - scratch_197*scratch_229 - scratch_198*scratch_230 - scratch_199*scratch_231 - scratch_200*scratch_232 - scratch_201*scratch_233) - scratch_225*scratch_33 - scratch_226*scratch_46 - scratch_227*scratch_61 - scratch_228*scratch_77 - scratch_229*scratch_95;
    x(4) = -scratch_115*scratch_238 - scratch_137*scratch_239 - scratch_161*scratch_240 - scratch_187*scratch_241 - scratch_214*(scratch_191 - scratch_192*scratch_26 - scratch_193*scratch_234 - scratch_195*scratch_235 - scratch_196*scratch_236 - scratch_197*scratch_237 - scratch_198*scratch_238 - scratch_199*scratch_239 - scratch_200*scratch_240 - scratch_201*scratch_241) + scratch_22 - scratch_234*scratch_46 - scratch_235*scratch_61 - scratch_236*scratch_77 - scratch_237*scratch_95 - scratch_26*scratch_33;
    x(5) = -scratch_115*scratch_245 - scratch_137*scratch_246 - scratch_161*scratch_247 - scratch_187*scratch_248 - scratch_214*(scratch_192 - scratch_193*scratch_38 - scratch_195*scratch_242 - scratch_196*scratch_243 - scratch_197*scratch_244 - scratch_198*scratch_245 - scratch_199*scratch_246 - scratch_200*scratch_247 - scratch_201*scratch_248) - scratch_242*scratch_61 - scratch_243*scratch_77 - scratch_244*scratch_95 + scratch_33 - scratch_38*scratch_46;
    x(6) = -scratch_115*scratch_251 - scratch_137*scratch_252 - scratch_161*scratch_253 - scratch_187*scratch_254 - scratch_214*(scratch_193 - scratch_195*scratch_52 - scratch_196*scratch_249 - scratch_197*scratch_250 - scratch_198*scratch_251 - scratch_199*scratch_252 - scratch_200*scratch_253 - scratch_201*scratch_254) - scratch_249*scratch_77 - scratch_250*scratch_95 + scratch_46 - scratch_52*scratch_61;
    x(7) = -scratch_115*scratch_256 - scratch_137*scratch_257 - scratch_161*scratch_258 - scratch_187*scratch_259 - scratch_214*(scratch_195 - scratch_196*scratch_69 - scratch_197*scratch_255 - scratch_198*scratch_256 - scratch_199*scratch_257 - scratch_200*scratch_258 - scratch_201*scratch_259) - scratch_255*scratch_95 + scratch_61 - scratch_69*scratch_77;
    x(8) = -scratch_115*scratch_260 - scratch_137*scratch_261 - scratch_161*scratch_262 - scratch_187*scratch_263 - scratch_214*(scratch_196 - scratch_197*scratch_86 - scratch_198*scratch_260 - scratch_199*scratch_261 - scratch_200*scratch_262 - scratch_201*scratch_263) + scratch_77 - scratch_86*scratch_95;
    x(9) = -scratch_105*scratch_115 - scratch_137*scratch_264 - scratch_161*scratch_265 - scratch_187*scratch_266 - scratch_214*(-scratch_105*scratch_198 + scratch_197 - scratch_199*scratch_264 - scratch_200*scratch_265 - scratch_201*scratch_266) + scratch_95;
    x(10) = scratch_115 - scratch_126*scratch_137 - scratch_161*scratch_267 - scratch_187*scratch_268 - scratch_214*(-scratch_126*scratch_199 + scratch_198 - scratch_200*scratch_267 - scratch_201*scratch_268);
    x(11) = scratch_137 - scratch_149*scratch_161 - scratch_187*scratch_269 - scratch_214*(-scratch_149*scratch_200 + scratch_199 - scratch_201*scratch_269);
    x(12) = scratch_161 - scratch_174*scratch_187 - scratch_214*(-scratch_174*scratch_201 + scratch_200);
    x(13) = scratch_187 - scratch_201*scratch_214;
    x(14) = scratch_214;
    x(15) = (-A(15,1)*scratch_1 + b(15) - scratch_10*scratch_270 - scratch_11*scratch_271 - scratch_115*scratch_278 - scratch_137*scratch_279 - scratch_161*scratch_280 - scratch_187*scratch_281 - scratch_214*(A(15,14) - A(15,1)*scratch_188 - scratch_189*scratch_270 - scratch_190*scratch_271 - scratch_191*scratch_272 - scratch_192*scratch_273 - scratch_193*scratch_274 - scratch_195*scratch_275 - scratch_196*scratch_276 - scratch_197*scratch_277 - scratch_198*scratch_278 - scratch_199*scratch_279 - scratch_200*scratch_280 - scratch_201*scratch_281) - scratch_22*scratch_272 - scratch_273*scratch_33 - scratch_274*scratch_46 - scratch_275*scratch_61 - scratch_276*scratch_77 - scratch_277*scratch_95)/A(15,15);

    for (int i = 1; i <= 15; ++i) {
        b(i) = x(i);
    }

}
#endif
