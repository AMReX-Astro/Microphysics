#ifndef _sparse_jacobian_type_H_
#define _sparse_jacobian_type_H_

#include <AMReX_Array.H>
#include <network_properties.H>
#include <burn_type.H>

using namespace amrex;
using namespace Species;

struct SparseJacobian : Array1D<Real, 1, 107>
{
    static constexpr int flatten(const int i, const int j) {
        if (i==He4) {
            if (j==He4)
                return 1;
            else if (j==C12)
                return 1+1;
            else if (j==O16)
                return 1+2;
            else if (j==Ne20)
                return 1+3;
            else if (j==Mg24)
                return 1+4;
            else if (j==Si28)
                return 1+5;
            else if (j==S32)
                return 1+6;
            else if (j==Ar36)
                return 1+7;
            else if (j==Ca40)
                return 1+8;
            else if (j==Ti44)
                return 1+9;
            else if (j==Cr48)
                return 1+10;
            else if (j==Fe52)
                return 1+11;
            else if (j==Ni56)
                return 1+12;
            else if (j==net_itemp)
                return 1+13;
        } else if (i==C12) {
            if (j==He4)
                return 1+14;
            else if (j==C12)
                return 1+15;
            else if (j==O16)
                return 1+16;
            else if (j==net_itemp)
                return 1+17;
        } else if (i==O16) {
            if (j==He4)
                return 1+18;
            else if (j==C12)
                return 1+19;
            else if (j==O16)
                return 1+20;
            else if (j==Ne20)
                return 1+21;
            else if (j==net_itemp)
                return 1+22;
        } else if (i==Ne20) {
            if (j==He4)
                return 1+23;
            else if (j==C12)
                return 1+24;
            else if (j==O16)
                return 1+25;
            else if (j==Ne20)
                return 1+26;
            else if (j==Mg24)
                return 1+27;
            else if (j==net_itemp)
                return 1+28;
        } else if (i==Mg24) {
            if (j==He4)
                return 1+29;
            else if (j==C12)
                return 1+30;
            else if (j==O16)
                return 1+31;
            else if (j==Ne20)
                return 1+32;
            else if (j==Mg24)
                return 1+33;
            else if (j==Si28)
                return 1+34;
            else if (j==net_itemp)
                return 1+35;
        } else if (i==Si28) {
            if (j==He4)
                return 1+36;
            else if (j==C12)
                return 1+37;
            else if (j==O16)
                return 1+38;
            else if (j==Mg24)
                return 1+39;
            else if (j==Si28)
                return 1+40;
            else if (j==S32)
                return 1+41;
            else if (j==net_itemp)
                return 1+42;
        } else if (i==S32) {
            if (j==He4)
                return 1+43;
            else if (j==O16)
                    return 1+44;
                else if (j==Si28)
                    return 1+45;
                else if (j==S32)
                    return 1+46;
                else if (j==Ar36)
                    return 1+47;
                else if (j==net_itemp)
                    return 1+48;
            } else if (i==Ar36) {
                if (j==He4)
                    return 1+49;
                else if (j==S32)
                    return 1+50;
                else if (j==Ar36)
                    return 1+51;
                else if (j==Ca40)
                    return 1+52;
                else if (j==net_itemp)
                    return 1+53;
            } else if (i==Ca40) {
                if (j==He4)
                    return 1+54;
                else if (j==Ar36)
                    return 1+55;
                else if (j==Ca40)
                    return 1+56;
                else if (j==Ti44)
                    return 1+57;
                else if (j==net_itemp)
                    return 1+58;
            } else if (i==Ti44) {
                if (j==He4)
                    return 1+59;
                else if (j==Ca40)
                    return 1+60;
                else if (j==Ti44)
                    return 1+61;
                else if (j==Cr48)
                    return 1+62;
                else if (j==net_itemp)
                    return 1+63;
            } else if (i==Cr48) {
                if (j==He4)
                    return 1+64;
                else if (j==Ti44)
                    return 1+65;
                else if (j==Cr48)
                    return 1+66;
                else if (j==Fe52)
                    return 1+67;
                else if (j==net_itemp)
                    return 1+68;
            } else if (i==Fe52) {
                if (j==He4)
                    return 1+69;
                else if (j==Cr48)
                    return 1+70;
                else if (j==Fe52)
                    return 1+71;
                else if (j==Ni56)
                    return 1+72;
                else if (j==net_itemp)
                    return 1+73;
            } else if (i==Ni56) {
                if (j==He4)
                    return 1+74;
                else if (j==Fe52)
                    return 1+75;
                else if (j==Ni56)
                    return 1+76;
                else if (j==net_itemp)
                    return 1+77;
            } else if (i==net_itemp) {
                if (j==He4)
                    return 1+78;
                else if (j==C12)
                    return 1+79;
                else if (j==O16)
                    return 1+80;
                else if (j==Ne20)
                    return 1+81;
                else if (j==Mg24)
                    return 1+82;
                else if (j==Si28)
                    return 1+83;
                else if (j==S32)
                    return 1+84;
                else if (j==Ar36)
                    return 1+85;
                else if (j==Ca40)
                    return 1+86;
                else if (j==Ti44)
                    return 1+87;
                else if (j==Cr48)
                    return 1+88;
                else if (j==Fe52)
                    return 1+89;
                else if (j==Ni56)
                    return 1+90;
                else if (j==net_itemp)
                    return 1+91;
            } else if (i==net_ienuc) {
                if (j==He4)
                    return 1+92;
                else if (j==C12)
                    return 1+93;
                else if (j==O16)
                    return 1+94;
                else if (j==Ne20)
                    return 1+95;
                else if (j==Mg24)
                    return 1+96;
                else if (j==Si28)
                    return 1+97;
                else if (j==S32)
                    return 1+98;
                else if (j==Ar36)
                    return 1+99;
                else if (j==Ca40)
                    return 1+100;
                else if (j==Ti44)
                    return 1+101;
                else if (j==Cr48)
                    return 1+102;
                else if (j==Fe52)
                    return 1+103;
                else if (j==Ni56)
                    return 1+104;
                else if (j==net_itemp)
                    return 1+105;
                else if (j==net_ienuc)
                    return 1+106;
            }

            return 0;
        }


        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void zero () const noexcept {
            for (int i = 0; i < 107; ++i)
                arr[i] = 0.0_rt;
        }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void mul (const Real x) const noexcept {
            for (int i = 0; i < 107; ++i)
                arr[i] *= x;
        }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void add_identity () const noexcept {
            for (int i = 1; i <= neqs; ++i) {
                const int k = SparseJacobian::flatten(i,i);
                arr[k-1] += 1.0_rt;
            }
        }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        const Real& operator() (const int i, const int j) const noexcept {
            const int k = SparseJacobian::flatten(i,j);
            return arr[k-1];
        }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        Real& operator() (const int i, const int j) noexcept {
            const int k = SparseJacobian::flatten(i,j);
            return arr[k-1];
        }
    };

#endif