// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
#pragma once

// Binary64 pow, log10, and tanh for the Loom translation. The upstream AMDGPU
// lowering does not support scalar.powf, log10f, or tanhf on f64. These are
// statement-for-statement transcriptions of AMD's OCML (ROCm-Device-Libs at
// llvm-project amd/device-libs d366fa84f3fd):
//   pow    powD.cl -> powD_base.h (COMPILING_POW), with eplnD.cl and
//          expepD.cl (EXTRA_ACCURACY is not defined, so expep is exp(x.hi)
//          corrected by x.lo), and the ep.h double-double helpers expanded
//          in place.
//   log10  log10D.cl -> logD_base.h (COMPILING_LOG10), with ep.h.
//   tanh   tanhD.cl with epexpepD.cl and ep.h.
// They are what hipcc calls for std::pow(double, double), std::log10(double),
// and std::tanh(double) (__ocml_pow_f64, __ocml_log10_f64, __ocml_tanh_f64).
// The HIP build links OCML with the default controls: __oclc_finite_only_opt
// = 0, __oclc_unsafe_math_opt = 0, f64 denormals on. So the non-finite
// branches below are the ones OCML takes, and no DAZ path applies to f64.
//
// The same conventions as f64_math.H apply, which this file reuses (fma and
// exp) and therefore extends in namespace chemistry::f64_math:
//   rint(v)      (v + 0x1.8p52) - 0x1.8p52, exact for |v| < 2^51.
//   ldexp(v, n)  Multiplication by a power of two built from exponent bits.
//   frexp(v)     The exponent and mantissa fields, after scaling subnormals by
//                2^54.
//   MATH_FAST_RCP
//                IEEE f64 division. OCML refines the hardware v_rcp_f64 with
//                two fma Newton steps, which gives the correctly rounded
//                reciprocal except possibly in rare cases; GPU verification
//                decides whether any input differs.
//   trunc, copysign
//                Loom scalar operations (v_trunc_f64; bit operations).
// One operation per statement keeps OCML's evaluation order and prevents
// contraction. Every conditional selects between values that are already
// computed. Where OCML computes a value that a later select discards and the
// transcription's integer exponent arithmetic would be out of range for it,
// the input is clamped first; the comment at each clamp says why no returned
// value changes. Remove this file and bind math::pow/log10/tanh back to the
// scalar operations once Loom lowers them for f64.
//
// The algorithms and coefficients are adapted from ROCm-Device-Libs, which
// carries the following license.
//
// University of Illinois/NCSA Open Source License
//
// Copyright (c) 2014-2016, Advanced Micro Devices, Inc.
// All rights reserved.
//
// Developed by:
//
//     AMD Research and AMD HSA Software Development
//
//     Advanced Micro Devices, Inc.
//
//     www.amd.com
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// with the Software without restriction, including without limitation the
// rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
// sell copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
//     * Redistributions of source code must retain the above copyright notice,
//       this list of conditions and the following disclaimers.
//
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimers in the
//       documentation and/or other materials provided with the distribution.
//
//     * Neither the names of the LLVM Team, University of Illinois at
//       Urbana-Champaign, nor the names of its contributors may be used to
//       endorse or promote products derived from this Software without
//       specific prior written permission.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS WITH
// THE SOFTWARE.

#include <loomcxx/scalar.h>

#include "f64_math.H"

namespace chemistry {
namespace f64_math {

// double2 as returned by the ep.h helpers.
struct ep_double {
    double hi;
    double lo;
};

// BUILTIN_TRUNC_F64 (v_trunc_f64).
DEVICE double trunc(double x) {
    return loom::scalar::truncf(x);
}

// BUILTIN_COPYSIGN_F64: the magnitude of x with the sign of y.
DEVICE double copysign(double x, double y) {
    const u64 x_bits = __builtin_bit_cast(u64, x);
    const u64 y_bits = __builtin_bit_cast(u64, y);
    const u64 magnitude = x_bits & 0x7fffffffffffffffULL;
    const u64 sign = y_bits & 0x8000000000000000ULL;
    const u64 result_bits = magnitude | sign;
    return __builtin_bit_cast(double, result_bits);
}

// eplnD.cl: ln(a) as a double-double, for finite a > 0. Other inputs produce
// unused values that the callers' selects replace.
DEVICE ep_double epln(double a) {
    const u64 bits = __builtin_bit_cast(u64, a);
    const u64 magnitude = bits & 0x7fffffffffffffffULL;

    // m = frexp_mant(a) in [0.5, 1) and its exponent.
    const bool subnormal = magnitude < 0x0010000000000000ULL;
    const double scaled_subnormal = a * 0x1p54;
    const double normalized_input = subnormal ? scaled_subnormal : a;
    const u64 normalized_bits = __builtin_bit_cast(u64, normalized_input);
    const u64 shifted_exponent = normalized_bits >> 52;
    const u64 exponent_bits = shifted_exponent & 2047;
    const int biased_exponent = (int)exponent_bits;
    const int exponent_offset = subnormal ? 1076 : 1022;
    const int frexp_exponent = biased_exponent - exponent_offset;
    const u64 mantissa_bits = normalized_bits & 0x000fffffffffffffULL;
    const u64 half_mantissa_bits = mantissa_bits | 0x3fe0000000000000ULL;
    const double frexp_mantissa = __builtin_bit_cast(double, half_mantissa_bits);

    // b = m < 2/3; m = ldexp(m, b); e = frexp_exp(a) - b.
    const bool below_two_thirds = frexp_mantissa < 0x1.5555555555555p-1;
    const double doubled_mantissa = frexp_mantissa * 2.0;
    const double m = below_two_thirds ? doubled_mantissa : frexp_mantissa;
    const int b = below_two_thirds ? 1 : 0;
    const int e = frexp_exponent - b;

    // x = div(m - 1.0, fadd(1.0, m))
    const double numerator = m - 1.0;
    const double denominator_hi = 1.0 + m;
    const double denominator_hi_minus_one = denominator_hi - 1.0;
    const double denominator_lo = m - denominator_hi_minus_one;
    const double reciprocal = 1.0 / denominator_hi;
    const double q_hi = numerator * reciprocal;
    const double product_hi = q_hi * denominator_hi;
    const double negative_product_hi = -product_hi;
    const double product_lo_hi = fma(q_hi, denominator_hi, negative_product_hi);
    const double product_lo = fma(q_hi, denominator_lo, product_lo_hi);
    const double normalized_product_hi = product_hi + product_lo;
    const double product_hi_change = normalized_product_hi - product_hi;
    const double normalized_product_lo = product_lo - product_hi_change;
    const double difference_hi = numerator - normalized_product_hi;
    const double numerator_minus_difference = numerator - difference_hi;
    const double difference_lo_hi =
        numerator_minus_difference - normalized_product_hi;
    const double difference_lo = difference_lo_hi - normalized_product_lo;
    const double difference = difference_hi + difference_lo;
    const double q_lo = difference * reciprocal;
    const double x_hi = q_hi + q_lo;
    const double x_hi_change = x_hi - q_hi;
    const double x_lo = q_lo - x_hi_change;

    // s = sqr(x)
    const double square_hi = x_hi * x_hi;
    const double negative_square_hi = -square_hi;
    const double square_lo_hi = fma(x_hi, x_hi, negative_square_hi);
    const double twice_x_lo = 2.0 * x_lo;
    const double square_lo = fma(x_hi, twice_x_lo, square_lo_hi);
    const double s_hi = square_hi + square_lo;
    const double s_hi_change = s_hi - square_hi;
    const double s_lo = square_lo - s_hi_change;

    const double t = s_hi;
    const double p8 = fma(t, 0x1.dee674222de17p-4, 0x1.a6564968915a9p-4);
    const double p7 = fma(t, p8, 0x1.e25e43abe935ap-4);
    const double p6 = fma(t, p7, 0x1.110ef47e6c9c2p-3);
    const double p5 = fma(t, p6, 0x1.3b13bcfa74449p-3);
    const double p4 = fma(t, p5, 0x1.745d171bf3c30p-3);
    const double p3 = fma(t, p4, 0x1.c71c71c7792cep-3);
    const double p2 = fma(t, p3, 0x1.24924924920dap-2);
    const double p = fma(t, p2, 0x1.999999999999cp-2);

    // mul(s, p)
    const double sp_product_hi = s_hi * p;
    const double negative_sp_product_hi = -sp_product_hi;
    const double sp_product_lo_hi = fma(s_hi, p, negative_sp_product_hi);
    const double sp_product_lo = fma(s_lo, p, sp_product_lo_hi);
    const double sp_hi = sp_product_hi + sp_product_lo;
    const double sp_hi_change = sp_hi - sp_product_hi;
    const double sp_lo = sp_product_lo - sp_hi_change;

    // c = fadd(con(0x1.5555555555555p-1, 0x1.543b0d5df274dp-55), mul(s, p))
    const double c3_hi = 0x1.5555555555555p-1;
    const double c3_lo = 0x1.543b0d5df274dp-55;
    const double c_sum_hi = c3_hi + sp_hi;
    const double c_sum_hi_change = c_sum_hi - c3_hi;
    const double c_sum_lo_hi = sp_hi - c_sum_hi_change;
    const double c_low_sum = c3_lo + sp_lo;
    const double c_sum_lo = c_sum_lo_hi + c_low_sum;
    const double c_hi = c_sum_hi + c_sum_lo;
    const double c_hi_change = c_hi - c_sum_hi;
    const double c_lo = c_sum_lo - c_hi_change;

    // sx = mul(s, x)
    const double sx_product_hi = s_hi * x_hi;
    const double negative_sx_product_hi = -sx_product_hi;
    const double sx_product_lo_hi = fma(s_hi, x_hi, negative_sx_product_hi);
    const double sx_product_lo_mid = fma(s_hi, x_lo, sx_product_lo_hi);
    const double sx_product_lo = fma(s_lo, x_hi, sx_product_lo_mid);
    const double sx_hi = sx_product_hi + sx_product_lo;
    const double sx_hi_change = sx_hi - sx_product_hi;
    const double sx_lo = sx_product_lo - sx_hi_change;

    // tail = mul(sx, c)
    const double tail_product_hi = sx_hi * c_hi;
    const double negative_tail_product_hi = -tail_product_hi;
    const double tail_product_lo_hi =
        fma(sx_hi, c_hi, negative_tail_product_hi);
    const double tail_product_lo_mid = fma(sx_hi, c_lo, tail_product_lo_hi);
    const double tail_product_lo = fma(sx_lo, c_hi, tail_product_lo_mid);
    const double tail_hi = tail_product_hi + tail_product_lo;
    const double tail_hi_change = tail_hi - tail_product_hi;
    const double tail_lo = tail_product_lo - tail_hi_change;

    // r = fadd(ldx(x, 1), tail)
    const double twice_x_hi = x_hi * 2.0;
    const double twice_x_lo_ldx = x_lo * 2.0;
    const double r_sum_hi = twice_x_hi + tail_hi;
    const double r_sum_hi_change = r_sum_hi - twice_x_hi;
    const double r_sum_lo_hi = tail_hi - r_sum_hi_change;
    const double r_low_sum = twice_x_lo_ldx + tail_lo;
    const double r_sum_lo = r_sum_lo_hi + r_low_sum;
    const double r_hi = r_sum_hi + r_sum_lo;
    const double r_hi_change = r_hi - r_sum_hi;
    const double r_lo = r_sum_lo - r_hi_change;

    // l = mul(con(ln2_hi, ln2_lo), (double)e)
    const double ln2_hi = 0x1.62e42fefa39efp-1;
    const double ln2_lo = 0x1.abc9e3b39803fp-56;
    const double de = (double)e;
    const double scaled_hi = ln2_hi * de;
    const double negative_scaled_hi = -scaled_hi;
    const double scaled_lo_hi = fma(ln2_hi, de, negative_scaled_hi);
    const double scaled_lo = fma(ln2_lo, de, scaled_lo_hi);
    const double l_hi = scaled_hi + scaled_lo;
    const double l_hi_change = l_hi - scaled_hi;
    const double l_lo = scaled_lo - l_hi_change;

    // add(l, r)
    const double high_sum = l_hi + r_hi;
    const double high_r_part = high_sum - l_hi;
    const double high_l_part = high_sum - high_r_part;
    const double high_l_error = l_hi - high_l_part;
    const double high_r_error = r_hi - high_r_part;
    const double high_error = high_l_error + high_r_error;
    const double low_sum = l_lo + r_lo;
    const double low_r_part = low_sum - l_lo;
    const double low_l_part = low_sum - low_r_part;
    const double low_l_error = l_lo - low_l_part;
    const double low_r_error = r_lo - low_r_part;
    const double low_error = low_l_error + low_r_error;
    const double total_lo_hi = high_error + low_sum;
    const double total_hi = high_sum + total_lo_hi;
    const double total_hi_change = total_hi - high_sum;
    const double total_lo_lo = total_lo_hi - total_hi_change;
    const double total_lo = total_lo_lo + low_error;
    const double result_hi = total_hi + total_lo;
    const double result_hi_change = result_hi - total_hi;
    const double result_lo = total_lo - result_hi_change;
    const ep_double result = {result_hi, result_lo};
    return result;
}

// expepD.cl without EXTRA_ACCURACY: exp(x.hi + x.lo).
DEVICE double expep(double x_hi, double x_lo) {
    const double z = exp(x_hi);
    const double zz = fma(z, x_lo, z);
    const u64 z_bits = __builtin_bit_cast(u64, z);
    const u64 z_magnitude = z_bits & 0x7fffffffffffffffULL;
    const bool z_is_infinite = z_magnitude == 0x7ff0000000000000ULL;
    return z_is_infinite ? z : zz;
}

// powD_base.h is_integer
DEVICE bool is_integer(double ay) {
    const double truncated = trunc(ay);
    return truncated == ay;
}

// powD_base.h is_even_integer
DEVICE bool is_even_integer(double ay) {
    const double half = 0.5 * ay;
    return is_integer(half);
}

// powD_base.h is_odd_integer
DEVICE bool is_odd_integer(double ay) {
    const bool integer = is_integer(ay);
    const bool even = is_even_integer(ay);
    const bool not_even = !even;
    return integer ? not_even : integer;
}

// powD_base.h with COMPILING_POW
DEVICE double pow(double x_in, double y_in) {
    // if (x == 1.0) y = 1.0; if (y == 0.0) x = 1.0;
    const bool x_is_one = x_in == 1.0;
    const double y = x_is_one ? 1.0 : y_in;
    const bool y_is_zero = y == 0.0;
    const double x = y_is_zero ? 1.0 : x_in;

    const u64 x_bits = __builtin_bit_cast(u64, x);
    const u64 x_magnitude = x_bits & 0x7fffffffffffffffULL;
    const double ax = __builtin_bit_cast(double, x_magnitude);

    // expylnx = expep(omul(y, epln(ax)))
    const ep_double ln_ax = epln(ax);
    const double product_hi = y * ln_ax.hi;
    const double negative_product_hi = -product_hi;
    const double product_lo_hi = fma(y, ln_ax.hi, negative_product_hi);
    const double product_lo = fma(y, ln_ax.lo, product_lo_hi);
    // onrm
    const double sum = product_hi + product_lo;
    const double sum_change = sum - product_hi;
    const double sum_error = product_lo - sum_change;
    const u64 product_hi_bits = __builtin_bit_cast(u64, product_hi);
    const u64 product_hi_magnitude = product_hi_bits & 0x7fffffffffffffffULL;
    const bool product_hi_is_infinite =
        product_hi_magnitude == 0x7ff0000000000000ULL;
    const double ylnx_hi = product_hi_is_infinite ? product_hi : sum;
    const u64 ylnx_hi_bits = __builtin_bit_cast(u64, ylnx_hi);
    const u64 ylnx_hi_magnitude = ylnx_hi_bits & 0x7fffffffffffffffULL;
    const bool ylnx_hi_is_infinite =
        ylnx_hi_magnitude == 0x7ff0000000000000ULL;
    const double ylnx_lo = ylnx_hi_is_infinite ? 0.0 : sum_error;
    const double expylnx = expep(ylnx_hi, ylnx_lo);

    const bool is_odd_y = is_odd_integer(y);

    // ret = copysign(expylnx, is_odd_y ? x : 1.0)
    const double sign_source = is_odd_y ? x : 1.0;
    const double base_result = copysign(expylnx, sign_source);

    const double quiet_nan = __builtin_bit_cast(double, 0x7ff8000000000000ULL);
    const double infinity = __builtin_bit_cast(double, 0x7ff0000000000000ULL);

    // if (x < 0.0 && !is_integer(y)) ret = QNAN
    const bool x_negative = x < 0.0;
    const bool y_integer = is_integer(y);
    const bool y_not_integer = !y_integer;
    const bool negative_base_domain_error =
        x_negative ? y_not_integer : x_negative;
    const double with_domain = negative_base_domain_error ? quiet_nan
                                                           : base_result;

    // if (isinf(ay)) ret = ax == 1.0 ? ax :
    //     ((ax < 1.0) ^ y_is_neg_inf ? 0.0 : ay)
    const u64 y_bits = __builtin_bit_cast(u64, y);
    const u64 y_magnitude = y_bits & 0x7fffffffffffffffULL;
    const double ay = __builtin_bit_cast(double, y_magnitude);
    const bool ay_is_infinite = y_magnitude == 0x7ff0000000000000ULL;
    const bool y_is_neg_inf = y != ay;
    const bool ax_below_one = ax < 1.0;
    const bool decays = ax_below_one != y_is_neg_inf;
    const double infinite_y_away_from_one = decays ? 0.0 : ay;
    const bool ax_is_one = ax == 1.0;
    const double infinite_y_result =
        ax_is_one ? ax : infinite_y_away_from_one;
    const double with_infinite_y =
        ay_is_infinite ? infinite_y_result : with_domain;

    // if (isinf(ax) || x == 0.0) ret = copysign(
    //     (x == 0.0) ^ (y < 0.0) ? 0.0 : PINF, is_odd_y ? x : 0.0)
    const bool ax_is_infinite = x_magnitude == 0x7ff0000000000000ULL;
    const bool x_is_zero = x == 0.0;
    const bool y_negative = y < 0.0;
    const bool vanishes = x_is_zero != y_negative;
    const double zero_or_infinite_magnitude = vanishes ? 0.0 : infinity;
    const double zero_or_infinite_sign_source = is_odd_y ? x : 0.0;
    const double zero_or_infinite_result =
        copysign(zero_or_infinite_magnitude, zero_or_infinite_sign_source);
    const bool x_zero_or_infinite = ax_is_infinite ? ax_is_infinite : x_is_zero;
    const double with_zero_or_infinite_x =
        x_zero_or_infinite ? zero_or_infinite_result : with_infinite_y;

    // if (isunordered(x, y)) ret = QNAN
    const bool x_is_nan = x_magnitude > 0x7ff0000000000000ULL;
    const bool y_is_nan = y_magnitude > 0x7ff0000000000000ULL;
    const bool unordered = x_is_nan ? x_is_nan : y_is_nan;
    return unordered ? quiet_nan : with_zero_or_infinite_x;
}

// logD_base.h with COMPILING_LOG10, with the ep.h double-double helpers
// expanded in place. Identical to f64_math.H's log up to r = fadd(ldx(x, 1),
// s * x.hi * p).
DEVICE double log10(double a) {
    const u64 bits = __builtin_bit_cast(u64, a);
    const u64 magnitude = bits & 0x7fffffffffffffffULL;

    // m = frexp_mant(a) in [0.5, 1) and its exponent. Signed, zero, and
    // non-finite inputs produce unused values that the selects below replace.
    const bool subnormal = magnitude < 0x0010000000000000ULL;
    const double scaled_subnormal = a * 0x1p54;
    const double normalized_input = subnormal ? scaled_subnormal : a;
    const u64 normalized_bits = __builtin_bit_cast(u64, normalized_input);
    const u64 shifted_exponent = normalized_bits >> 52;
    const u64 exponent_bits = shifted_exponent & 2047;
    const int biased_exponent = (int)exponent_bits;
    const int exponent_offset = subnormal ? 1076 : 1022;
    const int frexp_exponent = biased_exponent - exponent_offset;
    const u64 mantissa_bits = normalized_bits & 0x000fffffffffffffULL;
    const u64 half_mantissa_bits = mantissa_bits | 0x3fe0000000000000ULL;
    const double frexp_mantissa = __builtin_bit_cast(double, half_mantissa_bits);

    // b = m < 2/3; m = ldexp(m, b); e = frexp_exp(a) - b.
    const bool below_two_thirds = frexp_mantissa < 0x1.5555555555555p-1;
    const double doubled_mantissa = frexp_mantissa * 2.0;
    const double m = below_two_thirds ? doubled_mantissa : frexp_mantissa;
    const int b = below_two_thirds ? 1 : 0;
    const int e = frexp_exponent - b;

    // x = div(m - 1.0, fadd(1.0, m))
    const double numerator = m - 1.0;
    const double denominator_hi = 1.0 + m;
    const double denominator_hi_minus_one = denominator_hi - 1.0;
    const double denominator_lo = m - denominator_hi_minus_one;
    const double reciprocal = 1.0 / denominator_hi;
    const double q_hi = numerator * reciprocal;
    const double product_hi = q_hi * denominator_hi;
    const double negative_product_hi = -product_hi;
    const double product_lo_hi = fma(q_hi, denominator_hi, negative_product_hi);
    const double product_lo = fma(q_hi, denominator_lo, product_lo_hi);
    const double normalized_product_hi = product_hi + product_lo;
    const double product_hi_change = normalized_product_hi - product_hi;
    const double normalized_product_lo = product_lo - product_hi_change;
    const double difference_hi = numerator - normalized_product_hi;
    const double numerator_minus_difference = numerator - difference_hi;
    const double difference_lo_hi =
        numerator_minus_difference - normalized_product_hi;
    const double difference_lo = difference_lo_hi - normalized_product_lo;
    const double difference = difference_hi + difference_lo;
    const double q_lo = difference * reciprocal;
    const double x_hi = q_hi + q_lo;
    const double x_hi_change = x_hi - q_hi;
    const double x_lo = q_lo - x_hi_change;

    const double s = x_hi * x_hi;
    const double p5 = fma(s, 0x1.3ab76bf559e2bp-3, 0x1.385386b47b09ap-3);
    const double p4 = fma(s, p5, 0x1.7474dd7f4df2ep-3);
    const double p3 = fma(s, p4, 0x1.c71c016291751p-3);
    const double p2 = fma(s, p3, 0x1.249249b27acf1p-2);
    const double p1 = fma(s, p2, 0x1.99999998ef7b6p-2);
    const double p = fma(s, p1, 0x1.5555555555780p-1);

    // r = fadd(ldx(x, 1), s * x.hi * p)
    const double twice_x_hi = x_hi * 2.0;
    const double twice_x_lo = x_lo * 2.0;
    const double s_x_hi = s * x_hi;
    const double tail = s_x_hi * p;
    const double sum_hi = twice_x_hi + tail;
    const double sum_hi_change = sum_hi - twice_x_hi;
    const double sum_lo_hi = tail - sum_hi_change;
    const double sum_lo = sum_lo_hi + twice_x_lo;
    const double r_hi = sum_hi + sum_lo;
    const double r_hi_change = r_hi - sum_hi;
    const double r_lo = sum_lo - r_hi_change;

    // l = mul(con(log10(2)_hi, log10(2)_lo), (double)e)
    const double log10_2_hi = 0x1.34413509f79ffp-2;
    const double log10_2_lo = -0x1.9dc1da994fd21p-59;
    const double de = (double)e;
    const double scaled_hi = log10_2_hi * de;
    const double negative_scaled_hi = -scaled_hi;
    const double scaled_lo_hi = fma(log10_2_hi, de, negative_scaled_hi);
    const double scaled_lo = fma(log10_2_lo, de, scaled_lo_hi);
    const double l_hi = scaled_hi + scaled_lo;
    const double l_hi_change = l_hi - scaled_hi;
    const double l_lo = scaled_lo - l_hi_change;

    // v = mul(con(log10(e)_hi, log10(e)_lo), r)
    const double log10_e_hi = 0x1.bcb7b1526e50ep-2;
    const double log10_e_lo = 0x1.95355baaafad3p-57;
    const double v_product_hi = log10_e_hi * r_hi;
    const double negative_v_product_hi = -v_product_hi;
    const double v_product_lo_hi =
        fma(log10_e_hi, r_hi, negative_v_product_hi);
    const double v_product_lo_mid = fma(log10_e_hi, r_lo, v_product_lo_hi);
    const double v_product_lo = fma(log10_e_lo, r_hi, v_product_lo_mid);
    const double v_hi = v_product_hi + v_product_lo;
    const double v_hi_change = v_hi - v_product_hi;
    const double v_lo = v_product_lo - v_hi_change;

    // add(l, v)
    const double high_sum = l_hi + v_hi;
    const double high_v_part = high_sum - l_hi;
    const double high_l_part = high_sum - high_v_part;
    const double high_l_error = l_hi - high_l_part;
    const double high_v_error = v_hi - high_v_part;
    const double high_error = high_l_error + high_v_error;
    const double low_sum = l_lo + v_lo;
    const double low_v_part = low_sum - l_lo;
    const double low_l_part = low_sum - low_v_part;
    const double low_l_error = l_lo - low_l_part;
    const double low_v_error = v_lo - low_v_part;
    const double low_error = low_l_error + low_v_error;
    const double total_lo_hi = high_error + low_sum;
    const double total_hi = high_sum + total_lo_hi;
    const double total_hi_change = total_hi - high_sum;
    const double total_lo_lo = total_lo_hi - total_hi_change;
    const double total_lo = total_lo_lo + low_error;
    const double result = total_hi + total_lo;

    // OCML's NaN result comes from propagation through frexp. Setting the
    // quiet bit produces the same value without it.
    const bool is_infinite = magnitude == 0x7ff0000000000000ULL;
    const bool is_nan = magnitude > 0x7ff0000000000000ULL;
    const bool negative = a < 0.0;
    const bool is_zero = a == 0.0;
    const double negative_infinity =
        __builtin_bit_cast(double, 0xfff0000000000000ULL);
    const double canonical_nan =
        __builtin_bit_cast(double, 0x7ff8000000000000ULL);
    const double quiet_nan =
        __builtin_bit_cast(double, bits | 0x0008000000000000ULL);
    const double with_infinity = is_infinite ? a : result;
    const double with_negative = negative ? canonical_nan : with_infinity;
    const double with_zero = is_zero ? negative_infinity : with_negative;
    return is_nan ? quiet_nan : with_zero;
}

// epexpepD.cl: exp(x.hi + x.lo) as a double-double, for 0 <= x.hi <= 32.
DEVICE ep_double epexpep(double x_hi, double x_lo) {
    // dn = rint(x.hi * log2(e)) with the shift of f64_math.H's exp; dn is in
    // [0, 47], and its low 32 bits encode (int)dn.
    const double shift = 0x1.8p52;
    const double scaled_input = x_hi * 0x1.71547652b82fep+0;
    const double shifted = scaled_input + shift;
    const double dn = shifted - shift;
    const int n = (int)__builtin_bit_cast(u64, shifted);

    // t = fsub(fsub(fadd(MAD(dn, -ln2_1, x.hi), x.lo), dn * ln2_2),
    //          dn * ln2_3)
    const double reduced = fma(dn, -0x1.62e42fefa3000p-1, x_hi);
    const double a_hi = reduced + x_lo;
    const double a_hi_change = a_hi - reduced;
    const double a_lo = x_lo - a_hi_change;
    const double ln2_2_part = dn * 0x1.3de6af278e000p-42;
    const double b_difference = a_hi - ln2_2_part;
    const double b_kept = a_hi - b_difference;
    const double b_difference_lo = b_kept - ln2_2_part;
    const double b_unnormalized_lo = b_difference_lo + a_lo;
    const double b_hi = b_difference + b_unnormalized_lo;
    const double b_hi_change = b_hi - b_difference;
    const double b_lo = b_unnormalized_lo - b_hi_change;
    const double ln2_3_part = dn * 0x1.9cc01f97b57a0p-83;
    const double c_difference = b_hi - ln2_3_part;
    const double c_kept = b_hi - c_difference;
    const double c_difference_lo = c_kept - ln2_3_part;
    const double c_unnormalized_lo = c_difference_lo + b_lo;
    const double t_hi = c_difference + c_unnormalized_lo;
    const double t_hi_change = t_hi - c_difference;
    const double t_lo = c_unnormalized_lo - t_hi_change;

    const double th = t_hi;
    const double p9 = fma(th, 0x1.ade156a5dcb37p-26, 0x1.28af3fca7ab0cp-22);
    const double p8 = fma(th, p9, 0x1.71dee623fde64p-19);
    const double p7 = fma(th, p8, 0x1.a01997c89e6b0p-16);
    const double p6 = fma(th, p7, 0x1.a01a014761f6ep-13);
    const double p5 = fma(th, p6, 0x1.6c16c1852b7b0p-10);
    const double p4 = fma(th, p5, 0x1.1111111122322p-7);
    const double p3 = fma(th, p4, 0x1.55555555502a1p-5);
    const double p2 = fma(th, p3, 0x1.5555555555511p-3);
    const double p = fma(th, p2, 0x1.000000000000bp-1);

    // sqr(t)
    const double square_hi = t_hi * t_hi;
    const double negative_square_hi = -square_hi;
    const double square_lo_hi = fma(t_hi, t_hi, negative_square_hi);
    const double twice_t_lo = 2.0 * t_lo;
    const double square_lo = fma(t_hi, twice_t_lo, square_lo_hi);
    const double sq_hi = square_hi + square_lo;
    const double sq_hi_change = sq_hi - square_hi;
    const double sq_lo = square_lo - sq_hi_change;

    // mul(sqr(t), p)
    const double m_product_hi = sq_hi * p;
    const double negative_m_product_hi = -m_product_hi;
    const double m_product_lo_hi = fma(sq_hi, p, negative_m_product_hi);
    const double m_product_lo = fma(sq_lo, p, m_product_lo_hi);
    const double m_hi = m_product_hi + m_product_lo;
    const double m_hi_change = m_hi - m_product_hi;
    const double m_lo = m_product_lo - m_hi_change;

    // u = fadd(t, mul(sqr(t), p))
    const double u_sum_hi = t_hi + m_hi;
    const double u_sum_hi_change = u_sum_hi - t_hi;
    const double u_sum_lo_hi = m_hi - u_sum_hi_change;
    const double u_low_sum = t_lo + m_lo;
    const double u_sum_lo = u_sum_lo_hi + u_low_sum;
    const double u_hi = u_sum_hi + u_sum_lo;
    const double u_hi_change = u_hi - u_sum_hi;
    const double u_lo = u_sum_lo - u_hi_change;

    // r = fadd(1.0, u)
    const double r_sum_hi = 1.0 + u_hi;
    const double r_sum_hi_change = r_sum_hi - 1.0;
    const double r_sum_lo_hi = u_hi - r_sum_hi_change;
    const double r_sum_lo = r_sum_lo_hi + u_lo;
    const double r_hi = r_sum_hi + r_sum_lo;
    const double r_hi_change = r_hi - r_sum_hi;
    const double r_lo = r_sum_lo - r_hi_change;

    // ldx(r, (int)dn): 2^dn is normal and both products are exact.
    const int scale_exponent = n + 1023;
    const u64 scale_bits = (u64)scale_exponent << 52;
    const double scale = __builtin_bit_cast(double, scale_bits);
    const double result_hi = r_hi * scale;
    const double result_lo = r_lo * scale;
    const ep_double result = {result_hi, result_lo};
    return result;
}

// tanhD.cl
DEVICE double tanh(double x) {
    const u64 bits = __builtin_bit_cast(u64, x);
    const u64 magnitude = bits & 0x7fffffffffffffffULL;
    const double y = __builtin_bit_cast(double, magnitude);

    // Every y > 19.0625 returns 1.0 below, so clamping y to 32 changes no
    // returned value and keeps dn in epexpep within [0, 47]. NaN is selected
    // at the end.
    const bool y_in_range = y <= 32.0;
    const double ey = y_in_range ? y : 32.0;

    // e = epexpep(con(y, 0.0))
    const ep_double e = epexpep(ey, 0.0);

    // ei = rcp(e)
    const double qhi = 1.0 / e.hi;
    // r = fsub(1.0, mul(qhi, e))
    const double a_product_hi = qhi * e.hi;
    const double negative_a_product_hi = -a_product_hi;
    const double a_product_lo_hi = fma(qhi, e.hi, negative_a_product_hi);
    const double a_product_lo = fma(qhi, e.lo, a_product_lo_hi);
    const double a_hi = a_product_hi + a_product_lo;
    const double a_hi_change = a_hi - a_product_hi;
    const double a_lo = a_product_lo - a_hi_change;
    const double r1_difference = 1.0 - a_hi;
    const double r1_kept = 1.0 - r1_difference;
    const double r1_difference_lo = r1_kept - a_hi;
    const double r1_unnormalized_lo = r1_difference_lo - a_lo;
    const double r1_hi = r1_difference + r1_unnormalized_lo;
    const double r1_hi_change = r1_hi - r1_difference;
    const double r1_lo = r1_unnormalized_lo - r1_hi_change;
    const double qmi = r1_hi * qhi;
    // r = fsub(r, mul(qmi, e))
    const double b_product_hi = qmi * e.hi;
    const double negative_b_product_hi = -b_product_hi;
    const double b_product_lo_hi = fma(qmi, e.hi, negative_b_product_hi);
    const double b_product_lo = fma(qmi, e.lo, b_product_lo_hi);
    const double b_hi = b_product_hi + b_product_lo;
    const double b_hi_change = b_hi - b_product_hi;
    const double b_lo = b_product_lo - b_hi_change;
    const double r2_difference = r1_hi - b_hi;
    const double r2_kept = r1_hi - r2_difference;
    const double r2_difference_lo = r2_kept - b_hi;
    const double r2_lo_with_a = r2_difference_lo + r1_lo;
    const double r2_unnormalized_lo = r2_lo_with_a - b_lo;
    const double r2_hi = r2_difference + r2_unnormalized_lo;
    const double qlo = r2_hi * qhi;
    // q = fadd(qhi, qmi); q.lo += qlo; nrm(q)
    const double q_sum = qhi + qmi;
    const double q_sum_change = q_sum - qhi;
    const double q_sum_lo = qmi - q_sum_change;
    const double q_lo_total = q_sum_lo + qlo;
    const double ei_hi = q_sum + q_lo_total;
    const double ei_hi_change = ei_hi - q_sum;
    const double ei_lo = q_lo_total - ei_hi_change;

    // n = fsub(e, ei)
    const double n_difference = e.hi - ei_hi;
    const double n_kept = e.hi - n_difference;
    const double n_difference_lo = n_kept - ei_hi;
    const double n_lo_with_a = n_difference_lo + e.lo;
    const double n_unnormalized_lo = n_lo_with_a - ei_lo;
    const double n_hi = n_difference + n_unnormalized_lo;
    const double n_hi_change = n_hi - n_difference;
    const double n_lo = n_unnormalized_lo - n_hi_change;

    // d = fadd(e, ei)
    const double d_sum = e.hi + ei_hi;
    const double d_sum_change = d_sum - e.hi;
    const double d_sum_lo_hi = ei_hi - d_sum_change;
    const double d_low_sum = e.lo + ei_lo;
    const double d_sum_lo = d_sum_lo_hi + d_low_sum;
    const double d_hi = d_sum + d_sum_lo;
    const double d_hi_change = d_hi - d_sum;
    const double d_lo = d_sum_lo - d_hi_change;

    // t = fdiv(n, d)
    const double reciprocal = 1.0 / d_hi;
    const double t_qhi = n_hi * reciprocal;
    const double c_product_hi = t_qhi * d_hi;
    const double negative_c_product_hi = -c_product_hi;
    const double c_product_lo_hi = fma(t_qhi, d_hi, negative_c_product_hi);
    const double c_product_lo = fma(t_qhi, d_lo, c_product_lo_hi);
    const double c_hi = c_product_hi + c_product_lo;
    const double c_hi_change = c_hi - c_product_hi;
    const double c_lo = c_product_lo - c_hi_change;
    const double f_difference = n_hi - c_hi;
    const double f_kept = n_hi - f_difference;
    const double f_difference_lo = f_kept - c_hi;
    const double f_lo_minus_product = f_difference_lo - c_lo;
    const double f_lo = f_lo_minus_product + n_lo;
    const double f_sum = f_difference + f_lo;
    const double t_qlo = f_sum * reciprocal;
    const double z = t_qhi + t_qlo;

    // z = y > 19.0625 ? 1.0 : z; z = y < 0x1.0p-27 ? y : z;
    const bool saturated = y > 19.0625;
    const double with_saturation = saturated ? 1.0 : z;
    const bool tiny = y < 0x1.0p-27;
    const double with_tiny = tiny ? y : with_saturation;
    const double result = copysign(with_tiny, x);

    // OCML's NaN result comes from propagation through epexpep and keeps the
    // sign of x. Setting the quiet bit produces the same value.
    const bool is_nan = magnitude > 0x7ff0000000000000ULL;
    const double quiet_nan =
        __builtin_bit_cast(double, bits | 0x0008000000000000ULL);
    return is_nan ? quiet_nan : result;
}
} // namespace f64_math
} // namespace chemistry
