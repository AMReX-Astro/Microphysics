#!/usr/bin/env python3

import os
import warnings
import math
from scipy import constants

'''
This scripts writes out fundamental_constants.H
where most constants are from scipy.constants
'''

# Find PATH to write the file
MICROPHYSICS_HOME = os.environ.get("MICROPHYSICS_HOME")
if MICROPHYSICS_HOME is None:
    warnings.warn("MICROPHYSICS_HOME is not set. Writing to current directory.", stacklevel=2)
    fname = "fundamental_constants.H"
else:
    fname = os.path.join(MICROPHYSICS_HOME, "constants/fundamental_constants.H")


text = \
f"""
#ifndef FUNDAMENTAL_CONSTANTS_H
#define FUNDAMENTAL_CONSTANTS_H

#include <AMReX.H>
#include <AMReX_REAL.H>
#include <AMReX_Array.H>

//
// *** THIS FILE IS AUTO GENERATED VIA write_fundamental_constants.py ***
//
// Non-legacy constants are obtained from scipy.constants.
// Legacy constants are out-dated constants used for aprox networks only.
//

namespace C
{{
    // speed of light in vacuum
    constexpr amrex::Real c_light = {constants.value("speed of light in vacuum") / constants.centi};  // cm/s

    // newton's gravitational constant
    constexpr amrex::Real Gconst = 6.67428e-8;  // cm^3/g/s^2

    // new value; if uncommented initial models will need to be re-HSE'ed
    // constexpr amrex::Real Gconst = 6.67384e-8;  // cm^3/g/s^2

    // boltzmann's constant
    constexpr amrex::Real k_B = {constants.value("Boltzmann constant") / constants.erg};  // erg/K

    // planck's constant over 2pi
    constexpr amrex::Real hbar = {constants.value("Planck constant") / (2.0 * math.pi * constants.erg)};  // erg

    // planck's constant
    constexpr amrex::Real hplanck = {constants.value("Planck constant") / constants.erg};  // erg s

    // avogradro's Number
    constexpr amrex::Real n_A = {constants.value("Avogadro constant")};  // mol^-1

    // convert eV to erg
    constexpr amrex::Real ev2erg = {constants.eV / constants.erg};

    // convert MeV to eV
    constexpr amrex::Real MeV2eV = 1.0e6;

    // convert MeV to erg
    constexpr amrex::Real MeV2erg = {constants.mega * constants.eV / constants.erg};

    // convert MeV to grams
    constexpr amrex::Real MeV2gr  = (MeV2eV * ev2erg) / (c_light * c_light);

    // conversion factor for nuclear energy generation rate
    constexpr amrex::Real enuc_conv2 = -n_A * c_light * c_light;

    // mass of proton
    constexpr amrex::Real m_p = {constants.value("proton mass") * constants.kilo};  // g

    // mass of neutron
    constexpr amrex::Real m_n = {constants.value("neutron mass") * constants.kilo};  // g

    // mass of electron
    constexpr amrex::Real m_e = {constants.value("electron mass") * constants.kilo};  // g

    // atomic mass unit
    constexpr amrex::Real m_u = {constants.value("unified atomic mass unit") * constants.kilo}; // g

    // electron charge
    // NIST: q_e = 1.602176565e-19 C
    //
    // C is the SI unit Coulomb; in cgs we have the definition:
    //     1 C = 0.1 * |c_light| * 1 statC
    // where statC is the cgs unit statCoulomb; 1 statC = 1 erg^1/2 cm^1/2
    // and |c_light| is the speed of light in cgs (but without units)
    constexpr amrex::Real q_e = {constants.value("elementary charge") * (constants.c * 100) / 10};  // erg^1/2 cm^1/2

    // stefan-boltzmann constant
    constexpr amrex::Real sigma_SB = {constants.value("Stefan-Boltzmann constant") * constants.centi**2 / constants.erg};  // erg/s/cm^2/K^4

    // radiation constant
    constexpr amrex::Real a_rad = 4.0*sigma_SB/c_light;

    // Number of centimeters in a parsec and an AU.
    // Note that since the length of an AU is defined exactly
    // by IAU convention, the length of a parsec is also
    // defined exactly as (6.48e5 / pi) AU.
    constexpr amrex::Real AU = {constants.au / constants.centi};  // cm
    constexpr amrex::Real parsec = {constants.parsec / constants.centi};  // cm

    // Hubble constant (in s^{{-1}}, converted from 100 (km/s)/Mpc by dividing by 3.08568025e19km/Mpc)
    constexpr amrex::Real Hubble_const = 32.407764868e-19;

    // solar mass (from http://asa.usno.navy.mil/SecK/Constants.html)
    constexpr amrex::Real M_solar = 1.9884e33;

    // solar radius
    constexpr amrex::Real R_solar = 6.957e10;

    namespace Legacy
    {{
        // These are the values of the constants used in the original aprox networks
        constexpr amrex::Real m_n = 1.67492721184e-24;
        constexpr amrex::Real m_p = 1.67262163783e-24;
        constexpr amrex::Real m_e = 9.1093821545e-28;

        constexpr amrex::Real eV2erg  = 1.60217648740e-12;
        constexpr amrex::Real MeV2erg = eV2erg * 1.0e6;
        constexpr amrex::Real MeV2gr  = MeV2erg / (c_light * c_light);

        constexpr amrex::Real n_A = 6.0221417930e23;

        // conversion factor for nuclear energy generation rate
        constexpr amrex::Real enuc_conv2 = -n_A * c_light * c_light;
    }}
}}

#endif
"""


with open(fname, "w") as f:
    f.write(text)
