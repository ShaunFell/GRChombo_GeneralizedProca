/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef USERVARIABLES_HPP
#define USERVARIABLES_HPP

#include "DiagnosticVariables.hpp"

/// This enum gives the index of every variable stored in the grid
enum
{

    // Note that it is important that the first enum value is set to 1 more than
    c_phi, //scalar part of proca field

    c_Avec1, //spatial part of proca field
    c_Avec2,
    c_Avec3,

    c_Evec1, //conjugate momentum of proca field
    c_Evec2,
    c_Evec3,

    c_Z, //constraint violation damping scalar field

    //metric field
    c_chi,

    c_h11,
    c_h12,
    c_h13,
    c_h22,
    c_h23,
    c_h33,

    c_K,

    c_A11,
    c_A12,
    c_A13,
    c_A22,
    c_A23,
    c_A33,

    c_lapse,

    c_shift1,
    c_shift2,
    c_shift3,


    NUM_VARS
};

namespace UserVariables
{
static const std::array<std::string, NUM_VARS> variable_names = {
    "phi",

    "Avec1", "Avec2", "Avec3",

    "Evec1", "Evec2", "Evec3",

    "Z",
    
    "chi",

    "h11",    "h12",    "h13",    "h22", "h23", "h33",

    "K",

    "A11",    "A12",    "A13",    "A22", "A23", "A33",

    "lapse",

    "shift1", "shift2", "shift3"
    };

} // namespace UserVariables

#include "UserVariables.inc.hpp"

#endif /* USERVARIABLES_HPP */
