/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef DIAGNOSTICVARIABLES_HPP
#define DIAGNOSTICVARIABLES_HPP

// assign an enum to each variable
enum
{
    c_Ham,

    c_Ham_abs_sum,

    c_Mom1,
    c_Mom2,
    c_Mom3,

    c_Mom_abs_sum1,
    c_Mom_abs_sum2,
    c_Mom_abs_sum3,

    c_Weyl4_Re,
    c_Weyl4_Im,

    c_gauss,
    c_Asquared,
    c_gnn,
    c_g,
    c_Z_out,


    NUM_DIAGNOSTIC_VARS
};

namespace DiagnosticVariables
{
static const std::array<std::string, NUM_DIAGNOSTIC_VARS> variable_names = {
    "Ham",
    
    "Ham_abs_sum",

    "Mom1", "Mom2", "Mom3",

    "Mom_abs_sum1", "Mom_abs_sum2", "Mom_abs_sum3",

    "Weyl4_Re", "Weyl4_Im",
    
    "Gauss", "Asquared", "gnn", "g", "Z"};
}

#endif /* DIAGNOSTICVARIABLES_HPP */
