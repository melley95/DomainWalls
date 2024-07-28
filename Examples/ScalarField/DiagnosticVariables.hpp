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

    c_rho_grad,
    c_rho_kin,
    c_rho_pot,

    c_Mom,



    NUM_DIAGNOSTIC_VARS
};

namespace DiagnosticVariables
{
static const std::array<std::string, NUM_DIAGNOSTIC_VARS> variable_names = {
    "Ham",

    "rho_grad", "rho_kin", "rho_pot",

    "Mom"};
}

#endif /* DIAGNOSTICVARIABLES_HPP */
