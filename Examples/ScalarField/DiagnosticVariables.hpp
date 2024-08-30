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
    c_Mom,
  



    c_rho_grad,
    c_rho_kin,
    c_rho_pot,



    c_dxphi,
    c_Pi_out,

    c_rhoLL, // the Landau Lifshitz rho
    c_source,
    c_flux,







    NUM_DIAGNOSTIC_VARS
};

namespace DiagnosticVariables
{
static const std::array<std::string, NUM_DIAGNOSTIC_VARS> variable_names = {
    "Ham", "Ham_abs_sum", "Mom", 

    "rho_grad", "rho_kin", "rho_pot",

    "dxphi", "Pi",

    "rhoLL",    "source",    "flux"};
}

#endif /* DIAGNOSTICVARIABLES_HPP */
