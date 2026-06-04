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

    c_Mom1,
    c_Mom2,

    // sqrt(Mom_1^2 + Mom_2^2)
    c_Mom,

    c_Weyl4_Re,
    c_Weyl4_Im,


    c_rho,

    c_Sx,
    c_Sy,

    c_rhoLL,
    c_source,

    c_Sxx,
    c_Sxy,
    c_Syy,

    c_Sww,

    c_det_gamma,


    c_tr_A2,

    c_ricci_scalar,

    c_Madm,

   

    c_Px_adm,
    c_Py_adm,
    c_Pz_adm,

    // c_Sxx_TF,
    // c_Sxy_TF,
    // c_Syy_TF,
    // c_Sww_TF,

    NUM_DIAGNOSTIC_VARS
};

namespace DiagnosticVariables
{
static const std::array<std::string, NUM_DIAGNOSTIC_VARS> variable_names = {
    "Ham",

    "Mom1",     "Mom2",     "Mom",

    "Weyl4_Re", "Weyl4_Im",

 //   "M_adm",    "P_adm",

    "rho",  "Sx",  "Sy",

     "rhoLL", "source", 
     
     "Sxx",    "Sxy",    "Syy",
     
     "Sww",

     "det_gamma",


     "tr_A2",

     "ricci_scalar",

     "M_adm", 

     "Px_adm", "Py_adm", "Pz_adm"
    
    
    };     
// "S",        "Sxx_TF", "Sxy_TF", "Syy_TF", "Sww_TF"
} // namespace DiagnosticVariables

#endif /* DIAGNOSTICVARIABLES_HPP */
