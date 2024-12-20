#ifndef DIAGNOSTICVARIABLES_HPP_
#define DIAGNOSTICVARIABLES_HPP_

// assign an enum to each variable
enum
{
    c_Ham,

    c_Mom1,
    c_Mom2,
    c_Mom3,
    c_Mom,

    c_rho,

    c_S1,
    c_S2,
    c_S3,

    c_Ham_abs,
    c_Mom1_abs,
    c_Mom2_abs,
    c_Mom3_abs,
    c_Mom_abs,

    c_Ham_norm,
    c_Mom_norm,

    NUM_DIAGNOSTIC_VARS
};

namespace DiagnosticVariables
{
static const std::array<std::string, NUM_DIAGNOSTIC_VARS> variable_names = {
    "Ham",

    "Mom1",     "Mom2",     "Mom3",     "Mom",

    "rho",

    "S1",       "S2",       "S3",

    "Ham_abs",  "Mom1_abs", "Mom2_abs", "Mom3_abs", "Mom_abs",

    "Ham_norm", "Mom_norm"};
}

#endif /* DIAGNOSTICVARIABLES_HPP_ */