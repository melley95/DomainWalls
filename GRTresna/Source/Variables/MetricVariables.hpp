/* GRTresna
 * Copyright 2024 The GRTL Collaboration.
 * Please refer to LICENSE in GRTresna's root directory.
 */

#ifndef METRICVARIABLES_HPP
#define METRICVARIABLES_HPP

#include "ArrayTools.hpp"
#include "ParityDefinitions.hpp"

// assign an enum to each metric variable
enum
{
    c_psi_reg,

    c_V1_0,
    c_V2_0,
    c_V3_0,
    c_U_0,

    c_K_0,

    c_A11_0,
    c_A12_0,
    c_A13_0,
    c_A22_0,
    c_A23_0,
    c_A33_0,

    NUM_METRIC_VARS
};

namespace MetricVariables
{
static const std::array<std::string, NUM_METRIC_VARS> variable_names = {
    "psi_reg",

    "V1_0",    "V2_0",  "V3_0",  "U_0",

    "K_0",

    "A11_0",   "A12_0", "A13_0", "A22_0", "A23_0", "A33_0"};

static constexpr std::array<int, NUM_METRIC_VARS> const vars_parity = {
    EVEN, ODD_X,  ODD_Y,  ODD_Z, EVEN,   EVEN,
    EVEN, ODD_XY, ODD_XZ, EVEN,  ODD_YZ, EVEN};

} // namespace MetricVariables

#endif /* METRICVARIABLES_HPP */