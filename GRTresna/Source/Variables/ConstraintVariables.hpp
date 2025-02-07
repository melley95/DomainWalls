/* GRTresna
 * Copyright 2024 The GRTL Collaboration.
 * Please refer to LICENSE in GRTresna's root directory.
 */

#ifndef CONSTRAINTVARIABLES_HPP
#define CONSTRAINTVARIABLES_HPP

#include "ArrayTools.hpp"
#include "ParityDefinitions.hpp"

// assign an enum to each constraint variable
enum
{
    c_psi,

    c_V1,
    c_V2,
    c_V3,
    c_U,

    NUM_CONSTRAINT_VARS
};

namespace ConstraintVariables
{
static const std::array<std::string, NUM_CONSTRAINT_VARS> variable_names = {
    "psi", "V1", "V2", "V3", "U"};

static constexpr std::array<int, NUM_CONSTRAINT_VARS> const vars_parity = {
    EVEN, ODD_X, ODD_Y, ODD_Z, EVEN};

} // namespace ConstraintVariables

#endif /* CONSTRAINTVARIABLES_HPP */