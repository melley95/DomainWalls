#ifndef SCALARFIELDUSERVARIABLES_HPP_
#define SCALARFIELDUSERVARIABLES_HPP_

#include "CoreVariables.hpp"

// Matter Vars
enum
{
    c_phi_0 = NUM_METRIC_VARS,
    c_Pi_0,
    NUM_MULTIGRID_VARS
};

namespace MatterVariables
{

static const std::array<std::string, NUM_MULTIGRID_VARS - NUM_METRIC_VARS>
    variable_names = {"phi_0", "Pi_0"};

static constexpr std::array<int, NUM_MULTIGRID_VARS - NUM_METRIC_VARS> const
    vars_parity = {0, 0};

} // namespace MatterVariables

#endif // SCALARFIELDUSERVARIABLES_HPP_