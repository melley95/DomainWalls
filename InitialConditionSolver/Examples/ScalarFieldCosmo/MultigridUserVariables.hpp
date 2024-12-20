#ifndef USERVARIABLES_HPP
#define USERVARIABLES_HPP

#include "CoreVariables.hpp"
#include "ScalarFieldUserVariables.hpp"

namespace MultigridUserVariables
{
static const std::array<std::string, NUM_METRIC_VARS>
    metric_variable_names = MetricVariables::variable_names;
static const std::array<std::string, NUM_MULTIGRID_VARS - NUM_METRIC_VARS>
    matter_variable_names = MatterVariables::variable_names;
} // namespace MultigridUserVariables

#endif /* USERVARIABLES_HPP */