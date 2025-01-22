/* GRTresna
 * Copyright 2024 The GRTL Collaboration.
 * Please refer to LICENSE in GRTresna's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP_
#define SIMULATIONPARAMETERS_HPP_

#include "BoundaryConditions.hpp"
#include "GRParmParse.hpp"
#include "Grids.hpp"
#include "PsiAndAijFunctions.hpp"

struct base_params;

template <class method_t, class matter_t> class SimulationParameters
{
  public:
    SimulationParameters() {}

    SimulationParameters(GRParmParse &pp) { read_params(pp); }

    void read_params(GRParmParse &pp)
    {
        read_base_params(pp);
        method_t::read_params(pp, method_params);
        matter_t::read_params(pp, matter_params);
        Grids::read_params(pp, grid_params);
        PsiAndAijFunctions::read_params(pp, psi_and_Aij_params);
    }

    void read_base_params(GRParmParse &pp);

    struct BaseParams;

    typename method_t::params_t method_params;
    typename matter_t::params_t matter_params;
    typename BoundaryConditions::params_t boundary_params;
    typename Grids::params_t grid_params;
    typename PsiAndAijFunctions::params_t psi_and_Aij_params;
    BaseParams base_params;

  private:
};

#include "SimulationParameters.impl.hpp"

#endif /* SIMULATIONPARAMETERS_HPP_ */