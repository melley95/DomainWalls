#ifndef SIMULATIONPARAMETERS_HPP_
#define SIMULATIONPARAMETERS_HPP_

#include "BoundaryConditions.hpp"
#include "Metric.hpp"
#include "GRParmParse.hpp"
#include "Grids.hpp"

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
        Metric::read_params(pp, metric_params);
    }

    void read_base_params(GRParmParse &pp);

    struct BaseParams;

    typename method_t::params_t method_params;
    typename matter_t::params_t matter_params;
    typename BoundaryConditions::params_t boundary_params;
    typename Grids::params_t grid_params;
    typename Metric::params_t metric_params;
    BaseParams base_params;

  private:
};

#include "SimulationParameters.impl.hpp"

#endif /* SIMULATIONPARAMETERS_HPP_ */