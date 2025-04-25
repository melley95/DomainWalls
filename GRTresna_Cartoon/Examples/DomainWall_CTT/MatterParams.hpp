#ifndef MATTERPARAMS_HPP_
#define MATTERPARAMS_HPP_

#include "REAL.H"

namespace MatterParams
{

struct params_t
{
    Real lambda;
    Real eta;
    Real b;
    Real e;

 
};

inline void read_params(GRParmParse &pp, params_t &matter_params)
{
    pp.get("eta", matter_params.eta);
    pp.get("lambda", matter_params.lambda);
    pp.get("R0", matter_params.b);
    pp.get("e", matter_params.e);

}

}; // namespace MatterParams

#endif