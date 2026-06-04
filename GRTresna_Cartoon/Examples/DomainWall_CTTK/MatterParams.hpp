#ifndef MATTERPARAMS_HPP_
#define MATTERPARAMS_HPP_

#include "REAL.H"

namespace MatterParams
{

struct params_t
{
    Real lambda;
    Real eta;
    Real a;
    Real b;


 
};

inline void read_params(GRParmParse &pp, params_t &matter_params)
{
    pp.get("eta", matter_params.eta);
    pp.get("lambda", matter_params.lambda);
    pp.get("a", matter_params.a);
    pp.get("b", matter_params.b);


}

}; // namespace MatterParams

#endif