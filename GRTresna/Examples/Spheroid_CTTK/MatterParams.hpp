#ifndef MATTERPARAMS_HPP_
#define MATTERPARAMS_HPP_
#include "REAL.H"
namespace MatterParams
{
struct params_t
{
    Real lambda;
    Real eta;
    Real R0;
    Real eps1;
 
};
inline void read_params(GRParmParse &pp, params_t &matter_params)
{
    pp.get("eta", matter_params.eta);
    pp.get("lambda", matter_params.lambda);
    pp.get("R0", matter_params.R0);
    pp.get("eps1", matter_params.eps1);
}
}; // namespace MatterParams
#endif
