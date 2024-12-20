#ifndef MATTERPARAMS_HPP_
#define MATTERPARAMS_HPP_

#include "REAL.H"

namespace MatterParams
{

struct params_t
{
    Real phi_0;
    Real dphi;
    Real dphi_wavelength;
    Real pi_0;
    Real dpi;
    Real dpi_wavelength;
    Real scalar_mass;
};

inline void read_params(GRParmParse &pp, params_t &matter_params)
{
    pp.get("phi_0", matter_params.phi_0);
    pp.get("dphi", matter_params.dphi);
    pp.get("dphi_wavelength", matter_params.dphi_wavelength);
    pp.get("pi_0", matter_params.pi_0);
    pp.get("dpi", matter_params.dpi);
    pp.get("dpi_wavelength", matter_params.dpi_wavelength);
    pp.get("scalar_mass", matter_params.scalar_mass);
}

}; // namespace MatterParams

#endif