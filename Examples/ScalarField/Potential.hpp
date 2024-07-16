/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SCALARPOTENTIAL_HPP_
#define SCALARPOTENTIAL_HPP_

#include "simd.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total no. components
#include "VarsTools.hpp"
#include "simd.hpp"

class Potential
{
  

  public:
     struct params_t
    {
        double eta;
        double lambda;
    };

    //! The constructor
    Potential(const params_t a_params)
        : m_params(a_params)
    {
    }

    //! Set the potential function for the scalar field here
    template <class data_t, template <typename> class vars_t>
    void compute_potential(data_t &V_of_phi, data_t &dVdphi,
                           const vars_t<data_t> &vars) const
    {
        
        V_of_phi = m_params.lambda*0.25*pow(pow(vars.phi, 2.0)-pow(m_params.eta, 2.0), 2.0);

      
        dVdphi = m_params.lambda*vars.phi*(pow(vars.phi, 2.0)-pow(m_params.eta, 2.0));
    }

 protected:
    const params_t m_params;
};

#endif /* POTENTIAL_HPP_ */
