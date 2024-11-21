/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef POTENTIAL_HPP_
#define POTENTIAL_HPP_

#include "simd.hpp"

class Potential
{
  public:
    struct params_t
    {
        // double scalar_mass;
       double eta;
       double lambda;
    };

  private:
    params_t m_params;

  public:
    //! The constructor
    Potential(params_t a_params) : m_params(a_params) {}

    //! Set the potential function for the scalar field here
    template <class data_t, template <typename> class vars_t>
    void compute_potential(data_t &V_of_phi, data_t &dVdphi,
                           const vars_t<data_t> &vars) const
    {
         V_of_phi = m_params.lambda*0.25*pow(pow(vars.phi, 2.0)-pow(m_params.eta, 2.0), 2.0);

      
         dVdphi = m_params.lambda*vars.phi*(pow(vars.phi, 2.0)-pow(m_params.eta, 2.0));
    }
};

#endif /* POTENTIAL_HPP_ */
