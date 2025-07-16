/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef MOVINGPUNCTUREGAUGECOSMO_HPP_
#define MOVINGPUNCTUREGAUGECOSMO_HPP_

#include "DimensionDefinitions.hpp"
#include "Tensor.hpp"
#include "MovingPunctureGauge.hpp"


/// This is an example of a gauge class that can be used in the CCZ4RHS compute
/// class
/**
 * This class implements the shock avoiding condition of https://arxiv.org/pdf/2207.06376.pdf
 **/

class MovingPunctureGaugeCosmo
{
  public:
    using params_t = MovingPunctureGauge::params_t; 
   // matter_t my_matter;
    params_t m_params;
    

    
    MovingPunctureGaugeCosmo(const MovingPunctureGauge::params_t &a_params) : m_params(a_params) {}

    template <class data_t, template <typename> class vars_t,
              template <typename> class diff2_vars_t>
    inline void rhs_gauge(vars_t<data_t> &rhs, const vars_t<data_t> &vars,
                          const vars_t<Tensor<1, data_t>> &d1,
                          const diff2_vars_t<Tensor<2, data_t>> &d2,
                          const vars_t<data_t> &advec) const
    {


     /* const int nS =
      GR_SPACEDIM - CH_SPACEDIM; //!< Dimensions of the transverse sphere

      auto h_UU = TensorAlgebra::compute_inverse_sym(vars.h);
      auto h_UU_ww = 1. / vars.hww;
      auto chris = TensorAlgebra::compute_christoffel(d1.h, h_UU);

      Tensor<1, data_t> chris_ww;
      FOR(i)
      {
        chris_ww[i] =
            one_over_cartoon_coord * (delta(i, dI) - h_UU[i][dI] * vars.hww);
        FOR(j) chris_ww[i] -= 0.5 * h_UU[i][j] * d1.hww[j];
       }


       auto A_UU = TensorAlgebra::raise_all(vars.A, h_UU);
       data_t tr_A2 = TensorAlgebra::compute_trace(vars.A, A_UU) +
                      nS * vars.Aww * vars.Aww * h_UU_ww * h_UU_ww;

          // Matter contributions
      data_t V_of_phi = 0.0;
      data_t dVdphi = 0.0;
      m_potential.compute_potential(V_of_phi, dVdphi, vars);
    
       emtensorCartoon_t<data_t> emtensor = compute_SF_EM_tensor(vars, d1, h_UU, h_UU_ww, chris, nS, V_of_phi);
       */
     


        rhs.lapse = m_params.lapse_advec_coeff * advec.lapse -
                        (1.0/3.0)*(vars.lapse * vars.lapse + m_params.lapse_coeff) *
                        (vars.K  + sqrt(24.0*M_PI*emtensor.rho + (3.0/2.0)*tr_A2)- 2 * vars.Theta);
        FOR(i)
        {
            rhs.shift[i] = m_params.shift_advec_coeff * advec.shift[i] +
                           m_params.shift_Gamma_coeff * vars.B[i];
            rhs.B[i] = m_params.shift_advec_coeff * advec.B[i] -
                       m_params.shift_advec_coeff * advec.Gamma[i] +
                       (3.0/4.0)*rhs.Gamma[i]  - m_params.eta * vars.B[i];
            FOR(j)
            {
              rhs.B[i] += vars.lapse*h_UU[i][j](d1.K[j] + 16.0*M_PI*emtensor.Si[j]);  // Should have G_Netwon here
            }
        }
    }



};

#endif /* MOVINGPUNCTUREGAUGECOSMO_HPP_ */