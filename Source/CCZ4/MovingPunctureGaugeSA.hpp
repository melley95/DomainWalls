/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef MOVINGPUNCTUREGAUGESA_HPP_
#define MOVINGPUNCTUREGAUGESA_HPP_

#include "DimensionDefinitions.hpp"
#include "Tensor.hpp"
#include "MovingPunctureGauge.hpp"

/// This is an example of a gauge class that can be used in the CCZ4RHS compute
/// class
/**
 * This class implements the shock avoiding condition of https://arxiv.org/pdf/2207.06376.pdf
 **/
class MovingPunctureGaugeSA
{
  public:
    using params_t = MovingPunctureGauge::params_t; 

    params_t m_params;

    MovingPunctureGaugeSA(const MovingPunctureGauge::params_t &a_params) : m_params(a_params) {}

    template <class data_t, template <typename> class vars_t,
              template <typename> class diff2_vars_t>
    inline void rhs_gauge(vars_t<data_t> &rhs, const vars_t<data_t> &vars,
                          const vars_t<Tensor<1, data_t>> &d1,
                          const diff2_vars_t<Tensor<2, data_t>> &d2,
                          const vars_t<data_t> &advec) const
    {
        rhs.lapse = m_params.lapse_advec_coeff * advec.lapse -
                        (vars.lapse * vars.lapse + m_params.lapse_coeff) *
                        (vars.K - 2 * vars.Theta);
        FOR(i)
        {
            rhs.shift[i] = m_params.shift_advec_coeff * advec.shift[i] +
                           m_params.shift_Gamma_coeff * vars.B[i];
            rhs.B[i] = m_params.shift_advec_coeff * advec.B[i] -
                       m_params.shift_advec_coeff * advec.Gamma[i] +
                       rhs.Gamma[i] - m_params.eta * vars.B[i];
        }
    }
};

#endif /* MOVINGPUNCTUREGAUGESA_HPP_ */