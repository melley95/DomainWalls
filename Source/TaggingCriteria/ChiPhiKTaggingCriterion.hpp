/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef CHIPHIKTAGGINGCRITERION_HPP_
#define CHIPHIKTAGGINGCRITERION_HPP_

#include "Cell.hpp"
 // #include "Coordinates.hpp"
#include "DimensionDefinitions.hpp"
#include "FourthOrderDerivatives.hpp"
// #include "ScalarField.hpp"
#include "Tensor.hpp"

class ChiPhiKTaggingCriterion
{
  protected:
    const double m_dx;
    const FourthOrderDerivatives m_deriv;
    const double m_threshold_chi;
    const double m_threshold_phi;
    const double m_threshold_K;

   // template <class data_t>
   // using MatterVars = typename ScalarField<>::template Vars<data_t>;

    /// Vars object for chi
  /*  template <class data_t> struct Vars
    {
        data_t chi; //!< Conformal factor

        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function)
        {
            using namespace VarsTools; // define_enum_mapping is part of
                                       // VarsTools
            define_enum_mapping(mapping_function, c_chi, chi);
        }
    };
*/
  public:
    ChiPhiKTaggingCriterion(const double dx, const double threshold_chi,
                              const double threshold_phi, const double threshold_K)
        : m_dx(dx), m_deriv(dx), m_threshold_chi(threshold_chi),
          m_threshold_phi(threshold_phi), m_threshold_K(threshold_K){};

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
      Tensor<1, data_t> d1_phi;
      FOR(idir) m_deriv.diff1(d1_phi, current_cell, idir, c_phi);

      Tensor<1, data_t> d1_chi;
      FOR(idir) m_deriv.diff1(d1_chi, current_cell, idir, c_chi);

      Tensor<1, data_t> d1_K;
      FOR(idir) m_deriv.diff1(d1_K, current_cell, idir, c_K);

        data_t mod_d1_chi = 0;
        data_t mod_d1_phi = 0;
        data_t mod_d1_K = 0;

        FOR(idir)
        {
            mod_d1_phi += d1_phi[idir] * d1_phi[idir];
            mod_d1_chi += d1_chi[idir] * d1_chi[idir];
            mod_d1_K += d1_K[idir] * d1_K[idir];
        }

        data_t criterion = m_dx * (sqrt(mod_d1_phi) / m_threshold_phi +
                                   sqrt(mod_d1_chi) / m_threshold_chi +
                                   sqrt(mod_d1_K) / m_threshold_K);

        // Write back into the flattened Chombo box
        current_cell.store_vars(criterion, 0);

        
    }
};

#endif /* CHIPHIKTAGGINGCRITERION_HPP_ */
