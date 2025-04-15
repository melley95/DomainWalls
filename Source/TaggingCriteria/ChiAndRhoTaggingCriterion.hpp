/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef CHIANDRHOTAGGINGCRITERION_HPP_
#define CHIANDRHOTAGGINGCRITERION_HPP_

#include "Cell.hpp"
 // #include "Coordinates.hpp"
#include "DimensionDefinitions.hpp"
#include "FourthOrderDerivatives.hpp"
// #include "ScalarField.hpp"
#include "Tensor.hpp"

class ChiAndRhoTaggingCriterion
{
  protected:
    const double m_dx;
    const FourthOrderDerivatives m_deriv;
    const double m_threshold_chi;
    const double m_threshold_rho;

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
    ChiAndRhoTaggingCriterion(const double dx, const double threshold_chi,
                              const double threshold_rho)
        : m_dx(dx), m_deriv(dx), m_threshold_chi(threshold_chi),
          m_threshold_rho(threshold_rho){};

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
      data_t rho = current_cell.load_vars(c_rho);

      Tensor<1, data_t> d1_chi;
      FOR(idir) m_deriv.diff1(d1_chi, current_cell, idir, c_chi);

        data_t mod_d1_chi = 0;
     
        FOR(idir)
        {
          
            mod_d1_chi += d1_chi[idir] * d1_chi[idir];
        }

        data_t criterion = m_dx * m_dx * rho / m_threshold_rho +  m_dx *  sqrt(mod_d1_chi) / m_threshold_chi;

                               \
        // Write back into the flattened Chombo box
        current_cell.store_vars(criterion, 0);

        
    }
};

#endif /* CHIANDRHOTAGGINGCRITERION_HPP_ */
