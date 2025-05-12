/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef CHIANDPHITAGGINGCRITERION_HPP_
#define CHIANDPHITAGGINGCRITERION_HPP_

#include "Cell.hpp"
 // #include "Coordinates.hpp"
#include "DimensionDefinitions.hpp"
#include "FourthOrderDerivatives.hpp"
// #include "ScalarField.hpp"
#include "Tensor.hpp"

class ChiAndPhiTaggingCriterion
{
  protected:
    const double m_dx;
    const double m_dt;
    const FourthOrderDerivatives m_deriv;
    const double m_threshold_chi;
    const double m_threshold_phi;

    const double m_time;
    const double m_rebound_time;
    const double m_threshold_chi_rebound;
    const double m_threshold_phi_rebound;
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
    ChiAndPhiTaggingCriterion(const double dx, const double dt , const double threshold_chi,
                              const double threshold_phi, const double time, const double rebound_time, const double threshold_chi_rebound, const double threshold_phi_rebound)
        : m_dx(dx), m_dt(dt), m_deriv(dx), m_threshold_chi(threshold_chi),
          m_threshold_phi(threshold_phi), m_time(time), m_rebound_time{rebound_time}, m_threshold_chi_rebound(threshold_chi_rebound), m_threshold_phi_rebound(threshold_phi_rebound){};

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
      Tensor<1, data_t> d1_phi;
      FOR(idir) m_deriv.diff1(d1_phi, current_cell, idir, c_phi);

      Tensor<1, data_t> d1_chi;
      FOR(idir) m_deriv.diff1(d1_chi, current_cell, idir, c_chi);

        data_t pi = current_cell.load_vars(c_Pi);

        data_t mod_d1_chi = 0;
        data_t mod_d1_phi = 0;
   


        FOR(idir)
        {
            mod_d1_phi += d1_phi[idir] * d1_phi[idir];
            
            mod_d1_chi += d1_chi[idir] * d1_chi[idir];
        }
        data_t criterion = 0.0;
        if (m_time < m_rebound_time){


        criterion = m_dx * (sqrt(mod_d1_phi) / m_threshold_phi) + m_dt * (sqrt(pi*pi) / m_threshold_phi)
                                + m_dx  * (sqrt(mod_d1_chi) / m_threshold_chi);
        }

        else{

        criterion = m_dx * (sqrt(mod_d1_phi) / m_threshold_phi_rebound) + m_dt * (sqrt(pi*pi) / m_threshold_phi_rebound)
          + m_dx  * (sqrt(mod_d1_chi) / m_threshold_chi_rebound);  

        }

        // Write back into the flattened Chombo box
        current_cell.store_vars(criterion, 0);

        
    }
};

#endif /* CHIANDPHITAGGINGCRITERION_HPP_ */
