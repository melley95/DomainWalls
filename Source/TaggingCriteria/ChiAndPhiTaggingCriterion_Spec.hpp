/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef CHIANDPHITAGGINGCRITERION_SPEC_HPP_
#define CHIANDPHITAGGINGCRITERION_SPEC_HPP_

#include "Cell.hpp"
 // #include "Coordinates.hpp"
#include "DimensionDefinitions.hpp"
#include "FourthOrderDerivatives.hpp"
// #include "ScalarField.hpp"
#include "Tensor.hpp"

class ChiAndPhiTaggingCriterion_Spec
{
  protected:
    const double m_dx;
    const FourthOrderDerivatives m_deriv;
    const double m_threshold_chi;
    const double m_threshold_phi;

    const double m_origin_x;
    const double m_origin_y;

    const double m_L;

    const int m_level;
    const int m_max_level;
  

 

    std::array<double, CH_SPACEDIM> m_center;
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
    ChiAndPhiTaggingCriterion_Spec(const double dx, const double threshold_chi,
                              const double threshold_phi, const double origin_x,
                            const double origin_y, const double L, std::array<double, CH_SPACEDIM> center, const int level,
                          const int max_level)
        : m_dx(dx), m_deriv(dx), m_threshold_chi(threshold_chi),
          m_threshold_phi(threshold_phi), m_origin_x(origin_x), m_origin_y(origin_y), m_L(L),
          m_center(center), m_level(level), m_max_level(max_level) {};

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
      Tensor<1, data_t> d1_phi;
      FOR(idir) m_deriv.diff1(d1_phi, current_cell, idir, c_phi);

      Tensor<1, data_t> d1_chi;
      FOR(idir) m_deriv.diff1(d1_chi, current_cell, idir, c_chi);

        data_t mod_d1_chi = 0;
        data_t mod_d1_phi = 0;

        FOR(idir)
        {
            mod_d1_phi += d1_phi[idir] * d1_phi[idir];
            mod_d1_chi += d1_chi[idir] * d1_chi[idir];
        }

        data_t criterion = m_dx * (sqrt(mod_d1_phi) / m_threshold_phi +
                                   sqrt(mod_d1_chi) / m_threshold_chi);



        const Coordinates<double> coords(current_cell, m_dx, m_center);

        

        data_t x = coords.x;
        double y = coords.y;

        data_t L_x= pow((x - m_origin_x),2.0);
        double L_y = pow((y - m_origin_y),2.0);

        data_t dist = sqrt(L_x + L_y);

        auto regrid = simd_compare_lt(
          dist, m_L);
        
        if (m_level <= m_max_level){
        
        criterion = simd_conditional(regrid, 100.0, criterion);

        }

        else if (m_level > m_max_level){
        criterion = simd_conditional(regrid, 0.0, criterion);

        }
        


        // Write back into the flattened Chombo box
        current_cell.store_vars(criterion, 0);

        
    }
};

#endif /* CHIANDPHITAGGINGCRITERION_SPEC_HPP_ */
