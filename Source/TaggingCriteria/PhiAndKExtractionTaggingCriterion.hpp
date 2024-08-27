/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef PHIANDKEXTRACTIONTAGGINGCRITERION_HPP_
#define PHIANDKEXTRACTIONTAGGINGCRITERION_HPP_

#include "Cell.hpp"
#include "DimensionDefinitions.hpp"
#include "FourthOrderDerivatives.hpp"
#include "SphericalExtraction.hpp"
#include "Tensor.hpp"

class PhiAndKExtractionTaggingCriterion
{
  protected:
    const double m_dx;
    const FourthOrderDerivatives m_deriv;
    const double m_threshold_phi;
    const double m_threshold_K;

    const spherical_extraction_params_t m_params;
    const int m_level;
    const bool m_activate_extraction;


  public:
    PhiAndKExtractionTaggingCriterion(double dx, double threshold_phi, double threshold_K, const int a_level, 
                            const spherical_extraction_params_t a_params,
                            const bool activate_extraction = false)
     : m_dx(dx), m_deriv(dx), m_threshold_phi(threshold_phi),
       m_threshold_K(threshold_K), m_level(a_level), m_params(a_params), 
       m_activate_extraction(activate_extraction){};

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        Tensor<1, data_t> d1_phi;
        FOR(idir) m_deriv.diff1(d1_phi, current_cell, idir, c_phi);

        Tensor<1, data_t> d1_K;
        FOR(idir) m_deriv.diff1(d1_K, current_cell, idir, c_K);

        data_t mod_d1_phi = 0;
        data_t mod_d1_K = 0;
        FOR(idir)
        {
            mod_d1_phi += d1_phi[idir] * d1_phi[idir];
            mod_d1_K += d1_K[idir] * d1_K[idir];
        }

        data_t criterion = m_dx * (sqrt(mod_d1_phi) / m_threshold_phi +
                                   sqrt(mod_d1_K) / m_threshold_K);


        if (m_activate_extraction)
        {
            for (int iradius = 0; iradius < m_params.num_extraction_radii;
                 ++iradius)
            {
                // regrid if within extraction level and not at required
                // refinement
                if (m_level < m_params.extraction_levels[iradius])
                {
                    const Coordinates<data_t> coords(current_cell, m_dx,
                                                     m_params.center);
                    const data_t r = coords.get_radius();
                    // add a 20% buffer to extraction zone so not too near to
                    // boundary
                    auto regrid = simd_compare_lt(
                        r, 1.2 * m_params.extraction_radii[iradius]);
                    criterion = simd_conditional(regrid, 100.0, criterion);
                }
            }
        }
        // Write back into the flattened Chombo box
        current_cell.store_vars(criterion, 0);
    }
};

#endif /* PHIANDKEXTRACTIONTAGGINGCRITERION_HPP_ */
