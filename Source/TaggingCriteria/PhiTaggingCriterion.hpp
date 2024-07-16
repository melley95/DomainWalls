/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef PHITAGGINGCRITERION_HPP_
#define PHITAGGINGCRITERION_HPP_

#include "Cell.hpp"
#include "DimensionDefinitions.hpp"

#include "Tensor.hpp"

class PhiTaggingCriterion
{
//  protected:
 //   const double m_dx;
   // const FourthOrderDerivatives m_deriv;

  public:
   // PhiTaggingCriterion();

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        auto phi = current_cell.load_vars(c_phi);
        //Tensor<1, data_t> d1_chi;
        //FOR(idir) m_deriv.diff1(d1_chi, current_cell, idir, c_chi);

        //data_t mod_d1_chi = 0;
        //FOR(idir) mod_d1_chi += d1_chi[idir] * d1_chi[idir];
        data_t criterion = -phi;

        // Write back into the flattened Chombo box
        current_cell.store_vars(criterion, 0);
    }
};

#endif /* PHITAGGINGCRITERION_HPP_ */
