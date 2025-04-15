/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef RHOTAGGINGCRITERION_HPP_
#define RHOTAGGINGCRITERION_HPP_

#include "Cell.hpp"
#include "DimensionDefinitions.hpp"

#include "Tensor.hpp"

class RhoTaggingCriterion
{
  protected:
    const double m_time;
    const int m_level;
    std::array<double, 10> m_ref_times;
    std::array<int, 10> m_ref_levels;
    std::array<double, 10> m_ref_rho;

  public:
    RhoTaggingCriterion(double time, int level, std::array<double, 10> ref_times, std::array<double, 10> ref_rho, std::array<int, 10> ref_levels)
    : m_time(time),  m_level(level), m_ref_times(ref_times), m_ref_rho(ref_rho),
      m_ref_levels(ref_levels){};

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        data_t rho = current_cell.load_vars(c_rho);
        //Tensor<1, data_t> d1_chi;
        //FOR(idir) m_deriv.diff1(d1_chi, current_cell, idir, c_chi);

        //data_t mod_d1_chi = 0;
        //FOR(idir) mod_d1_chi += d1_chi[idir] * d1_chi[idir];
        data_t criterion = 0.0;
        int size = 9;
        for (int i = 0; i < size; i++) {
            if (m_time > m_ref_times[i]){

                    
                    if(m_level == m_ref_levels[i]){
                      auto crit =  simd_compare_gt(rho, m_ref_rho[i]);

                      criterion = simd_conditional(crit, 100.0, criterion);
              

                   
                    
                        
                    
                

            }
        }
        
      }

      

    //    if (simd_compare_gt(phi, -0.015) && simd_compare_lt(phi, -0.015)){
     //   criterion = 100.0; 
     //   }
        // Write back into the flattened Chombo box
        current_cell.store_vars(criterion, 0);
    }
};

#endif /* RHOTAGGINGCRITERION_HPP_ */
