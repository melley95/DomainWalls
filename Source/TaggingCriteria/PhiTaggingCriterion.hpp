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
  protected:
    const double m_dx;
    const double m_time;
    const int m_level;

    std::array<double, 10> m_ref_times;
    std::array<int, 10> m_ref_levels;
    
    const FourthOrderDerivatives m_deriv;
    const double m_threshold_chi;
    const double m_threshold_K;



  public:
    PhiTaggingCriterion(const double dx, double time, int level, std::array<double, 10> ref_times, std::array<int, 10> ref_levels, const double threshold_chi, const double threshold_K)
    : m_dx(dx), m_deriv(dx), m_time(time), m_level(level), m_ref_times(ref_times),
      m_ref_levels(ref_levels), m_threshold_chi(threshold_chi), m_threshold_K(threshold_K) {};

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        data_t phi = current_cell.load_vars(c_phi);
        data_t chi = current_cell.load_vars(c_chi);

        //Tensor<1, data_t> d1_chi;
        //FOR(idir) m_deriv.diff1(d1_chi, current_cell, idir, c_chi);

        //data_t mod_d1_chi = 0;
        //FOR(idir) mod_d1_chi += d1_chi[idir] * d1_chi[idir];
        data_t criterion = 0.0;
        int size = 9;
        for (int i = 0; i < size; i++) {
            if (m_time > m_ref_times[i]){
             
                    if(m_level < m_ref_levels[i]){

                        auto crit1 = simd_compare_gt(phi, -0.012);
                        criterion = simd_conditional(crit1, 100.0, criterion);
                        auto crit2 = simd_compare_lt(phi, 0.012);
                        criterion = simd_conditional(crit2, criterion, 0.0);

                        auto crit3 = simd_compare_lt(phi, -0.017);
                        criterion = simd_conditional(crit3, 100.0, criterion);

                        auto crit4 = simd_compare_gt(phi, 0.017);
                        criterion = simd_conditional(crit4, 100.0, criterion);

                        auto crit5 = simd_compare_lt(chi, m_threshold_chi);
                        criterion = simd_conditional(crit5, 100.0, criterion);
                        
                    }
                

            }
        }
        
      Tensor<1, data_t> d1_K;
      FOR(idir) m_deriv.diff1(d1_K, current_cell, idir, c_K);
      data_t mod_d1_K = 0;

      FOR(idir)
        {
            mod_d1_K += d1_K[idir] * d1_K[idir];
        }

      criterion += m_dx * (sqrt(mod_d1_K) / m_threshold_K);

    /*  Tensor<1, data_t> d1_chi;
      FOR(idir) m_deriv.diff1(d1_chi, current_cell, idir, c_chi);

        data_t mod_d1_chi = 0;
        

        FOR(idir)
        {
            mod_d1_K += d1_K[idir] * d1_K[idir];
            mod_d1_chi += d1_chi[idir] * d1_chi[idir];
        }

       
        */

        
    
       

        
      

    //    if (simd_compare_gt(phi, -0.015) && simd_compare_lt(phi, -0.015)){
     //   criterion = 100.0; 
     //   }
        // Write back into the flattened Chombo box
        current_cell.store_vars(criterion, 0);
    }
};

#endif /* PHITAGGINGCRITERION_HPP_ */
