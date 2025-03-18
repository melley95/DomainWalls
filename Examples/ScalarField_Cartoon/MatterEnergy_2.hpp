/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

 #ifndef MATTERENERGY_2_HPP_
 #define MATTERENERGY_2_HPP_
 
 #include "CCZ4CartoonGeometry.hpp"
 #include "CCZ4CartoonVars.hpp"
 #include "Cell.hpp"
 #include "CoordinateTransformations.hpp"
 #include "Coordinates.hpp"
 #include "FourthOrderDerivatives.hpp"
 #include "GRInterval.hpp"
 #include "CCZ4Cartoon.hpp"
 #include "Tensor.hpp"
 #include "UserVariables.hpp" //This files needs NUM_VARS - total number of components
 #include "VarsTools.hpp"
 #include "simd.hpp"
 #include "CCZ4Vars.hpp"
 #include "MovingPunctureGauge.hpp"
 
 //! Calculates the MatterEnergy rho with type matter_t and writes it to the grid
 template <class potential_t, class gauge_t = MovingPunctureGauge> class MatterEnergy_2
 {
     // Use the variable definition in CCZ4
     template <class data_t>
     using Vars = typename CCZ4Cartoon<potential_t>::template Vars<data_t>;
     template <class data_t>
    using Diff2Vars = CCZ4CartoonVars::Diff2VarsNoGauge<data_t>;
 
   protected:
     const FourthOrderDerivatives
         m_deriv; //!< An object for calculating derivatives of the variables
     const std::array<double, CH_SPACEDIM> m_center;
     const double m_dx;
     potential_t m_potential;
     gauge_t m_gauge;
 
   public:

     MatterEnergy_2(gauge_t a_gauge, potential_t a_potential, double a_dx,
                  std::array<double, CH_SPACEDIM> a_center)
         : m_gauge(a_gauge), m_potential(a_potential), m_dx(a_dx), m_deriv(a_dx), m_center(a_center)
     {
     }
 
     template <class data_t> void compute(Cell<data_t> current_cell) const
     {
        // copy data from chombo gridpoint into local variables, derivs
        const auto vars = current_cell.template load_vars<Vars>();
        const auto d1 = m_deriv.template diff1<Vars>(current_cell);
        const auto d2 = m_deriv.template diff2<Diff2Vars>(current_cell);
        const auto advec =
            m_deriv.template advection<Vars>(current_cell, vars.shift);

        Vars<data_t> rhs;

        m_gauge.rhs_gauge(rhs, vars, d1, d2, advec);
         // set up the metric quantities needed
         using namespace TensorAlgebra;
 
 
         auto h_UU = TensorAlgebra::compute_inverse_sym(vars.h);
         auto h_UU_ww = 1. / vars.hww;
         auto gamma_UU_ww = vars.chi * h_UU_ww ; // CHECK
         auto chris = TensorAlgebra::compute_christoffel(d1.h, h_UU);
         const int dI = CH_SPACEDIM - 1;
         const int nS = GR_SPACEDIM - CH_SPACEDIM; //!< Dimensions of the transverse sphere
 
         Coordinates<data_t> coords{current_cell, m_dx /*, m_center*/};
         const double cartoon_coord = coords.y;
        
        
 
             // Matter contributions
         data_t V_of_phi = 0.0;
         data_t dVdphi = 0.0;
         m_potential.compute_potential(V_of_phi, dVdphi, vars);
 
         emtensorCartoon_t<data_t> emtensor = compute_SF_EM_tensor(vars, d1, h_UU, h_UU_ww, chris, nS, V_of_phi); 
             
         const data_t det_gamma = pow(vars.chi, -1.5) * 2.0 * M_PI * abs(coords.y);
         Tensor<2, data_t> vars_gamma, vars_K_tensor;
         const data_t vars_K_ww = 1.0 / vars.chi * (vars.Aww + 1.0 / 3.0 * vars.hww * vars.K);
         FOR2(i, j)
         {
             vars_gamma[i][j] = vars.h[i][j] / vars.chi;
             vars_K_tensor[i][j] =
                 1.0 / vars.chi *
                 (vars.A[i][j] + 1.0 / 3.0 * vars.h[i][j] * vars.K);  //THIS RIGHT??
         }
         const auto gamma_UU = compute_inverse_sym(vars_gamma);
         const Tensor<3, data_t> chris_phys =
             compute_phys_chris(d1.chi, vars.chi, vars.h, h_UU, chris.ULL);
         // coordinate quantities
         // The unit covector normal to the spherical surface
  /*       Coordinates<data_t> coords(current_cell, m_dx, m_center);
         Tensor<1, data_t> si_L, Ni_L;
         data_t R = coords.get_radius();
         si_L[0] = coords.x / R;
         si_L[1] = coords.y / R;
         si_L[2] = coords.z / R;
         // normalise this to 1 using spatial metric
         data_t mod_N2 = 0.0;
         FOR2(i, j) { mod_N2 += gamma_UU[i][j] * si_L[i] * si_L[j]; }
         FOR1(i) { Ni_L[i] = si_L[i] / sqrt(mod_N2); }
 
         // the area element of the sphere - check matches up
         data_t rxy2 =
             simd_max(coords.x * coords.x + coords.y * coords.y, 1e-12);
         data_t r2sintheta = sqrt(rxy2) * R;
         using namespace CoordinateTransformations;
         Tensor<2, data_t> spherical_gamma =
             cartesian_to_spherical_LL(vars_gamma, coords.x, coords.y, coords.z);
         data_t sqrt_det_Sigma = area_element_sphere(spherical_gamma);
 */
         // calculate according to Landau Lifshitz method
         data_t rho1 = emtensor.rho * det_gamma;
 
   /*      // Energy flux
         data_t flux1 = 0.0;
         FOR1(i)
         {
             flux1 += -si_L[i] * vars.shift[i] * emtensor.rho;
 
             FOR1(j)
             {
                 flux1 += vars.lapse * si_L[i] * emtensor.Si[j] * gamma_UU[i][j];
             }
         }
         flux1 *= det_gamma;
 */
         // calculate the E source
         data_t source1 = -emtensor.rho * rhs.lapse - vars.lapse * vars.lapse * vars_K_ww * emtensor.Sww * gamma_UU_ww * gamma_UU_ww
                          +   vars.lapse * emtensor.Sww * gamma_UU_ww * vars.shift[1]/cartoon_coord ;
                       
         FOR1(i)
         {
            source1 += emtensor.Si[i]*rhs.shift[i];
           
             FOR1(j)
             {
         
                 FOR1(k)
                 {
                    source1 += 0.5 * vars.lapse *  emtensor.Sij[i][j] *
                    (gamma_UU[i][k] * d1.shift[j][k] + gamma_UU[j][k] * d1.shift[i][k]);
                 
 
                     FOR1(l)
                     {
                         source1 -= vars.lapse * vars.lapse * emtensor.Sij[i][j] *
                                    gamma_UU[i][k] * gamma_UU[j][l] *
                                    vars_K_tensor[k][l];
                         source1 += 0.5 * vars.lapse *  emtensor.Sij[i][j] * (gamma_UU[i][k] * vars.shift[l]*chris_phys[j][l][k] 
                                    + gamma_UU[j][k] * vars.shift[l]*chris_phys[i][l][k]);
                     }
                 }
             }
         }
         source1 *= det_gamma;
 
         // assign values of MatterEnergy in output box
         current_cell.store_vars(rho1, c_rhoLL);
        // current_cell.store_vars(flux1, c_flux);
         current_cell.store_vars(source1, c_source);
     }
 };
 
 #endif /* MATTERENERGY_2_HPP_ */