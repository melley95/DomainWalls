/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(INITIALSCALARDATA_2D_HPP_)
#error "This file should only be included through InitialScalarData_2D.hpp"
#endif

#ifndef INITIALSCALARDATA_2D_IMPL_HPP_
#define INITIALSCALARDATA_2D_IMPL_HPP_




inline InitialScalarData_2D::InitialScalarData_2D(params_t a_init_SF_params, double a_dx)
        : m_init_SF_params(a_init_SF_params), m_dx(a_dx)
    {
    }

// Compute the value of the initial vars on the grid
template <class data_t>
void InitialScalarData_2D::compute(Cell<data_t> current_cell) const
{
    CCZ4CartoonVars::VarsWithGauge<data_t> vars;
    VarsTools::assign(vars,
                      0.); // Set only the non-zero components explicitly below

    // Get coords and radius
    Coordinates<data_t> coords(current_cell, m_dx, m_init_SF_params.centerSF);
    data_t x = coords.x;
    data_t y = coords.y;
    

  
    data_t r2 = simd_max(x*x + y*y, 1e-12);

    data_t cos2phi = y*y/r2;
    data_t sin2phi = x*x/r2;
 

    data_t R;

    R = m_init_SF_params.R0/sqrt(cos2phi+pow(m_init_SF_params.eps1, -2)*sin2phi);

   // Real a = 35.0;
    // Real b = 70.0;

     //R = (a*b)/(a*x/sqrt(x*x+y*y)+b*y/sqrt(x*x+y*y));

    // data_t phi = tanh((sqrt(x*x+y*y)-R)/sqrt(2.0));

    data_t phi = m_init_SF_params.eta*tanh(sqrt(2.0*m_init_SF_params.lambda)*m_init_SF_params.eta*(sqrt(r2)-R)/2.0);

    data_t Pi = 0;

    vars.phi = phi;
    vars.Pi = Pi;
    vars.lapse = 1.;
    vars.chi = 1.;
    FOR(i) { vars.h[i][i] = 1.; }
    vars.hww = 1.;
    // Store the initial values of the variables
    current_cell.store_vars(vars);
}

#endif /* INITIALSCALARDATA_2D_IMPL_HPP_ */
