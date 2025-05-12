/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(SPHEROID_HPP_)
#error "This file should only be included through Spheroid.hpp"
#endif

#ifndef SPHEROID_IMPL_HPP_
#define SPHEROID_IMPL_HPP_




inline Spheroid::Spheroid(params_t a_init_SF_params, double a_dx)
        : m_init_SF_params(a_init_SF_params), m_dx(a_dx)
    {
    }

// Compute the value of the initial vars on the grid
template <class data_t>
void Spheroid::compute(Cell<data_t> current_cell) const
{
    CCZ4CartoonVars::VarsWithGauge<data_t> vars;
    VarsTools::assign(vars,
                      0.); // Set only the non-zero components explicitly below

    // Get coords and radius
    Coordinates<data_t> coords(current_cell, m_dx, m_init_SF_params.centerSF);
    data_t x_p = coords.x;
    data_t y_p = coords.y;
    

  
   // data_t r2 = simd_max(x*x + y*y, 1e-12);

   // data_t cos2phi = y*y/r2;
  //  data_t sin2phi = x*x/r2;
 

    data_t R;

    //R = m_init_SF_params.R0/sqrt(cos2phi+pow(m_init_SF_params.eps1, -2)*sin2phi);

    data_t a = m_init_SF_params.a;
    data_t b = m_init_SF_params.b;


    

    data_t min_theta = 0.0;
    data_t min_dist = std::numeric_limits<double>::max();
    data_t best_theta = 0.0;

    int steps = 10000;
    for (int i = 0; i <= steps; ++i) {
        double theta = 2.0 * M_PI * i / steps;
        data_t x = a * std::cos(theta);
        data_t y = b * std::sin(theta);
        data_t dist= (x - x_p) * (x - x_p) + (y - y_p) * (y - y_p);
       
        if (dist < min_dist) {
            min_dist = dist;
            best_theta = theta;
        }
    }
    data_t sign = 1;
    if (x_p*x_p/(a*a) + y_p*y_p/(b*b) < 1 ){
        sign = -1;
    }

    

    data_t phi = m_init_SF_params.eta*tanh(sqrt(0.5*m_init_SF_params.lambda)*m_init_SF_params.eta*(sign*sqrt(min_dist)));

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

#endif /* SPHEROID_IMPL_HPP_ */
