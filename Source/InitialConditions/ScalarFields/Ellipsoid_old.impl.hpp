/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(ELLIPSOID_HPP_)
#error "This file should only be included through Spheroid.hpp"
#endif

#ifndef ELLIPSOID_IMPL_HPP_
#define ELLIPSOID_IMPL_HPP_




inline Ellipsoid::Ellipsoid(params_t a_init_SF_params, double a_dx)
        : m_init_SF_params(a_init_SF_params), m_dx(a_dx)
    {
    }

// Compute the value of the initial vars on the grid
template <class data_t>
void Ellipsoid::compute(Cell<data_t> current_cell) const
{
    MatterCCZ4RHS<ScalarField<>>::Vars<data_t> vars;
    VarsTools::assign(vars,
                      0.); // Set only the non-zero components explicitly below

    // Get coords and radius
    Coordinates<data_t> coords(current_cell, m_dx, m_init_SF_params.centerSF);
    data_t x_p = coords.x;
    data_t y_p = coords.y;
    data_t z_p = coords.z;
    



    data_t a = m_init_SF_params.a;
    data_t b = m_init_SF_params.b;
    data_t c = m_init_SF_params.c;

    

    data_t min_theta = 0.0;
    data_t min_dist = std::numeric_limits<double>::max();
    data_t best_theta = 0.0;

    int steps = 1000;
    for (int i = 0; i <= steps; ++i) {
        for (int j = 0; j <= steps; ++j){
            double theta =   M_PI * i / (2.0*steps);
            double alpha =   M_PI * j / (2.0*steps);
            data_t x = a * std::cos(theta) * std::sin(alpha);
            data_t y = b * std::sin(theta) * std::sin(alpha);
            data_t z = c * std::cos(alpha);
            data_t dist= (x - x_p) * (x - x_p) + (y - y_p) * (y - y_p) + (z - z_p) * (z - z_p);
       
        if (dist < min_dist) {
            min_dist = dist;
            best_theta = theta;
        }
    }
}
    data_t sign = 1;
    if (x_p*x_p/(a*a) + y_p*y_p/(b*b) + z_p*z_p/(c*c) < 1 ){
        sign = -1;
    }

    

    data_t phi = m_init_SF_params.eta*tanh(sqrt(0.5*m_init_SF_params.lambda)*m_init_SF_params.eta*(sign*sqrt(min_dist)));

    data_t Pi = 0;

    vars.phi = phi;
    vars.Pi = Pi;
    vars.lapse = 1.;
    vars.chi = 1.;
    FOR(i) { vars.h[i][i] = 1.; }

    // Store the initial values of the variables
    current_cell.store_vars(vars);
}

#endif /* ELLIPSOID_IMPL_HPP_ */
