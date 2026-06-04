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
    


    data_t a = m_init_SF_params.a;
    data_t b = m_init_SF_params.b;


    

    data_t sign = 1;
    if (x_p*x_p/(a*a) + y_p*y_p/(b*b) < 1 ){
        sign = -1;
    }


    data_t f = x_p*x_p/(a*a) + y_p*y_p/(b*b) - 1.0;
    data_t r0 = a*a/(b*b);
    data_t sbar = get_root(r0, x_p/a, y_p/b, f);
    data_t x_1 = r0*x_p/(sbar + r0);
    data_t y_1 = y_p/(sbar + 1.0);
    
    data_t min_dist = sqrt((x_1-x_p)*(x_1-x_p) + (y_1-y_p)*(y_1-y_p)); 

    

    data_t phi = m_init_SF_params.eta*tanh(sqrt(0.5*m_init_SF_params.lambda)*m_init_SF_params.eta*(sign*min_dist));

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

template <class data_t> 
data_t Spheroid::get_root(data_t r0, data_t z0, data_t z1, data_t g) const
{
    data_t n0 = r0*z0;
    data_t s0, s1;
    data_t eps = 1e-12;

    if (g < 0) // inside ellipse
    {
        s0 = -r0 + eps;
        s1 = 0.0;
    }
    else // outside ellipse
    {
        s0 = 0.0;
        s1 = 1.0;

        // grow upper bound until root is bracketed
        while (true)
        {
            data_t ratio0 = n0/(s1 + r0);
            data_t ratio1 = z1/(s1 + 1.0);
            data_t g1 = ratio0*ratio0 + ratio1*ratio1 - 1.0;

            if (g1 < 0.0) break;
            s1 *= 2.0;
        }
    }

    data_t s = 0;
    int max_iterations = 1000;

    for (int i = 0; i < max_iterations; i++)
    {
        s = 0.5*(s0 + s1);
        if (s == s0 || s == s1) break;

        data_t ratio0 = n0/(s + r0);
        data_t ratio1 = z1/(s + 1.0);
        data_t gs = ratio0*ratio0 + ratio1*ratio1 - 1.0;

        if (gs > 0.0)
            s0 = s;
        else if (gs < 0.0)
            s1 = s;
        else
            break;
    }

    return s;
}

#endif /* SPHEROID_IMPL_HPP_ */
