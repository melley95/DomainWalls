/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(ELLIPSOID_HPP_)
#error "This file should only be included through Ellipsoid.hpp"
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


    

    data_t sign = 1;
    if (x_p*x_p/(a*a) + y_p*y_p/(b*b)  + z_p*z_p/(c*c) < 1 ){
        sign = -1;
    }


    data_t f = x_p*x_p/(a*a) + y_p*y_p/(b*b)  + z_p*z_p/(c*c) - 1.0;
    data_t r0 = a*a/(c*c);
    data_t r1 = b*b/(c*c);
    data_t sbar = get_root(r0, r1, x_p/a, y_p/b, z_p/c, f);
    data_t x_1 = r0*x_p/(sbar + r0);
    data_t y_1 = r1*y_p/(sbar + r1);
    data_t z_1 = z_p/(sbar + 1.0);
    
    data_t min_dist = sqrt((x_1-x_p)*(x_1-x_p) + (y_1-y_p)*(y_1-y_p) + (z_1-z_p)*(z_1-z_p)); 

    

    data_t phi = m_init_SF_params.eta*tanh(sqrt(0.5*m_init_SF_params.lambda)*m_init_SF_params.eta*(sign*min_dist));

    data_t Pi = 0;

    vars.phi = phi;
    vars.Pi = Pi;
    vars.lapse = 1.;
    vars.chi = 1.;
    FOR(i) { vars.h[i][i] = 1.; }

    // Store the initial values of the variables
    current_cell.store_vars(vars);
}

template <class data_t> 
data_t Ellipsoid::get_root(data_t r0, data_t r1, data_t z0, data_t z1,  data_t z2, data_t g) const
{
    data_t n0 = r0*z0;
    data_t n1 = r1*z1;
    data_t s0 = z2 - 1.0;
    data_t s1 = -1 + sqrt(r0*r0*z0*z0 + r1*r1*z1*z1 + z2*z2);
    data_t s = 0;
    data_t max_iterations = 1000; // std::numeric_limits<Real>::digits - std::numeric_limits<Real>::min exponent; // ME: change??

    for (int i =0; i < max_iterations; i++){
        s = (s0 + s1) / 2.0;
        if (s == s0 || s == s1) {
            break;
        }
        data_t ratio0 = n0/(s + r0);
        data_t ratio1 = n1/(s+r1);
        data_t ratio2 = z2/(s+1);
        g = ratio0*ratio0 + ratio1*ratio1 + ratio2*ratio2  -1.0;
        if (g > 0){
            s0 = s;
        }
        else if(g < 0)
        {
            s1 = s;
        }
        else{
            break;
        }

    }
    return s;
}

#endif /* ELLIPSOID_IMPL_HPP_ */
