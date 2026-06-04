#include "ScalarField.hpp"

Real ScalarField::my_potential_function(const Real &phi_here) const
{
    return m_matter_params.lambda*0.25*pow(pow(phi_here, 2.0)-pow(m_matter_params.eta, 2.0), 2.0);
}

Real ScalarField::my_phi_function(const RealVect &loc) const
{
    Real r2 = max(loc[0] * loc[0] + loc[1] * loc[1], 1e-12);

    Real x_p = loc[0];
    Real y_p = loc[1];

    Real a = m_matter_params.a;
    Real b = m_matter_params.b;
    

    Real sign = 1;
    if (x_p*x_p/(a*a) + y_p*y_p/(b*b) < 1 ){
        sign = -1;
    }

    Real f = x_p*x_p/(a*a) + y_p*y_p/(b*b) - 1.0;
    Real r0 = a*a/(b*b);

    
    Real sbar = get_root(r0, x_p/a, y_p/b, f);
    Real x_1 = r0*x_p/(sbar + r0);
    Real y_1 = y_p/(sbar + 1.0);
    
    Real min_dist = sqrt((x_1-x_p)*(x_1-x_p) + (y_1-y_p)*(y_1-y_p)); 

    

    Real phi = m_matter_params.eta*tanh(sqrt(0.5*m_matter_params.lambda)*m_matter_params.eta*(sign*min_dist));
    


    return phi;
}

Real ScalarField::my_Pi_function(const RealVect &loc) const
{
    
    return 0.0;
}


Real ScalarField::get_root(Real r0, Real z0, Real z1, Real g) const
{
    Real n0 = r0*z0;
    Real s0 = z1 - 1;
    Real s1 = -1 + sqrt(r0*r0*z0*z0 + z1*z1);
    Real s = 0;
    Real max_iterations = 1000; // std::numeric_limits<Real>::digits - std::numeric_limits<Real>::min exponent; // ME: change??

    for (int i =0; i < max_iterations; i++){
        s = (s0 + s1) / 2.0;
        if (s == s0 || s == s1) {
            break;
        }
        Real ratio0 = n0/(s + r0);
        Real ratio1 = z1/(s+1);
        g = ratio0*ratio0 + ratio1*ratio1 -1.0;
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
