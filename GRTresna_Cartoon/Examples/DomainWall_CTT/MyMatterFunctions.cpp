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

    Real b = m_matter_params.b;
    Real a = b/sqrt(1.0 - pow(m_matter_params.e, 2.0));

    Real min_theta = 0.0;
    Real min_dist = std::numeric_limits<double>::max();
    Real best_theta = 0.0;

    int steps = 10000;
    for (int i = 0; i <= steps; ++i) {
        double theta = 2.0 * M_PI * i / steps;
        Real x = a * std::cos(theta);
        Real y = b * std::sin(theta);
        Real dist= (x - x_p) * (x - x_p) + (y - y_p) * (y - y_p);
       
        if (dist < min_dist) {
            min_dist = dist;
            best_theta = theta;
        }
    }
    int fac = 1;
    if (x_p*x_p/(a*a) + y_p*y_p/(b*b) < 1 ){
        fac = -1;
    }
    


    return m_matter_params.eta*tanh(sqrt(0.5*m_matter_params.lambda)*m_matter_params.eta*(fac*sqrt(min_dist)));
}

Real ScalarField::my_Pi_function(const RealVect &loc) const
{
    
    return 0.0;
}
