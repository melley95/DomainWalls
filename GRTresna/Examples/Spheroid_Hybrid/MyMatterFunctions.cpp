#include "ScalarField.hpp"
Real ScalarField::my_potential_function(const Real &phi_here) const
{
    return m_matter_params.lambda*0.25*pow(pow(phi_here, 2.0)-pow(m_matter_params.eta, 2.0), 2.0);
}
Real ScalarField::my_phi_function(const RealVect &loc) const
{
    Real r2 = max(loc[0] * loc[0] + loc[1] * loc[1] + loc[2] * loc[2], 1e-12);

    Real x = loc[0];
    Real y = loc[1];
    Real z = loc[2];

    Real a = m_matter_params.a;
    Real b = 1.0/sqrt(1.0 - pow(m_matter_params.e, 2.0));

    Real R = (a*b)/sqrt((a*a*x*x + b*b*y*y + a*a*z*z)/r2);

    return m_matter_params.eta*tanh(sqrt(0.5*m_matter_params.lambda)*m_matter_params.eta*(sqrt(r2)-R));
}
Real ScalarField::my_Pi_function(const RealVect &loc) const
{
    
    return 0.0;
}