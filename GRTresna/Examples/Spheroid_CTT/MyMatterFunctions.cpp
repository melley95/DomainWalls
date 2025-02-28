#include "ScalarField.hpp"
Real ScalarField::my_potential_function(const Real &phi_here) const
{
    return m_matter_params.lambda*0.25*pow(pow(phi_here, 2.0)-pow(m_matter_params.eta, 2.0), 2.0);
}
Real ScalarField::my_phi_function(const RealVect &loc) const
{
    Real rho2 = max(loc[0] * loc[0] + loc[1] * loc[1] + loc[2] * loc[2],1e-12);
    Real r2 = max(loc[0] * loc[0] + loc[1] * loc[1], 1e-12);
    Real cos2phi = (loc[1] * loc[1])/r2;
    Real sin2phi = (loc[0] * loc[0])/r2;
    Real sin2theta = r2/ rho2;
    Real cos2theta = loc[2] * loc[2] / rho2;
    Real R = m_matter_params.R0/sqrt((cos2phi+pow(m_matter_params.eps1, -2)*sin2phi)*sin2theta+cos2theta);
    return m_matter_params.eta*tanh(sqrt(2.0*m_matter_params.lambda)*m_matter_params.eta*(sqrt(rho2)-R)/2.0);
}
Real ScalarField::my_Pi_function(const RealVect &loc) const
{
    
    return 0.0;
}