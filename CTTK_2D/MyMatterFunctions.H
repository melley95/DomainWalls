/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef MYMATTERFUNCTIONS_HPP
#define MYMATTERFUNCTIONS_HPP

#include "SetPsiAndAij.H"

Real my_phi_function(const RealVect loc, const PoissonParameters &a_params,
                     const RealVect &a_dx)
{

    Real xx = loc[0];
    Real yy = loc[1];
    Real zz = loc[2];

    Real rho2 = std::max(xx*xx + yy*yy + zz*zz,1e-12);
    Real r2 = std::max(xx*xx + yy*yy, 1e-12);

    Real cos2phi = xx*xx/r2;
    Real sin2phi = yy*yy/r2;
    Real sin2theta = r2/ rho2;
    Real cos2theta = zz*zz / rho2;

    Real RR = a_params.R0/sqrt((cos2phi+pow(a_params.eps1, -2)*sin2phi)*sin2theta+pow(a_params.eps2,-2)*cos2theta);
    
    
    Real out_phi = a_params.eta*tanh(sqrt(2.0*a_params.lambda)*a_params.eta*(sqrt(rho2)-RR)/2.0);


    
    return out_phi;
   
   
   
   // RealVect L = a_params.domainLength;


}

Real my_Pi_function(const RealVect loc, const PoissonParameters &a_params,
                    const RealVect &a_dx)
{

    
    // sine waves in 3 directions - use with periodic BCs
    return 0.0;
}

Real my_potential_function(const Real &phi_here,
                           const PoissonParameters &a_params)
{
    Real V_of_phi;
    

    V_of_phi = a_params.lambda*0.25*pow(pow(phi_here, 2.0)-pow(a_params.eta, 2.0), 2.0);
    




    return V_of_phi;
}

Real my_rho_function(const IntVect a_iv, const RealVect &a_dx,
                     FArrayBox &a_multigrid_vars_box,
                     const PoissonParameters &a_params)
{

    RealVect loc;
    get_loc(loc, a_iv, a_dx, a_params);

    Real psi_reg = a_multigrid_vars_box(a_iv, c_psi_reg);
    Real psi_bh = set_binary_bh_psi(loc, a_params);
    Real psi_0 = psi_reg + psi_bh;

    Real V_of_phi =
        my_potential_function(a_multigrid_vars_box(a_iv, c_phi_0), a_params);
    Real Pi_0 = a_multigrid_vars_box(a_iv, c_Pi_0);
    Tensor<1, Real, SpaceDim> d1_phi =
        get_d1(a_iv, a_multigrid_vars_box, a_dx, c_phi_0);
    Real d1_phi_squared = 0;
    FOR1(i) { d1_phi_squared += d1_phi[i] * d1_phi[i]; }

    return 0.5 * Pi_0 * Pi_0 + V_of_phi +
           0.5 * d1_phi_squared * pow(psi_0, -4.0);
}

void my_Si_function(Tensor<1, Real> &Si, const IntVect a_iv,
                    const RealVect &a_dx, FArrayBox &a_multigrid_vars_box,
                    const PoissonParameters &a_params)
{

    RealVect loc;
    get_loc(loc, a_iv, a_dx, a_params);

    Real Pi_0 = a_multigrid_vars_box(a_iv, c_Pi_0);
    Tensor<1, Real, SpaceDim> d1_phi =
        get_d1(a_iv, a_multigrid_vars_box, a_dx, c_phi_0);

    FOR1(i) { Si[i] = -Pi_0 * d1_phi[i]; }
}
#endif /* MYMATTERFUNCTIONS_HPP */
