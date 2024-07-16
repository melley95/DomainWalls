/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SPHERE_HPP_
#define SPHERE_HPP_

#include "Cell.hpp"
#include "Coordinates.hpp"
#include "MatterCCZ4RHS.hpp"
#include "ScalarField.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total no. components
#include "VarsTools.hpp"
#include "simd.hpp"


class Sphere
{
  public:
    //! A structure for the input params for scalar field properties and initial
    //! conditions
    struct params_t
    {
        
        std::array<double, CH_SPACEDIM>
            center;   //!< Centre of sphere
        double eta;
        double lambda;
        double R0;

    };

    //! The constructor
    Sphere(params_t a_params, double a_dx): m_dx(a_dx), m_params(a_params){}

    //! Function to compute the value of all the initial vars on the grid
    template <class data_t> void compute(Cell<data_t> current_cell) const
    {

        Coordinates<data_t> coords(current_cell, m_dx, m_params.center);

  

   

        // Store the initial values of the variables
        current_cell.store_vars(compute_phi(coords),c_phi);
        current_cell.store_vars(0.0 ,c_Pi);


    }



    template <class data_t>
    data_t compute_phi(Coordinates<data_t> coords) const
{
    data_t x = coords.x;
    data_t y = coords.y;
    data_t z = coords.z;
    data_t rho = sqrt(x*x+y*y+z*z);
    data_t out_phi = m_params.eta*tanh(sqrt(2.0*m_params.lambda)*m_params.eta*(rho-m_params.R0)/2.0);

    return out_phi;
}

  protected:
    double m_dx;
    const params_t m_params; //!< The matter initial condition params

};



#endif /* SPHERE_HPP_ */
