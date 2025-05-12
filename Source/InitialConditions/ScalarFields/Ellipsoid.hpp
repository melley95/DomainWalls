/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef ELLIPSOID_HPP_
#define ELLIPSOID_HPP_

#include "MatterCCZ4RHS.hpp"
#include "ScalarField.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total no. components
#include "VarsTools.hpp"
#include "parstream.H" //gives pout
#include "simd.hpp"
#include <array>
#include <vector>

#include <iostream>
#include <cmath>
#include <limits>

//! Class which solves for the initial data for a spherically symmetric boson
//! star with phi^4 coupling
class Ellipsoid
{

  public:
    struct params_t
    {
        double eta;
        double lambda;
        double a;
        double b;
        double c;
        std::array<double, CH_SPACEDIM>
            centerSF;
  
    };
    //! The constructor
    Ellipsoid(params_t a_init_SF_params, double a_dx);


    //! Function to compute the value of all the initial vars on the grid
    template <class data_t> void compute(Cell<data_t> current_cell) const;

  
 

  protected:

    template <class data_t>
    data_t get_root(data_t r0, data_t r1, data_t z0, data_t z1, data_t z2, data_t g) const;
    const params_t m_init_SF_params;
 
    double m_dx;
};

#include "Ellipsoid.impl.hpp"

#endif /* ELLIPSOID_HPP_ */
