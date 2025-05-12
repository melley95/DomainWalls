/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SPHEROID_HPP_
#define SPHEROID_HPP_

#include "CCZ4CartoonVars.hpp"
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
class Spheroid
{

  public:
    struct params_t
    {
        double eta;
        double lambda;
        double a;
        double b;
        std::array<double, CH_SPACEDIM>
            centerSF;
  
    };
    //! The constructor
    Spheroid(params_t a_init_SF_params, double a_dx);


    //! Function to compute the value of all the initial vars on the grid
    template <class data_t> void compute(Cell<data_t> current_cell) const;

  
 

  protected:

    template <class data_t>
    data_t get_root(data_t r0, data_t z0, data_t z1, data_t g) const;
    const params_t m_init_SF_params;
 
    double m_dx;
};

#include "Spheroid.impl.hpp"

#endif /* SPHEROID_HPP_ */
