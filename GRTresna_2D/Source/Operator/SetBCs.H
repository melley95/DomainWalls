/* GRTresna                                                                          
* Copyright 2024 The GRTL Collaboration.
* Please refer to LICENSE in GRTresna's root directory.
*/

#ifndef _SETBCS_H_
#define _SETBCS_H_

#include <iostream>
using std::cerr;

#include "AMRMultiGrid.H"
#include "BoundaryConditions.hpp"
#include "DimensionDefinitions.hpp"
#include "DisjointBoxLayout.H"
#include "FArrayBox.H"
#include "ParmParse.H"
#include "VariableCoeffPoissonOperator.H"

// BC class for the linear solver steps
class GlobalBCRS
{
  public:
    static BoundaryConditions::params_t
        s_boundary_params; // set boundaries in each dir
    static bool s_areBCsParsed;
};

// This is the function that the solver looks for to apply the BCs to the
// variables in a_state
extern void ParseBC(FArrayBox &a_state, const Box &a_valid,
                    const ProblemDomain &a_domain, Real a_dx,
                    bool a_homogeneous);

#endif
