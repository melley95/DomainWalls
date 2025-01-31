/* GRTresna
 * Copyright 2024 The GRTL Collaboration.
 * Please refer to LICENSE in GRTresna's root directory.
 */

#ifndef GRCHOMBOVARIABLES_HPP
#define GRCHOMBOVARIABLES_HPP

#include "ArrayTools.hpp"
#include "ParityDefinitions.hpp"

// assign an enum to each variable
enum
{
    c_chi,

    c_h11,
    c_h12,
    c_h22,
    c_hww,


    c_K,

    c_A11,
    c_A12,
    c_A22,

    c_Aww,


    c_Theta,

    c_Gamma1,
    c_Gamma2,


    c_lapse,

    c_shift1,
    c_shift2,


    c_B1,
    c_B2,


    c_phi, // matter field added
    c_Pi,  //(minus) conjugate momentum

    NUM_GRCHOMBO_VARS
};

namespace GRChomboVariables
{
static constexpr char const *variable_names[NUM_GRCHOMBO_VARS] = {
    "chi",

    "h11",    "h12",  "h22",

    "hww",

    "K",

    "A11",    "A12",   "A22",

    "Aww",

    "Theta",

    "Gamma1", "Gamma2", 

    "lapse",

    "shift1", "shift2", 

    "B1",     "B2",     

    "phi",    "Pi"};

static constexpr std::array<int, NUM_GRCHOMBO_VARS> const vars_parity = {
    EVEN,   
    EVEN,   ODD_XY,  EVEN,   
    EVEN,
    EVEN, 
    EVEN, ODD_XY, EVEN,
    EVEN,
    EVEN,  
    ODD_X, ODD_Y,  
    EVEN, 
    ODD_X, ODD_Y,  
    ODD_X, ODD_Y,  
    EVEN,  EVEN};

} // namespace GRChomboVariables

#endif /* GRCHOMBOVARIABLES_HPP */