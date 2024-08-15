/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP_
#define SIMULATIONPARAMETERS_HPP_

// General includes
#include "GRParmParse.hpp"
#include "SimulationParametersBase.hpp"

// Problem specific includes:
#include "Spheroid.hpp"
//#include "KerrBH.hpp"
#include "Potential.hpp"
#include "PhiAndKTaggingCriterion.hpp"

class SimulationParameters : public SimulationParametersBase
{
  public:
    SimulationParameters(GRParmParse &pp) : SimulationParametersBase(pp)
    {
        // read the problem specific params
        read_params(pp);
        check_params();
    }

    void read_params(GRParmParse &pp)
    {
        // Initial scalar field data
        initial_params.center =
            center; // already read in SimulationParametersBase
        pp.load("G_Newton", G_Newton,
                0.0); // for now the example neglects backreaction
      
        pp.load("R0", initial_params.R0, 0.0);
        pp.load("sf_eta", initial_params.eta, 0.0);
        pp.load("sf_lambda", initial_params.lambda, 0.0);

        pp.load("a0", initial_params.a, 1.0);
        pp.load("b0", initial_params.b, 1.0);
        pp.load("c0", initial_params.c, 1.0);

        pot_params.eta = initial_params.eta;
        pot_params.lambda = initial_params.lambda;

        pp.load("thresh_phi", threshold_phi, 0.0);
        pp.load("thresh_K", threshold_K, 0.0);

        pp.load("activate_extraction", activate_extraction, false);

        #ifdef USE_AHFINDER
        pp.load("AH_initial_guess", AH_initial_guess, 0.5);
        #endif




    }

    void check_params()
    {
   
    }

    // Initial data for matter and potential and BH
    double G_Newton;
    Potential::params_t pot_params;
    Spheroid::params_t initial_params;
    
    double threshold_phi;
    double threshold_K; 

    bool activate_extraction;
  //  PhiAndKTaggingCriterion::params_t tag_crit;

    #ifdef USE_AHFINDER
    double AH_initial_guess;
    #endif



};

#endif /* SIMULATIONPARAMETERS_HPP_ */
