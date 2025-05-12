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
#include "Ellipsoid.hpp"
//#include "KerrBH.hpp"
#include "Potential.hpp"
#include "PhiAndKExtractionTaggingCriterion.hpp"

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
     //   initial_params.center =
      //      center; // already read in SimulationParametersBase
        pp.load("G_Newton", G_Newton,
                0.0); // for now the example neglects backreaction
      
     
        pp.load("sf_eta", pot_params.eta, 0.0);
        pp.load("sf_lambda", pot_params.lambda, 0.0);

        pp.load("a", initial_params.a, 1.0);
        pp.load("b", initial_params.b, 1.0);
        pp.load("c", initial_params.c, 1.0);
        pp.load("center_SF", initial_params.centerSF, center);

        initial_params.eta = pot_params.eta;
        initial_params.lambda = pot_params.lambda;

       

        pp.load("thresh_phi", threshold_phi, 0.0);
        pp.load("thresh_chi", threshold_chi, 0.0);
        pp.load("r_limit", r_limit, 0.0);

        // pp.load("activate_extraction", activate_extraction, false);

        #ifdef USE_AHFINDER
        pp.load("AH_initial_guess", AH_initial_guess, 0.5);
        #endif

        pp.load("calculate_weyl", calc_weyl, false);

        pp.load("constraint_norms", calculate_constraint_norms, false);




    }

    void check_params()
    {
   
    }

    // Initial data for matter and potential and BH
    double G_Newton;
    Potential::params_t pot_params;
    Ellipsoid::params_t initial_params;
    
    double threshold_phi;
    double threshold_chi; 

    double r_limit;

    bool calc_weyl, calculate_constraint_norms; // activate_extraction , 


  //  PhiAndKTaggingCriterion::params_t tag_crit;

    #ifdef USE_AHFINDER
    double AH_initial_guess;
    #endif



};

#endif /* SIMULATIONPARAMETERS_HPP_ */
