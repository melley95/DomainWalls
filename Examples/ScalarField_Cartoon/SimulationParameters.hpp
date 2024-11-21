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
#include "ArrayTools.hpp"
#include "BoostedBH.hpp"
#include "InitialScalarData_2D.hpp"
#include "Potential.hpp"
#include "PhiAndKExtractionTaggingCriterion.hpp"


#ifdef USE_AHFINDER
#include "AHInitialGuess.hpp"
#endif

class SimulationParameters : public SimulationParametersBase
{
  public:
    SimulationParameters(GRParmParse &pp) : SimulationParametersBase(pp)
    {
        readParams(pp);
        check_params();
    }

    void readParams(GRParmParse &pp)
    {
        pp.load("G_Newton", m_G_Newton);

        // Initial data
        pp.load("massA", bh1_params.mass, 0.);
        pp.load("momentumA", bh1_params.momentum, {0., 0.});
        pp.load("massB", bh2_params.mass, 0.);
        pp.load("momentumB", bh2_params.momentum, {0., 0.});

    

      

        // Initial data InitialScalarData_2D
        // pp.load("sf_phi0", init_SF_params.phi0, .1);
        pp.load("sf_eta", init_SF_params.eta, 0.0);
        pp.load("lambda", init_SF_params.lambda, 0.0);
        pp.load("R0", init_SF_params.R0, 10.0);
        pp.load("eps1", init_SF_params.eps1, 1.0);
        pp.load("center_SF", init_SF_params.centerSF, center);

        // Potential params
       potential_params.eta = init_SF_params.eta;
       potential_params.lambda = init_SF_params.lambda;
 

        pp.load("thresh_phi", threshold_phi, 0.0);
        pp.load("thresh_K", threshold_K, 0.0);



        // Get the centers of the BHs either explicitly or as
        // an offset (not both, or they will be offset from center
        // provided)
        std::array<double, CH_SPACEDIM> centerA, centerB;
        std::array<double, CH_SPACEDIM> offsetA, offsetB;
        pp.load("centerA", centerA, center);
        pp.load("centerB", centerB, center);
        pp.load("offsetA", offsetA, {0.0, 0.0});
        pp.load("offsetB", offsetB, {0.0, 0.0});
        // Do we want Weyl extraction, constraint norm
        // calculation?
        pp.load("activate_extraction", activate_extraction, false);
        pp.load("calculate_constraint_norms", calculate_constraint_norms,
                false);

 


        FOR(idir)
        {
            bh1_params.center[idir] = centerA[idir] + offsetA[idir];
            bh2_params.center[idir] = centerB[idir] + offsetB[idir];
        }

        

    #ifdef USE_AHFINDER
        pp.load("AH_initial_guess", AH_initial_guess, 0.5);
        
    
    #endif
    }
    void check_params()
    {
        
   
       
    }
    // tagging
    bool activate_extraction, calculate_constraint_norms;
 

 


    // For PhiAndK regridding
    double threshold_phi;
    double threshold_K; 



    // Collection of parameters necessary for initial conditions
    BoostedBH::params_t bh1_params;
    BoostedBH::params_t bh2_params;

    InitialScalarData_2D::params_t init_SF_params;
    Potential::params_t potential_params;

    extraction_params_t extraction_params_ADM;

    double m_G_Newton;

    #ifdef USE_AHFINDER
    double AH_initial_guess;
    #endif
};
#endif /* SIMULATIONPARAMETERS_HPP_ */
