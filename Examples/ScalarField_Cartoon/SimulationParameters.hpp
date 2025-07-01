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
#include "Spheroid.hpp"
#include "Potential.hpp"



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
        pp.load("a", init_SF_params.a, 1.0);
        pp.load("b", init_SF_params.b, 1.0);
        pp.load("lambda", potential_params.lambda, 0.0);
        pp.load("sf_eta", potential_params.eta, 0.0);
        pp.load("center_SF", init_SF_params.centerSF, center);

        // Potential params
        init_SF_params.eta = potential_params.eta;
        init_SF_params.lambda = potential_params.lambda ;


 

   //     pp.load("thresh_rho", threshold_rho, 0.0);
    //   pp.load("thresh_K", threshold_K, 0.0);
        pp.load("thresh_chi", threshold_chi, 0.0);
        pp.load("thresh_phi", threshold_phi, 0.0);

   //     pp.load("rebound_time", rebound_time, 0.0);
    //    pp.load("rebound_thresh_phi", threshold_phi_rebound, 0.0);
     //   pp.load("rebound_thresh_chi", threshold_chi_rebound, 0.0);
      



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

        pp.load("AH1_pos", AH1_pos, {0.0, 0.0});
        pp.load("AH2_pos", AH2_pos, {0.0, 0.0});
        pp.load("AH3_pos", AH3_pos, {0.0, 0.0});

        pp.load("AH1_a", AH1_a, 1.0);
        pp.load("AH1_b", AH1_b, 1.0);

        pp.load("AH2_a", AH2_a, 1.0);
        pp.load("AH2_b", AH2_b, 1.0);

        pp.load("AH3_a", AH3_a, 1.0);
        pp.load("AH3_b", AH3_b, 1.0);


        pp.load("origin_x", origin_x, 0.0);
        pp.load("origin_y", origin_y, 0.0);
        pp.load("L", L, 0.0);

        pp.load("max_ref_level", max_ref_level, 0);
       
       // pp.load("AH2_shape", AH2_shape, {0.0, 0.0});
       // pp.load("AH3_shape", AH3_shape, {0.0, 0.0});


   //     pp.load("excise", excise, false);
    //    pp.load("r_excise", r_excise);

      



       
    //   pp.load("ref_times", ref_times);
     //  pp.load("ref_levels", ref_levels);
   //    pp.load("threshold_rho", threshold_rho);
       

        

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
   // double threshold_rho;
  // double threshold_K;
    double threshold_chi; 
    double threshold_phi;

 //   double rebound_time;
  //  double threshold_chi_rebound; 
   // double threshold_phi_rebound;



    // Collection of parameters necessary for initial conditions
    BoostedBH::params_t bh1_params;
    BoostedBH::params_t bh2_params;

    Spheroid::params_t init_SF_params;
    Potential::params_t potential_params;

    extraction_params_t extraction_params_ADM;

    double m_G_Newton;

    #ifdef USE_AHFINDER
    double AH_initial_guess;
    #endif

    bool excise;
    double r_excise;


    std::array<double,CH_SPACEDIM> AH1_pos;
    std::array<double,CH_SPACEDIM> AH2_pos;
    std::array<double,CH_SPACEDIM> AH3_pos;

    double AH1_a;
    double AH1_b;

    double AH2_a;
    double AH2_b;

    double AH3_a;
    double AH3_b;


    double origin_x;
    double origin_y;

    double L;

    int max_ref_level;







   // std::array<double, 10> ref_times;
   // std::array<int, 10> ref_levels;

   
  //  std::array<double, 10> ref_rho;

};
#endif /* SIMULATIONPARAMETERS_HPP_ */
