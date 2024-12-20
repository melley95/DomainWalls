/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// General includes common to most GR problems
#include "ScalarFieldLevel.hpp"
#include "BoxLoops.hpp"
#include "NanCheck.hpp"
#include "PositiveChiAndAlpha.hpp"
#include "SixthOrderDerivatives.hpp"
#include "TraceARemoval.hpp"

// For RHS update
#include "MatterCCZ4RHS.hpp"

// For constraints calculation
#include "NewMatterConstraints.hpp"

// For tag cells
#include "PhiAndKExtractionTaggingCriterion.hpp"

// Problem specific includes
#include "ComputePack.hpp"
#include "GammaCalculator.hpp"
#include "Spheroid.hpp"
#include "Potential.hpp"
#include "ScalarField.hpp"
#include "SetValue.hpp"
#include "Flat.hpp"

#include "MatterEnergy.hpp"
#include "FluxExtraction.hpp"
#include "AMRReductions.hpp"
#include "ExcisionDiagnostics.hpp"
#include "MovingPunctureGauge.hpp"

#include "ScalarExtraction.hpp"


#include "Weyl4.hpp"
#include "WeylExtraction.hpp"




// Things to do at each advance step, after the RK4 is calculated
void ScalarFieldLevel::specificAdvance()
{
    // Enforce trace free A_ij and positive chi and alpha
    BoxLoops::loop(
        make_compute_pack(TraceARemoval(),
                          PositiveChiAndAlpha(m_p.min_chi, m_p.min_lapse)),
        m_state_new, m_state_new, INCLUDE_GHOST_CELLS);

    // Check for nan's
    if (m_p.nan_check)
        BoxLoops::loop(
            NanCheck(m_dx, m_p.center, "NaNCheck in specific Advance"),
            m_state_new, m_state_new, EXCLUDE_GHOST_CELLS, disable_simd());
}

// Initial data for field and metric variables
void ScalarFieldLevel::initialData()
{
    CH_TIME("ScalarFieldLevel::initialData");
    if (m_verbosity)
        pout() << "ScalarFieldLevel::initialData " << m_level << endl;

    // First set everything to zero then initial conditions for scalar field -
    // here a Kerr BH and a scalar field profile
    BoxLoops::loop(
        make_compute_pack(SetValue(0.), Flat(),
                          Spheroid(m_p.initial_params, m_dx)),
        m_state_new, m_state_new, INCLUDE_GHOST_CELLS);

    fillAllGhosts();
    BoxLoops::loop(GammaCalculator(m_dx), m_state_new, m_state_new,
                   EXCLUDE_GHOST_CELLS);
}

#ifdef CH_USE_HDF5
// Things to do before outputting a checkpoint file
void ScalarFieldLevel::prePlotLevel()
{
    fillAllGhosts();
    Potential potential(m_p.pot_params);
    ScalarFieldWithPotential scalar_field(potential);
    BoxLoops::loop(
        MatterConstraints<ScalarFieldWithPotential>(
            scalar_field, m_dx, m_p.G_Newton, c_Ham, Interval(c_Mom1, c_Mom3)),
        m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);

    if (m_p.activate_extraction == 1 && m_p.calc_weyl == 1)
    {
        BoxLoops::loop(
                Constraints(m_dx, c_Ham, Interval(c_Mom1, c_Mom3)),
            m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
    }
}
#endif

// Things to do in RHS update, at each RK4 step
void ScalarFieldLevel::specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                       const double a_time)
{
    // Enforce trace free A_ij and positive chi and alpha
    BoxLoops::loop(
        make_compute_pack(TraceARemoval(),
                          PositiveChiAndAlpha(m_p.min_chi, m_p.min_lapse)),
        a_soln, a_soln, INCLUDE_GHOST_CELLS);

    // Calculate MatterCCZ4 right hand side with matter_t = ScalarField
    Potential potential(m_p.pot_params);
    ScalarFieldWithPotential scalar_field(potential);
    if (m_p.max_spatial_derivative_order == 4)
    {
        MatterCCZ4RHS<ScalarFieldWithPotential, MovingPunctureGauge,
                      FourthOrderDerivatives>
            my_ccz4_matter(scalar_field, m_p.ccz4_params, m_dx, m_p.sigma,
                           m_p.formulation, m_p.G_Newton);
        BoxLoops::loop(my_ccz4_matter, a_soln, a_rhs, EXCLUDE_GHOST_CELLS);
    }
    else if (m_p.max_spatial_derivative_order == 6)
    {
        MatterCCZ4RHS<ScalarFieldWithPotential, MovingPunctureGauge,
                      SixthOrderDerivatives>
            my_ccz4_matter(scalar_field, m_p.ccz4_params, m_dx, m_p.sigma,
                           m_p.formulation, m_p.G_Newton);
        BoxLoops::loop(my_ccz4_matter, a_soln, a_rhs, EXCLUDE_GHOST_CELLS);
    }
}

// Things to do at ODE update, after soln + rhs
void ScalarFieldLevel::specificUpdateODE(GRLevelData &a_soln,
                                         const GRLevelData &a_rhs, Real a_dt)
{
    // Enforce trace free A_ij
    BoxLoops::loop(TraceARemoval(), a_soln, a_soln, INCLUDE_GHOST_CELLS);
}

void ScalarFieldLevel::preTagCells()
{
    fillAllGhosts(VariableType::evolution, Interval(c_phi, c_phi));
    fillAllGhosts(VariableType::evolution, Interval(c_K, c_K));
}

void ScalarFieldLevel::computeTaggingCriterion(
    FArrayBox &tagging_criterion, const FArrayBox &current_state,
    const FArrayBox &current_state_diagnostics)
{
    BoxLoops::loop(
        PhiAndKExtractionTaggingCriterion(m_dx, m_p.threshold_phi, m_p.threshold_K, m_level, m_p.scalar_extraction_params, m_p.activate_extraction),
        current_state, tagging_criterion);
}
void ScalarFieldLevel::specificPostTimeStep()
{

     bool first_step = (m_time == 0.0);

     if (m_p.activate_scalar_extraction == 1)
    {
    
    int min_level = m_p.scalar_extraction_params.min_extraction_level();
    bool fill_ghosts = false;
    bool calculate_min_level = at_level_timestep_multiple(min_level);

    if (calculate_min_level){


    
    fillAllGhosts();
    Potential potential(m_p.pot_params);
    ScalarFieldWithPotential scalar_field(potential);
    BoxLoops::loop(
        MatterConstraints<ScalarFieldWithPotential>(
            scalar_field, m_dx, m_p.G_Newton, c_Ham, Interval(c_Mom1, c_Mom3)),
        m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
     BoxLoops::loop(MatterEnergy<ScalarFieldWithPotential>(scalar_field, m_dx, m_p.center),
                       m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);

   

    // Remove diagnostics of volume outside extraction regions -- TODO: multiple volumes

 
    BoxLoops::loop(
        ExcisionDiagnostics(m_dx, m_p.center, 0.0, 
                            m_p.scalar_extraction_params.extraction_radii[0]), //m_p.extraction_params.extraction_radii.size() -1 - i
        m_state_diagnostics, m_state_diagnostics, SKIP_GHOST_CELLS,
        disable_simd());

  
    if (m_level == min_level)
    {
  
    m_bh_amr.m_interpolator->refresh(fill_ghosts);

    AMRReductions<VariableType::diagnostic> amr_reductions(m_bh_amr);
    double rho_sum = amr_reductions.sum(c_rhoLL);
    double source_sum = amr_reductions.sum(c_source);




    SmallDataIO integral_file("data/VolumeIntegrals_r"+std::to_string((int)m_p.scalar_extraction_params.extraction_radii[0]), m_dt, m_time,
                                  m_restart_time, SmallDataIO::APPEND,
                                  first_step);
    // remove any duplicate data if this is post restart
        integral_file.remove_duplicate_time_data();
        std::vector<double> data_for_writing = {rho_sum, source_sum};
        // write data
        if (first_step)
        {
            integral_file.write_header_line({"rho", "source"});
        }
        integral_file.write_time_data_line(data_for_writing);


  
    // Now refresh the interpolator and do the interpolation
      //  m_bh_amr.m_interpolator->refresh(fill_ghosts);
        m_bh_amr.fill_multilevel_ghosts(VariableType::diagnostic,
                                        Interval(c_flux, c_flux));
        FluxExtraction flux_extraction(m_p.scalar_extraction_params, m_dt, m_time,
                                       first_step, m_restart_time);
        flux_extraction.execute_query(m_bh_amr.m_interpolator);

        ScalarExtraction phi_extraction(
                        m_p.scalar_extraction_params, m_dt, m_time, first_step,
                        m_restart_time);
        phi_extraction.execute_query(m_gr_amr.m_interpolator);

    }


     if (m_p.calc_weyl){
            
            int min_level_weyl = m_p.extraction_params.min_extraction_level();
            // Populate the Weyl Scalar values on the grid
            fillAllGhosts();
            BoxLoops::loop(
                Weyl4(m_p.extraction_params.center, m_dx, m_p.formulation),
                m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);

            // Do the extraction on the min extraction level
            if (m_level == min_level_weyl)
            {
                CH_TIME("WeylExtraction");
                // Now refresh the interpolator and do the interpolation
                // fill ghosts manually to minimise communication
                bool fill_ghosts = false;
                m_gr_amr.m_interpolator->refresh(fill_ghosts);
                m_gr_amr.fill_multilevel_ghosts(
                    VariableType::diagnostic, Interval(c_Weyl4_Re, c_Weyl4_Im),
                    min_level);
                WeylExtraction my_extraction(m_p.extraction_params, m_dt,
                                             m_time, first_step,
                                             m_restart_time);
                my_extraction.execute_query(m_gr_amr.m_interpolator);
            }
     }
    }

        if (m_p.calculate_constraint_norms)
    {
        fillAllGhosts();
        BoxLoops::loop(Constraints(m_dx, c_Ham, Interval(c_Mom1, c_Mom3)),
                       m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
        if (m_level == 0)
        {
            AMRReductions<VariableType::diagnostic> amr_reductions(m_gr_amr);
            double L2_Ham = amr_reductions.norm(c_Ham);
            double L2_Mom = amr_reductions.norm(Interval(c_Mom1, c_Mom3));
            SmallDataIO constraints_file(m_p.data_path + "constraint_norms",
                                         m_dt, m_time, m_restart_time,
                                         SmallDataIO::APPEND, first_step);
            constraints_file.remove_duplicate_time_data();
            if (first_step)
            {
                constraints_file.write_header_line({"L^2_Ham", "L^2_Mom"});
            }
            constraints_file.write_time_data_line({L2_Ham, L2_Mom});
        }
    }

    }
    

    #ifdef USE_AHFINDER
    // if print is on and there are Diagnostics to write, calculate them!
    if (m_p.AH_activate && m_level == m_p.AH_params.level_to_run)
        m_bh_amr.m_ah_finder.solve(m_dt, m_time, m_restart_time);
    #endif

}