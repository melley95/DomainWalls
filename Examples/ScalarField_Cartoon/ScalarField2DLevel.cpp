/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include "ScalarField2DLevel.hpp"
// #include "BinaryPunctureTaggingCriterion.hpp"
//  #include "BoostedPunctureTrackerTaggingCriterion.hpp"
#include "AMRReductions.hpp"
#include "BoxLoops.hpp"
#include "CCZ4Cartoon.hpp"
#include "ComputePack.hpp"
#include "ConstraintsCartoon.hpp"
#include "MovingPunctureGauge.hpp"
#include "NanCheck.hpp"
#include "PositiveChiAndAlpha.hpp"
#include "SetValue.hpp"
#include "TraceARemovalCartoon.hpp"
#include "WeylExtraction.hpp"
#include "WeylOmScalar.hpp"

#include "PhiAndKExtractionTaggingCriterion.hpp"

// Initial data
//#include "HeadOn2D.hpp"
#include "InitialScalarData_2D.hpp"
//#include "IsotropicBoostedBH_bk.hpp"


// Reference connection
#include "ADMQuantities.hpp"
#include "ADMQuantitiesExtraction.hpp"
#include "GammaCartoonCalculator.hpp"

void ScalarField2DLevel::specificAdvance()
{
    // Enforce the trace free A_ij condition and positive chi and alpha
    BoxLoops::loop(
        make_compute_pack(TraceARemovalCartoon(), PositiveChiAndAlpha()),
        m_state_new, m_state_new, INCLUDE_GHOST_CELLS);

    // Check for nan's
    if (m_p.nan_check)
        BoxLoops::loop(
            NanCheck(m_dx, m_p.center, "NaNCheck in specific Advance: "),
            m_state_new, m_state_new, EXCLUDE_GHOST_CELLS, disable_simd());
}

void ScalarField2DLevel::initialData()
{
    CH_TIME("ScalarField2DLevel::initialData");
    if (m_verbosity)
        pout() << "ScalarField2DLevel::initialData " << m_level << endl;

    // When changing class here, don't forget to change potential if necessary
    // double spacing = .01;
    // Oscilloton oscilloton(m_p.oscilloton_params, m_dx, spacing);
    // BoxLoops::loop(make_compute_pack(SetValue(0.0), oscilloton), m_state_new,
    //                m_state_new, INCLUDE_GHOST_CELLS, disable_simd());

    // When changing class here, don't forget to change potential if necessary
    InitialScalarData_2D initialscalardata_2D(m_p.init_SF_params, m_dx);
    BoxLoops::loop(make_compute_pack(SetValue(0.0), initialscalardata_2D),
                   m_state_new, m_state_new, INCLUDE_GHOST_CELLS,
                   disable_simd());

    // When changing class here, don't forget to change potential if necessary
    // ScalarBubble_2D scalarbubble2D(m_p.bubble_params, m_p.potential_params,
    //                                m_dx);
    // BoxLoops::loop(make_compute_pack(SetValue(0.0), scalarbubble2D),
    //                m_state_new, m_state_new, INCLUDE_GHOST_CELLS,
    //                disable_simd());

    fillAllGhosts();

    BoxLoops::loop(GammaCartoonCalculator(m_dx), m_state_new, m_state_new,
                   EXCLUDE_GHOST_CELLS, disable_simd());
}

// Things to do before a plot level - need to calculate the Weyl scalars
void ScalarField2DLevel::prePlotLevel()
{
    fillAllGhosts();
    Potential potential(m_p.potential_params);
    BoxLoops::loop(Constraints<Potential>(m_dx, potential, m_p.m_G_Newton),
                   m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
}

void ScalarField2DLevel::specificEvalRHS(GRLevelData &a_soln,
                                         GRLevelData &a_rhs,
                                         const double a_time)
{
    // Enforce positive chi and alpha and trace free A
    BoxLoops::loop(
        make_compute_pack(TraceARemovalCartoon(), PositiveChiAndAlpha()),
        a_soln, a_soln, INCLUDE_GHOST_CELLS);

    // Calculate CCZ4 right hand side
    Potential potential(m_p.potential_params);
    BoxLoops::loop(
        CCZ4Cartoon<MovingPunctureGauge, FourthOrderDerivatives, Potential>(
            m_p.ccz4_params, m_dx, m_p.sigma, potential, m_p.m_G_Newton,
            m_p.formulation),
        a_soln, a_rhs, EXCLUDE_GHOST_CELLS);
}

void ScalarField2DLevel::specificUpdateODE(GRLevelData &a_soln,
                                           const GRLevelData &a_rhs, Real a_dt)
{
    // Enforce the trace free A_ij condition
    BoxLoops::loop(TraceARemovalCartoon(), a_soln, a_soln, INCLUDE_GHOST_CELLS);
}

void ScalarField2DLevel::preTagCells()
{

}

void ScalarField2DLevel::computeTaggingCriterion(FArrayBox &tagging_criterion,
                                                 const FArrayBox &current_state
                                                 )
{

    BoxLoops::loop(
        PhiAndKExtractionTaggingCriterion(m_dx, m_p.threshold_phi, m_p.threshold_K, m_level, m_p.extraction_params, m_p.activate_extraction),
        current_state, tagging_criterion);


}

void ScalarField2DLevel::specificPostTimeStep()
{
    CH_TIME("ScalarField2DLevel::specificPostTimeStep");

    bool first_step =
        (m_time == 0.); // this form is used when 'specificPostTimeStep' was
                        // called during setup at t=0 from Main

 #ifdef USE_AHFINDER
    // if print is on and there are Diagnostics to write, calculate them!
    if (m_p.AH_activate && m_level == m_p.AH_params.level_to_run)
        m_bh_amr.m_ah_finder.solve(m_dt, m_time, m_restart_time);
 #endif

    if (m_p.calculate_constraint_norms)
    {
        fillAllGhosts();
        Potential potential(m_p.potential_params);
        BoxLoops::loop(Constraints<Potential>(m_dx, potential, m_p.m_G_Newton),
                       m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
        if (m_level == 0)
        {
            bool first_step = (m_time == 0.);
            AMRReductions<VariableType::diagnostic> amr_reductions(m_bh_amr);
            double L2_Ham = amr_reductions.sum(c_Ham);
            double L2_Mom = amr_reductions.sum(Interval(c_Mom, c_Mom));
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

    int adm_min_level = 0;
    bool calculate_adm = at_level_timestep_multiple(adm_min_level);
    if (calculate_adm)
    {
        Potential potential(m_p.potential_params);
        if (!m_p.calculate_constraint_norms)
        {
            fillAllGhosts();
            BoxLoops::loop(
                Constraints<Potential>(m_dx, potential, m_p.m_G_Newton),
                m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
        }
        if (m_level == 0)
        {
            AMRReductions<VariableType::diagnostic> amr_reductions(m_bh_amr);
            double M_ADM = amr_reductions.sum(c_rho_ADM);
            SmallDataIO M_ADM_file(m_p.data_path + "M_ADM", m_dt, m_time,
                                   m_restart_time, SmallDataIO::APPEND,
                                   first_step);
            M_ADM_file.remove_duplicate_time_data();
            if (first_step)
            {
                M_ADM_file.write_header_line({"M_ADM"});
            }
            M_ADM_file.write_time_data_line({M_ADM});
        }
    }

    if (m_p.activate_extraction == 1)
    {
        // Deleted Thomas' ADM extraction code for now, still available in
        // previous commits.
    }
}
