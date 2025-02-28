/* GRTresna
 * Copyright 2024 The GRTL Collaboration.
 * Please refer to LICENSE in GRTresna's root directory.
 */

#ifndef CTT_HPP_
#error "This file should only be included through CTT.hpp"
#endif

#include "DimensionDefinitions.hpp"
#include "GRParmParse.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"

template <typename matter_t> struct CTT<matter_t>::params_t
{
    int sign_of_K;
    bool use_compact_Vi_ansatz;
    Real regularised_part_psi;
    bool deactivate_zero_mode;
};

template <typename matter_t>
CTT<matter_t>::CTT(params_t a_method_params, matter_t *a_matter,
                                 PsiAndAijFunctions *a_psi_and_Aij_functions,
                                 int a_numLevels,
                                 const std::array<double, SpaceDim> a_center,
                                 Real a_G_Newton)
    : m_method_params(a_method_params), matter(a_matter), G_Newton(a_G_Newton),
      numLevels(a_numLevels), psi_and_Aij_functions(a_psi_and_Aij_functions),
      center(a_center)
{
}

template <typename matter_t>
void CTT<matter_t>::read_params(GRParmParse &pp,
                                       params_t &a_method_params)
{
    pp.load("sign_of_K", a_method_params.sign_of_K, 1);
    pp.load("use_compact_Vi_ansatz", a_method_params.use_compact_Vi_ansatz,
            false);
    pp.load("regularised_part_psi", a_method_params.regularised_part_psi, 1.0);
    pp.load("deactivate_zero_mode", a_method_params.deactivate_zero_mode,
            false);
}

template <typename matter_t>
void CTT<matter_t>::solve_analytic(
    LevelData<FArrayBox> *a_multigrid_vars, LevelData<FArrayBox> *a_rhs,
    const RealVect &a_dx)
{
    DerivativeOperators derivs(a_dx);
    // Iterate through the boxes in turn
    DataIterator dit = a_rhs->dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
    {
        FArrayBox &multigrid_vars_box = (*a_multigrid_vars)[dit()];
        FArrayBox &rhs_box = (*a_rhs)[dit()];
        Box unghosted_box = rhs_box.box();

        // Iterate through the interior of boxes
        // (ghosts need to be filled later due to gradient terms)
        BoxIterator bit(unghosted_box);
        for (bit.begin(); bit.ok(); ++bit)
        {
            // work out location on the grid
            IntVect iv = bit();
            RealVect loc;
            Grids::get_loc(loc, iv, a_dx, center);

            // Calculate the actual value of psi including BH part
            Real psi_reg = multigrid_vars_box(iv, c_psi_reg);

            Real psi_0 = psi_reg;
            Real laplacian_psi_reg;
            derivs.scalar_Laplacian(laplacian_psi_reg, iv, multigrid_vars_box,
                                    c_psi_reg);

         

      



            // be careful if at a point K = 0, may have discontinuity
            multigrid_vars_box(iv, c_K_0) =
                0.0;

            // set values for \bar Aij_0
            multigrid_vars_box(iv, c_A11_0) = 0.0;
            multigrid_vars_box(iv, c_A22_0) = 0.0;
            multigrid_vars_box(iv, c_A33_0) = 0.0;
            multigrid_vars_box(iv, c_A12_0) = 0.0;
            multigrid_vars_box(iv, c_A13_0) = 0.0;
            multigrid_vars_box(iv, c_A23_0) = 0.0;
        }
    }
}

template <typename matter_t>
void CTT<matter_t>::set_elliptic_terms(
    LevelData<FArrayBox> *a_multigrid_vars, LevelData<FArrayBox> *a_rhs,
    RefCountedPtr<LevelData<FArrayBox>> a_aCoef,
    RefCountedPtr<LevelData<FArrayBox>> a_bCoef, const RealVect &a_dx)
{
    DerivativeOperators derivs(a_dx);
    DataIterator dit = a_rhs->dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
    {
        FArrayBox &multigrid_vars_box = (*a_multigrid_vars)[dit()];
        FArrayBox &rhs_box = (*a_rhs)[dit()];
        FArrayBox &aCoef_box = (*a_aCoef)[dit()];
        FArrayBox &bCoef_box = (*a_bCoef)[dit()];
        // JCAurre: Initialise rhs=0, aCoef=0 and bCoef=1 for all constraint
        // variables
        for (int comp = 0; comp < NUM_CONSTRAINT_VARS; comp++)
        {
            rhs_box.setVal(0.0, comp);
            aCoef_box.setVal(0.0, comp);
            bCoef_box.setVal(1.0, comp);

            // this prevents small amounts of noise in the sources
            // activating the zero modes - (Garfinkle trick) see 2207.03125
            if (m_method_params.deactivate_zero_mode)
            {
                Real small_number = 1e-10;
                aCoef_box.setVal(-small_number, comp);
            }
        }
        Box unghosted_box = rhs_box.box();
        BoxIterator bit(unghosted_box);
        for (bit.begin(); bit.ok(); ++bit)
        {

            // work out location on the grid
            IntVect iv = bit();
            RealVect loc;
            Grids::get_loc(loc, iv, a_dx, center);

            // Calculate the actual value of psi including BH part
            Real psi_reg = multigrid_vars_box(iv, c_psi_reg);
   
            Real psi_0 = psi_reg;
            Real laplacian_psi_reg;
            derivs.scalar_Laplacian(laplacian_psi_reg, iv, multigrid_vars_box,
                                    c_psi_reg);

           
            // Compute emtensor components
            const auto emtensor =
                matter->compute_emtensor(iv, a_dx, multigrid_vars_box);

  

            // rhs terms, K is set to cancel matter terms only
            rhs_box(iv, c_psi) =
            - 2.0 * M_PI * G_Newton * emtensor.rho * pow(psi_reg, 5.0) - laplacian_psi_reg;

          
            // now set the values in the box
            rhs_box(iv, c_V1) = 0.0;
            rhs_box(iv, c_V2) = 0.0;
            rhs_box(iv, c_V3) = 0.0;

            // Periodic: Use ansatz B.3 in B&S (p547)
            // Non-periodic: Compact ansatz B.7 in B&S (p547)
        
            rhs_box(iv, c_U) =
                    0.0;
          

            // add the aCoef term
            aCoef_box(iv, c_psi) += 10.0 * M_PI * G_Newton * pow(psi_0, 4.0) * emtensor.rho;
        }
    }
}

template <typename matter_t>
void CTT<matter_t>::initialise_method_vars(
    LevelData<FArrayBox> &a_multigrid_vars, const RealVect &a_dx) const
{
    CH_assert(a_multigrid_vars.nComp() == NUM_MULTIGRID_VARS);

    DataIterator dit = a_multigrid_vars.dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
    {
        // These contain the vars in the boxes, set them all to zero
        FArrayBox &multigrid_vars_box = a_multigrid_vars[dit()];
        for (int comp = 0; comp < NUM_MULTIGRID_VARS; comp++)
        {
            multigrid_vars_box.setVal(0.0, comp);
        }

        // Iterate over the box and set non zero comps
        Box ghosted_box = multigrid_vars_box.box();
        BoxIterator bit(ghosted_box);
        for (bit.begin(); bit.ok(); ++bit)
        {

            // work out location on the grid
            IntVect iv = bit();
            RealVect loc;
            Grids::get_loc(loc, iv, a_dx, center);

            // note that we don't include the singular part of psi
            // for the BHs - this is added at the output data stage
            // and when we calculate psi_reg in the rhs etc
            // as it already satisfies Laplacian(psi) = 0
            multigrid_vars_box(iv, c_psi_reg) =
                m_method_params.regularised_part_psi;
        }
    }
}

template <typename matter_t>
void CTT<matter_t>::initialise_constraint_vars(
    LevelData<FArrayBox> &a_constraint_vars, const RealVect &a_dx) const
{

    DataIterator dit = a_constraint_vars.dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
    {
        // These contain the vars in the boxes, set them all to zero
        FArrayBox &constraint_vars_box = a_constraint_vars[dit()];

        for (int comp = 0; comp < NUM_CONSTRAINT_VARS; comp++)
        {
            constraint_vars_box.setVal(0.0, comp);
        }
    }
}
