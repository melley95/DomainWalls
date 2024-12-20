#ifndef CTTKHybrid_HPP_
#error "This file should only be included through CTTKHybrid.hpp"
#endif

#include "DimensionDefinitions.hpp"
#include "GRParmParse.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"

template <typename matter_t> struct CTTKHybrid<matter_t>::params_t
{
    int sign_of_K;
    bool method_compact;
    Real psi_reg;

    bool deactivate_zero_mode;
};

template <typename matter_t>
CTTKHybrid<matter_t>::CTTKHybrid(params_t a_method_params, matter_t *a_matter,
                                 Metric *a_metric, int a_numLevels,
                                 const std::array<double, SpaceDim> a_center,
                                 Real a_G_Newton)
    : m_method_params(a_method_params), matter(a_matter), G_Newton(a_G_Newton),
      numLevels(a_numLevels), metric(a_metric), center(a_center)
{
}

template <typename matter_t>
void CTTKHybrid<matter_t>::read_params(GRParmParse &pp,
                                       params_t &a_method_params)
{
    pp.get("sign_of_K", a_method_params.sign_of_K);

    a_method_params.method_compact = false;
    pp.get("method_compact", a_method_params.method_compact);

    pp.get("psi_reg", a_method_params.psi_reg);

    a_method_params.deactivate_zero_mode = false;
    pp.get("deactivate_zero_mode", a_method_params.deactivate_zero_mode);
}

template <typename matter_t>
void CTTKHybrid<matter_t>::solve_analytic(
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
            Real psi_bh = metric->compute_bowenyork_psi(loc);
            Real psi_0 = psi_reg + psi_bh;
            Real laplacian_psi_reg;
            derivs.scalar_Laplacian(laplacian_psi_reg, iv, multigrid_vars_box,
                                    c_psi_reg);

            // Assign values of Aij
            Tensor<2, Real> Aij_reg;
            metric->compute_ctt_Aij(Aij_reg, multigrid_vars_box, iv, a_dx, loc);
            Tensor<2, Real> Aij_bh;
            metric->compute_bowenyork_Aij(Aij_bh, loc);
            // This is \bar  A_ij \bar A^ij
            Real A2_0 = 0.0;
            FOR2(i, j)
            {
                A2_0 += (Aij_reg[i][j] + Aij_bh[i][j]) *
                        (Aij_reg[i][j] + Aij_bh[i][j]);
            }
            Real Aww_reg; // Cartoon term
            if (SpaceDim == 2)
            {
                metric->set_Aww_reg(Aww_reg, multigrid_vars_box, iv, a_dx, loc);
                A2_0 += Aww_reg * Aww_reg;
            }

            // Compute emtensor components
            const auto emtensor =
                matter->compute_emtensor(iv, a_dx, multigrid_vars_box);

            // Set value for K
            Real K_0_squared = 24.0 * M_PI * G_Newton * emtensor.rho;

            // be careful if at a point K = 0, may have discontinuity
            multigrid_vars_box(iv, c_K_0) =
                m_method_params.sign_of_K * sqrt(K_0_squared);

            // set values for \bar Aij_0
            multigrid_vars_box(iv, c_A11_0) = Aij_reg[0][0] + Aij_bh[0][0];
            multigrid_vars_box(iv, c_A22_0) = Aij_reg[1][1] + Aij_bh[1][1];
            multigrid_vars_box(iv, c_A12_0) = Aij_reg[0][1] + Aij_bh[0][1];
#if CH_SPACEDIM == 3
            multigrid_vars_box(iv, c_A13_0) = Aij_reg[0][2] + Aij_bh[0][2];
            multigrid_vars_box(iv, c_A23_0) = Aij_reg[1][2] + Aij_bh[1][2];
            multigrid_vars_box(iv, c_A33_0) = Aij_reg[2][2] + Aij_bh[2][2];
#elif CH_SPACEDIM == 2
            multigrid_vars_box(iv, c_Aww_0) = Aww_reg;
#endif
        }
    }
}

template <typename matter_t>
void CTTKHybrid<matter_t>::set_elliptic_terms(
    LevelData<FArrayBox> *a_multigrid_vars, LevelData<FArrayBox> *a_rhs,
    RefCountedPtr<LevelData<FArrayBox>> a_aCoef,
    RefCountedPtr<LevelData<FArrayBox>> a_bCoef,
    RefCountedPtr<LevelData<FArrayBox>> a_cCoef, const RealVect &a_dx)
{
    DerivativeOperators derivs(a_dx);
    DataIterator dit = a_rhs->dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
    {
        FArrayBox &multigrid_vars_box = (*a_multigrid_vars)[dit()];
        FArrayBox &rhs_box = (*a_rhs)[dit()];
        FArrayBox &aCoef_box = (*a_aCoef)[dit()];
        FArrayBox &bCoef_box = (*a_bCoef)[dit()];
        FArrayBox &cCoef_box = (*a_cCoef)[dit()];
        // JCAurre: Initialise rhs=0, aCoef=0 and bCoef=1 for all constraint
        // variables
        for (int comp = 0; comp < NUM_CONSTRAINT_VARS; comp++)
        {
            rhs_box.setVal(0.0, comp);
            aCoef_box.setVal(0.0, comp);
            bCoef_box.setVal(1.0, comp);
            cCoef_box.setVal(0.0, comp);

            // // this prevents small amounts of noise in the sources
            // // activating the zero modes - (Garfinkle trick)
            // if (m_method_params.deactivate_zero_mode)
            // {
            //     // Seems to work best to set this relative to the tolerance
            //     aCoef_box.setVal(-1e4 * tolerance, iconstraint);
            // }
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
            Real psi_bh = metric->compute_bowenyork_psi(loc);
            Real psi_0 = psi_reg + psi_bh;
            Real laplacian_psi_reg;
            derivs.scalar_Laplacian(laplacian_psi_reg, iv, multigrid_vars_box,
                                    c_psi_reg);

            // Get values of Aij
            Tensor<2, Real> Aij_reg;
            metric->compute_ctt_Aij(Aij_reg, multigrid_vars_box, iv, a_dx, loc);
            Tensor<2, Real> Aij_bh;
            metric->compute_bowenyork_Aij(Aij_bh, loc);
            // This is \bar  A_ij \bar A^ij
            Real A2_0 = 0.0;
            FOR2(i, j)
            {
                A2_0 += (Aij_reg[i][j] + Aij_bh[i][j]) *
                        (Aij_reg[i][j] + Aij_bh[i][j]);
            }
            Real Aww_reg; // Cartoon term
            if (SpaceDim == 2)
            {
                metric->set_Aww_reg(Aww_reg, multigrid_vars_box, iv, a_dx, loc);
                A2_0 += Aww_reg * Aww_reg;
            }

            // Compute emtensor components
            const auto emtensor =
                matter->compute_emtensor(iv, a_dx, multigrid_vars_box);

            Tensor<1, Real, SpaceDim> d1_K;
            derivs.get_d1(d1_K, iv, multigrid_vars_box, c_K_0);

            // rhs terms, K is set to cancel matter terms only
            rhs_box(iv, c_psi) =
                -0.125 * A2_0 * pow(psi_0, -7.0) - laplacian_psi_reg;

            // Get d_i V_i and laplacians
#if CH_SPACEDIM == 3
            Tensor<2, Real, SpaceDim> di_Vi;
            derivs.get_d1_vector(di_Vi, iv, multigrid_vars_box,
                                 Interval(c_V1_0, c_V3_0));

            Tensor<1, Real, SpaceDim> laplacian_V;
            derivs.vector_Laplacian(laplacian_V, iv, multigrid_vars_box,
                                    Interval(c_V1_0, c_V3_0));
#endif
#if CH_SPACEDIM == 2
            Tensor<2, Real, SpaceDim> di_Vi;
            derivs.get_d1_vector(di_Vi, iv, multigrid_vars_box,
                                 Interval(c_V1_0, c_V2_0));

            Tensor<1, Real, SpaceDim> laplacian_V;
            derivs.vector_Laplacian(laplacian_V, iv, multigrid_vars_box,
                                    Interval(c_V1_0, c_V2_0));
#endif
            Real laplacian_U;
            derivs.scalar_Laplacian(laplacian_U, iv, multigrid_vars_box, c_U_0);

            // now set the values in the box
            rhs_box(iv, c_V1) =
                pow(psi_0, 6.0) *
                (2.0 / 3.0 * d1_K[0] + 8.0 * M_PI * G_Newton * emtensor.Si[0]);
            rhs_box(iv, c_V2) =
                pow(psi_0, 6.0) *
                (2.0 / 3.0 * d1_K[1] + 8.0 * M_PI * G_Newton * emtensor.Si[1]);
#if CH_SPACEDIM == 3
            rhs_box(iv, c_V3) =
                pow(psi_0, 6.0) *
                (2.0 / 3.0 * d1_K[2] + 8.0 * M_PI * G_Newton * emtensor.Si[2]);
#endif
            // Cartoon terms
            if (SpaceDim == 2)
            {
                RealVect cartoon_loc;
                Grids::get_loc_cartoon(cartoon_loc, iv, a_dx);
                int cartoon_idx = 1;
                // rhs_box(iv, c_V1) +=                 //THIS NEEDS TO BE
                // ADAPTED
                //     -(d1_V1[cartoon_idx] / cartoon_loc[cartoon_idx]);
                // rhs_box(iv, c_V2) += -(
                //     d1_V2[cartoon_idx] / cartoon_loc[cartoon_idx] -
                //     multigrid_vars_box(iv, c_V2_0) /
                //         (cartoon_loc[cartoon_idx] *
                //         cartoon_loc[cartoon_idx]));
            }
            // Periodic: Use ansatz B.3 in B&S (p547) JCA TODO: We are not using
            // this U when constructing Aij Non-periodic: Compact ansatz B.7 in
            // B&S (p547)
            if (!m_method_params.method_compact)
            {
                rhs_box(iv, c_U) = -0.25 * (di_Vi[0][0] + di_Vi[1][1]);
#if CH_SPACEDIM == 3
                rhs_box(iv, c_U) += -0.25 * di_Vi[2][2];
#endif
                // Cartoon terms
                if (SpaceDim == 2)
                {
                    RealVect cartoon_loc;
                    Grids::get_loc_cartoon(cartoon_loc, iv, a_dx);
                    int cartoon_idx = 1;
                    // Tensor<1, Real, SpaceDim> d1_U = // THIS NEEDS TO BE
                    // ADAPTED
                    //     derivs.get_d1_vector(iv, multigrid_vars_box, a_dx,
                    //     c_U_0);
                    // rhs_box(iv, c_U) +=
                    //     -d1_U[cartoon_idx] / cartoon_loc[cartoon_idx] -
                    //     .25 * multigrid_vars_box(iv, c_V2_0) /
                    //         cartoon_loc[cartoon_idx];
                }
            }
            else
            {
                rhs_box(iv, c_U) = 0.0;
                FOR1(i)
                {
                    rhs_box(iv, c_U) +=
                        -loc[i] * (pow(psi_0, 6.0) *
                                   (2.0 / 3.0 * d1_K[i] +
                                    8.0 * M_PI * G_Newton * emtensor.Si[i]));
                }
            }

            if (m_method_params.deactivate_zero_mode)
            {
                rhs_box(iv, c_V1) += -laplacian_V[0];
                rhs_box(iv, c_V2) += -laplacian_V[1];
#if CH_SPACEDIM == 3
                rhs_box(iv, c_V3) += -laplacian_V[2];
#endif
                rhs_box(iv, c_U) += -laplacian_U;
            }

            // add the aCoef term
            aCoef_box(iv, c_psi) += -0.875 * A2_0 * pow(psi_0, -8.0);

            Grids::get_loc_cartoon(loc, iv, a_dx);
            Real yy = loc[1];
            // Cartoon term
            aCoef_box(iv, c_V2) += -1.0 / pow(yy, 2.0);

            cCoef_box(iv, c_V1) += 1.0 / yy;
            cCoef_box(iv, c_V2) += 1.0 / yy;
            cCoef_box(iv, c_U) += 1.0 / yy;
        }
    }
}

template <typename matter_t>
void CTTKHybrid<matter_t>::initialise_method_vars(
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
            multigrid_vars_box(iv, c_psi_reg) = m_method_params.psi_reg;
        }
    }
}

template <typename matter_t>
void CTTKHybrid<matter_t>::initialise_constraint_vars(
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
