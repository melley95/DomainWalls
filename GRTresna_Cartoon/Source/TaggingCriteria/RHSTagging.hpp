/* GRTresna
 * Copyright 2024 The GRTL Collaboration.
 * Please refer to LICENSE in GRTresna's root directory.
 */

#ifndef RHSTAGGING_HPP_
#define RHSTAGGING_HPP_

#include "DerivativeOperators.hpp"
#include "REAL.H"
#include "TaggingCriterion.hpp"

template <typename method_t, typename matter_t>
class RHSTagging : public TaggingCriterion
{
  public:
    RHSTagging(method_t *a_method, matter_t *a_matter, Real a_G_Newton)
        : TaggingCriterion(), method(a_method), matter(a_matter),
          G_Newton(a_G_Newton)
    {
    }

    ~RHSTagging() {}

    void set_regrid_condition(LevelData<FArrayBox> &a_condition,
                              LevelData<FArrayBox> &a_multigrid_vars,
                              const RealVect &a_dx,
                              const std::array<double, SpaceDim> center,
                              Real regrid_x, Real regrid_y);

  private:
    method_t const *method;
    matter_t const *matter;

    const Real G_Newton;
};

template <typename method_t, typename matter_t>
void RHSTagging<method_t, matter_t>::set_regrid_condition(
    LevelData<FArrayBox> &a_condition, LevelData<FArrayBox> &a_multigrid_vars,
    const RealVect &a_dx, const std::array<double, SpaceDim> center,
    Real regrid_x, Real regrid_y)
{
    DerivativeOperators derivs(a_dx);
    CH_assert(a_multigrid_vars.nComp() == NUM_MULTIGRID_VARS);

    method->initialise_method_vars(a_multigrid_vars, a_dx);
    matter->initialise_matter_vars(a_multigrid_vars, a_dx);

    DataIterator dit = a_condition.dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
    {
        FArrayBox &multigrid_vars_box = a_multigrid_vars[dit()];
        FArrayBox &condition_box = a_condition[dit()];
        condition_box.setVal(0.0, 0);

        Box unghosted_box = condition_box.box();
        BoxIterator bit(unghosted_box);
        for (bit.begin(); bit.ok(); ++bit)
        {
            // work out location on the grid
            IntVect iv = bit();
            RealVect loc;
            Grids::get_loc(loc, iv, a_dx, center);

            int cartoon_idx = 1;
            Real yy = loc[cartoon_idx];

            // Calculate the actual value of psi including BH part
            Real psi_reg = multigrid_vars_box(iv, c_psi_reg);
            
            Real psi_0 = psi_reg;

           

            Tensor<1, Real, SpaceDim> d1_psi_reg;
            derivs.get_d1(d1_psi_reg, iv, multigrid_vars_box, c_psi_reg);

            Real laplacian_psi_reg;
            derivs.scalar_Laplacian(laplacian_psi_reg, iv, multigrid_vars_box,
                                    c_psi_reg);

            // Get values of Aij
            Tensor<2, Real> Aij_reg;
            method->psi_and_Aij_functions->compute_ctt_Aij(
                Aij_reg, multigrid_vars_box, iv, a_dx, loc);
   //         Tensor<2, Real> Aij_bh;
     //       method->psi_and_Aij_functions->compute_bowenyork_Aij(Aij_bh, loc);
            // This is \bar  A_ij \bar A^ij
            Real A2_0 = 0.0;
            FOR2(i, j)
            {
                A2_0 += (Aij_reg[i][j]  *
                        Aij_reg[i][j] );
            }

            Real Aww_reg; // Cartoon term

            method->psi_and_Aij_functions->set_Aww_reg(Aww_reg, multigrid_vars_box, iv, a_dx, loc);
            A2_0 += Aww_reg * Aww_reg;

            // Compute emtensor components
            const auto emtensor =
                matter->compute_emtensor(iv, a_dx, multigrid_vars_box);


            
           

            Tensor<1, Real, SpaceDim> d1_phi;
            Real mod_d1_phi;
            derivs.get_d1(d1_phi, iv, multigrid_vars_box, c_phi_0);

            FOR(idir)
        {
            mod_d1_phi += d1_phi[idir] * d1_phi[idir];
          
        }



            if (regrid_x > 0 && regrid_y > 0)
            {
             //   Real rr = sqrt(loc[0] * loc[0] +loc[1] * loc[1]);
                if (loc[0] < regrid_x && loc[1] < regrid_y)
                {
                    condition_box(iv, 0) = 1.0;
                }
            }
            else
            {
                // the condition is similar to the rhs but we take abs
                // value of the contributions and add in effect of psi_0 via log
                condition_box(iv, 0) =             //a_dx[0] * (sqrt(mod_d1_phi));


                    2.0 * M_PI * G_Newton * emtensor.rho + abs(0.125 * A2_0) +
                    log(psi_0) + laplacian_psi_reg +
                    8.0 * M_PI * G_Newton *
                        (abs(emtensor.Si[0]) + abs(emtensor.Si[1]))  
                    + d1_psi_reg[cartoon_idx] / yy;
            }
        }
    }
}

#endif /* RHSTAGGING_HPP_*/