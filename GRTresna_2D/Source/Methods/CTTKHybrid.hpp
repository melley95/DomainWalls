/* GRTresna
 * Copyright 2024 The GRTL Collaboration.
 * Please refer to LICENSE in GRTresna's root directory.
 */

#ifndef CTTKHybrid_HPP_
#define CTTKHybrid_HPP_

#include "BiCGStabSolver.H"
#include "GRParmParse.hpp"
#include "Grids.hpp"
#include "MultilevelLinearOp.H"
#include "PsiAndAijFunctions.hpp"
#include "Tensor.hpp"

template <typename matter_t> class CTTKHybrid
{
  public:
    struct params_t;

    CTTKHybrid() {}

    CTTKHybrid(params_t a_method_params, matter_t *a_matter,
               PsiAndAijFunctions *a_psi_and_Aij_functions, int a_numLevels,
               const std::array<double, SpaceDim> a_center, Real a_G_Newton);

    static void read_params(GRParmParse &pp, params_t &a_method_params);

    void initialise_method_vars(LevelData<FArrayBox> &a_multigrid_vars,
                                const RealVect &a_dx) const;

    void initialise_constraint_vars(LevelData<FArrayBox> &a_constraint_vars,
                                    const RealVect &a_dx) const;

    void solve_analytic(LevelData<FArrayBox> *multigrid_vars,
                        LevelData<FArrayBox> *rhs, const RealVect &a_dx);

    void set_elliptic_terms(LevelData<FArrayBox> *a_multigrid_vars,
                            LevelData<FArrayBox> *a_rhs,
                            RefCountedPtr<LevelData<FArrayBox>> a_aCoef,
                            RefCountedPtr<LevelData<FArrayBox>> a_bCoef,
                            const RealVect &a_dx);

    params_t m_method_params;
    PsiAndAijFunctions::params_t m_psi_and_Aij_params;

    matter_t *matter;
    PsiAndAijFunctions *psi_and_Aij_functions;

    ~CTTKHybrid() {}

  private:
    int numLevels;

    const std::array<double, SpaceDim> center;

    const Real G_Newton;
};

#include "CTTKHybrid.impl.hpp"

#endif /* CTTKHybrid_HPP_ */