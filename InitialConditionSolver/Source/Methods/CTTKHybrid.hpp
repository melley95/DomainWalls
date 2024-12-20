#ifndef CTTKHybrid_HPP_
#define CTTKHybrid_HPP_

#include "BiCGStabSolver.H"
#include "Metric.hpp"
#include "GRParmParse.hpp"
#include "Grids.hpp"
#include "MultilevelLinearOp.H"
#include "Tensor.hpp"

template <typename matter_t> class CTTKHybrid
{
  public:
    struct params_t;

    CTTKHybrid() {}

    CTTKHybrid(params_t a_method_params, matter_t *a_matter,
               Metric *a_metric, int a_numLevels,
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
                            RefCountedPtr<LevelData<FArrayBox>> a_cCoef,
                            const RealVect &a_dx);

    params_t m_method_params;
    Metric::params_t m_metric_params;

    matter_t *matter;
    Metric *metric;

    ~CTTKHybrid() {}

  private:
    int numLevels;

    const std::array<double, SpaceDim> center;

    const Real G_Newton;
};

#include "CTTKHybrid.impl.hpp"

#endif /* CTTKHybrid_HPP_ */