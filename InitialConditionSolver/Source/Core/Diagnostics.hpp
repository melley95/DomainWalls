#ifndef DIAGNOSTICS_HPP_
#define DIAGNOSTICS_HPP_

template <typename method_t, typename matter_t> class Diagnostics
{
  public:
    Diagnostics(method_t *a_method, matter_t *a_matter, Metric *a_metric,
                const Real a_G_Newton,
                const std::array<double, SpaceDim> a_center);

    void compute_constraint_terms(LevelData<FArrayBox> *a_multigrid_vars,
                                  LevelData<FArrayBox> *a_diagnostic_vars,
                                  LevelData<FArrayBox> *a_rhs,
                                  const RealVect &a_dx) const;

    void normalise_constraints(LevelData<FArrayBox> *a_multigrid_vars,
                               LevelData<FArrayBox> *a_diagnostic_vars,
                               LevelData<FArrayBox> *a_rhs,
                               const RealVect &a_dx, IntVect &nCells) const;


    void initialise_diagnostic_vars(LevelData<FArrayBox> &a_diagnostic_vars,
                                    const RealVect &a_dx) const;

  private:
    method_t *method;
    matter_t *matter;

    Metric *metric;

    const std::array<double, SpaceDim> center;

    const Real G_Newton;
};

#include "Diagnostics.impl.hpp"

#endif /* DIAGNOSTICS_HPP_ */