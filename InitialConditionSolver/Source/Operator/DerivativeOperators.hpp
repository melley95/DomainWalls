#ifndef DERIVATIVEOPERATORS_HPP_
#define DERIVATIVEOPERATORS_HPP_

#include "DimensionDefinitions.hpp"
#include "FArrayBox.H"
#include "IntVect.H"
#include "Interval.H"
#include "MultigridUserVariables.hpp"
#include "REAL.H"
#include "RealVect.H"
#include "Tensor.hpp"
#include "UsingNamespace.H"

class DerivativeOperators
{
  public:
    DerivativeOperators(const RealVect &a_dx) : m_dx(a_dx){};

    void get_d1(Tensor<1, Real, SpaceDim> &d1, const IntVect &a_iv,
                const FArrayBox &a_vars_box, const int icomp);

    void get_d2(Tensor<2, Real, SpaceDim> &d2, const IntVect &a_iv,
                const FArrayBox &a_vars_box, const int icomp);

    void get_d1_vector(Tensor<2, Real, SpaceDim> &d1, const IntVect &a_iv,
                       const FArrayBox &a_vars_box, const Interval &a_interval);

    void get_d2_vector(Tensor<3, Real, SpaceDim> &d2, const IntVect &a_iv,
                       const FArrayBox &a_vars_box, const Interval &a_interval);

    void scalar_Laplacian(Real &laplacian, const IntVect &a_iv,
                          const FArrayBox &a_vars_box, const int a_comp);

    void vector_Laplacian(Tensor<1, Real, SpaceDim> &laplacian,
                          const IntVect &a_iv, const FArrayBox &a_vars_box,
                          const Interval &a_interval);

  private:
    const RealVect m_dx;
};

#endif /* DERIVATIVEOPERATORS_HPP_ */
