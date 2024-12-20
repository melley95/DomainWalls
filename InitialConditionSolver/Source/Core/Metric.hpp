#ifndef METRIC_HPP_
#define METRIC_HPP_

#include "GRParmParse.hpp"
#include "REAL.H"
#include "RealVect.H"
#include "TensorAlgebra.hpp"

#include "DerivativeOperators.hpp"
#include "FArrayBox.H"
#include "Interval.H"
#include "UsingNamespace.H"

class Metric
{
  public:
    struct params_t
    {
        Real bh1_bare_mass;
        Real bh2_bare_mass;
        RealVect bh1_spin;
        RealVect bh2_spin;
        RealVect bh1_momentum;
        RealVect bh2_momentum;
        RealVect bh1_offset;
        RealVect bh2_offset;
        bool method_compact;
        std::array<double, SpaceDim> center;
    };

    static void read_params(GRParmParse &pp, params_t &a_metric_params);

    explicit Metric(params_t a_metric_params)
        : m_metric_params(a_metric_params)
    {
    }

    void get_bh_coords(Real &bh_radius, RealVect &loc_bh, const RealVect &loc,
                       const RealVect &bh_offset);

    Real compute_bowenyork_psi(const RealVect &loc);

    void compute_bowenyork_Aij(Tensor<2, Real> &Aij, // const IntVect &iv,
                               const RealVect &loc);

    void compute_ctt_Aij(Tensor<2, Real> &Aij,
                         const FArrayBox &multigrid_vars_box, const IntVect &iv,
                         const RealVect &a_dx, const RealVect &loc) const;

    void set_Aww_reg(Real &Aww, const FArrayBox &multigrid_vars_box,
                    const IntVect &iv, const RealVect &a_dx, const RealVect &loc);
                    
    params_t m_metric_params;
};

#endif /* METRIC_HPP_ */
