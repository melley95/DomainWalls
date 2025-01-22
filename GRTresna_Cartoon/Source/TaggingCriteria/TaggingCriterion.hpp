/* GRTresna
 * Copyright 2024 The GRTL Collaboration.
 * Please refer to LICENSE in GRTresna's root directory.
 */

#ifndef TAGGINGCRITERION_HPP_
#define TAGGINGCRITERION_HPP_

class TaggingCriterion
{
  public:
    TaggingCriterion(){};

    virtual void set_regrid_condition(LevelData<FArrayBox> &a_condition,
                                      LevelData<FArrayBox> &a_multigrid_vars,
                                      const RealVect &a_dx,
                                      const std::array<double, SpaceDim> center,
                                      Real regrid_radius) = 0;

    virtual ~TaggingCriterion() = default;

  private:
};

#endif /* TAGGINGCRITERION_HPP_ */