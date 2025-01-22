/* GRTresna
 * Copyright 2024 The GRTL Collaboration.
 * Please refer to LICENSE in GRTresna's root directory.
 */

#ifndef TENSOR_HPP_
#define TENSOR_HPP_

#include "DimensionDefinitions.hpp"

/// This class implements a Tensor with given rank, element data type, and
/// dimension.  By default the dimension is equal to spatial dims.
template <int rank, class data_t, int size = CH_SPACEDIM> class Tensor
{
    template <int, class, int> friend class Tensor;
    typedef typename Tensor<rank - 1, data_t, size>::arr_t arr_t[size];
    arr_t arr;

  public:
    Tensor() {}

    template <typename... T> Tensor(T... data) : arr{data...} {}

    operator arr_t &() { return arr; }

    operator const arr_t &() const { return arr; }
};

template <class data_t, int size> class Tensor<0, data_t, size>
{
    template <int, class, int> friend class Tensor;
    typedef data_t arr_t;
    arr_t arr;

  public:
    Tensor() {}

    Tensor(data_t val) : arr(val) {}

    operator arr_t &() { return arr; }

    operator const arr_t &() const { return arr; }
};

#endif /* TENSOR_HPP_ */
