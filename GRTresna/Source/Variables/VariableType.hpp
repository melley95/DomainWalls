/* GRTresna
 * Copyright 2024 The GRTL Collaboration.
 * Please refer to LICENSE in GRTresna's root directory.
 */

#ifndef VARIABLETYPE_HPP
#define VARIABLETYPE_HPP

// enum for multigrid, grchombo or constraint vars
enum class VariableType
{
    multigrid, // The metric plus the matter vars
    grchombo,  // The output vars for GRChombo
    constraint // The vars that go into the linear solver steps
};

#endif /* VARIABLETYPE_HPP */