#ifndef GRCHOMBOUSERVARIABLES_HPP
#define GRCHOMBOUSERVARIABLES_HPP

#include "ArrayTools.hpp"

// assign an enum to each variable
enum
{
    c_chi,

    c_h11,
    c_h12,
#if CH_SPACEDIM == 3
    c_h13,
#endif
    c_h22,
#if CH_SPACEDIM == 3
    c_h23,
    c_h33,
#elif CH_SPACEDIM == 2
    c_hww,
#endif

    c_K,

    c_A11,
    c_A12,
#if CH_SPACEDIM == 3
    c_A13,
#endif
    c_A22,
#if CH_SPACEDIM == 3
    c_A23,
    c_A33,
#elif CH_SPACEDIM == 2
    c_Aww,
#endif

    c_Theta,

    c_Gamma1,
    c_Gamma2,
#if CH_SPACEDIM == 3
    c_Gamma3,
#endif

    c_lapse,

    c_shift1,
    c_shift2,
#if CH_SPACEDIM == 3
    c_shift3,
#endif

    c_B1,
    c_B2,
#if CH_SPACEDIM == 3
    c_B3,
#endif

    c_phi, // matter field added
    c_Pi,  //(minus) conjugate momentum

    NUM_GRCHOMBO_VARS
};

namespace GRChomboUserVariables
{
static constexpr char const *variable_names[NUM_GRCHOMBO_VARS] = {
    "chi",

    D_DECL("h11", "h12", "h13"),
    "h22",
#if CH_SPACEDIM == 3
    "h23",
    "h33",
#elif CH_SPACEDIM == 2
    "hww",
#endif

    "K",

    D_DECL("A11", "A12", "A13"),
    "A22",
#if CH_SPACEDIM == 3
    "A23",
    "A33",
#elif CH_SPACEDIM == 2
    "Aww",
#endif

    "Theta",

    D_DECL("Gamma1", "Gamma2", "Gamma3"),

    "lapse",

    D_DECL("shift1", "shift2", "shift3"),

    D_DECL("B1", "B2", "B3"),

    "phi",
    "Pi"};

static constexpr std::array<int, NUM_GRCHOMBO_VARS> const vars_parity = {
    0,
    D_DECL(0, 3, 6),
    0,
#if CH_SPACEDIM == 3
    5,
    0,
#elif CH_SPACEDIM == 2
    0,
#endif
    0,
    D_DECL(0, 3, 6),
    0,
#if CH_SPACEDIM == 3
    5,
    0,
#elif CH_SPACEDIM == 2
    0,
#endif
    0,
    D_DECL(1, 2, 3),
    0,
    D_DECL(1, 2, 3),
    D_DECL(1, 2, 3),
    0,
    0};

} // namespace GRChomboUserVariables

#endif /* GRCHOMBOUSERVARIABLES_HPP */