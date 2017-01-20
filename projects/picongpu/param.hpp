#pragma once

#include "solver_helper.hpp"
#include "phys/units/quantity.hpp"

#include <string>

using namespace phys::units;
using namespace phys::units::literals;

constexpr unsigned int simDim = 3u;

constexpr quantity< length_d, double > dtestExt = 0.2e-6_m;

struct Grid
{
    static constexpr quantity< length_d, double > dx = 0.2e-6_m;
    static constexpr quantity< length_d, double > dy = 0.2e-6_m;
    static constexpr quantity< length_d, double > dz = 0.2e-6_m;
};

struct Iteration
{
    static constexpr quantity< time_interval_d, double > t =
        solver_helper::optimal_step< simDim >(
            Grid::dx, Grid::dy, Grid::dz
        );
};
//constexpr const quantity< time_interval_d, double > Iteration::t;

