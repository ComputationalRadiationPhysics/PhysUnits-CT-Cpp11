#pragma once

#include "phys/units/quantity.hpp"
#include "phys/units/physical_constants.hpp"


namespace solver_helper
{
    using namespace phys::units;

    /* Yee */
    template< unsigned int dim >
    constexpr quantity< time_interval_d, double >
    optimal_step(
        quantity< length_d, double > const dx,
        quantity< length_d, double > const dy,
        quantity< length_d, double > const dz
    )
    {
        return sqrt(
            1.0 /
            (
                c * c *
                (
                    1.0 / ( dx * dx ) * (dim >= 1) +
                    1.0 / ( dy * dy ) * (dim >= 2) +
                    1.0 / ( dz * dz ) * (dim >= 3)
                )
            )
        );
    }
}
