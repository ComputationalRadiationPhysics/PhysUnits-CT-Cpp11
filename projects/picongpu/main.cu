/* general includes needed in param files
 * for lazy quantity definitions:
 *     constexpr quantity<length_d, T=PHYS_UNITS_REP_TYPE>
 */
//define PHYS_UNITS_REP_TYPE float
#include "phys/units/quantity.hpp"
// for printing
#include "phys/units/quantity_io.hpp"

/* constants */
#include "phys/units/physical_constants.hpp"

/* CT math */

/* used namespaces */
using namespace phys::units;
using namespace phys::units::literals;

/* define more literals for plasma physics */
namespace phys
{
namespace units
{
namespace literals
{
    // more literals
    QUANTITY_DEFINE_LITERALS( V, electric_potential_d )
    QUANTITY_DEFINE_LITERALS( T, magnetic_flux_density_d )
    //QUANTITY_DEFINE_LITERALS( K, thermodynamic_temperature_d )
    QUANTITY_DEFINE_SCALING_LITERALS( qe, electric_charge_d, -1.602176462e-19 )
    QUANTITY_DEFINE_SCALING_LITERALS( qp, electric_charge_d, 1.602176462e-19 )
    QUANTITY_DEFINE_SCALING_LITERALS( me, mass_d, 9.10938356e-31 )
    QUANTITY_DEFINE_SCALING_LITERALS( mp, mass_d, 1.672621898e-27 )
    QUANTITY_DEFINE_SCALING_LITERALS( eV, energy_d, 1.60217733e-19 )
}
}
}

/* input helper file(s) */
#include "solver_helper.hpp"

/* param file(s) */
#include "param.hpp"

/* PIConGPU unit system */
namespace picRatioSI
{
    constexpr quantity< length_d > length{
        // SI::DELTA_T * TYPICAL_SPEED (c)
        Iteration::t * c
    };
    constexpr quantity< mass_d > mass{
        // SI::BASE_MASS (* typical weighting)
        detail::magnitude_tag, 9.10938356e-31
    };
    constexpr quantity< time_interval_d > time_interval{
        // SI::DELTA_T
        Iteration::t
    };
    constexpr quantity< electric_current_d > electric_current{
        // SI::BASE_CHARGE (* typical weighting) / UNIT_TIME
        e / time_interval
    };
    constexpr quantity< thermodynamic_temperature_d > thermodynamic_temperature{
        detail::magnitude_tag, 1.0
    };
    constexpr quantity< amount_of_substance_d > amount_of_substance{
        detail::magnitude_tag, 1.0
    };
    constexpr quantity< luminous_intensity_d > luminous_intensity{
        detail::magnitude_tag, 1.0
    };
}

/* to PIConGPU unit system */
template< typename T_Dimensions >
struct ToPIC;

template<>
template< int D1, int D2, int D3, int D4, int D5, int D6, int D7 >
struct ToPIC< dimensions< D1, D2, D3, D4, D5, D6, D7 > >
{
    using dim = dimensions< D1, D2, D3, D4, D5, D6, D7 >;
    
    template<
        typename DX,
        typename X
    >
    constexpr auto
    operator()( quantity< DX, X > const & x ) const
    -> X //detail::Quotient<dim, DX, dim, X>
    {
        /*static_assert(
            std::is_same<DX, dim>::value,
            "Dimension mismatch in variable ... macro magic..."
        );*/
        return x /
            (
                nth_power< D1 >( picRatioSI::length ) *
                nth_power< D2 >( picRatioSI::mass ) *
                nth_power< D3 >( picRatioSI::time_interval ) *
                nth_power< D4 >( picRatioSI::electric_current ) *
                nth_power< D5 >( picRatioSI::thermodynamic_temperature ) *
                nth_power< D6 >( picRatioSI::amount_of_substance ) *
                nth_power< D7 >( picRatioSI::luminous_intensity )
            );
    }
};

#include <iostream>
//include <cuda_runtime.h>


/* a host funtion using the input */
void
printStuff()
{
    using namespace phys::units::io;
    
    constexpr auto ddd = Iteration::t;
    std::cout << ddd << std::endl;

    constexpr auto two_e = e + e + e - 1.0_qe;
    
    std::cout << ToPIC< electric_charge_d >()( two_e ) << std::endl;
    // must fail:
    //std::cout << ToPIC< length_d >()( two_e ) << std::endl;
}
/* a kernel using the input */
/*
__device__ void
dev_foo()
{
    constexpr auto ddd = B::t;
    printf("%d", ddd.to( second ) );
}*/


int main()
{
    printStuff();
    
    return 0;
}

