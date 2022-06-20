#ifndef UTILS_H
#define UTILS_H

#include "user_specification.h"
#include <stdexcept>
#include "stencils/stencil.h"

namespace Utils {
    double ConvertConservativesToFlux(const double (&conservatives)[FI::EN()], const double pressure, const std::underlying_type<FI::EquationEnum>::type equation) ;
    //     if (equation == FI::Equations[EIndex(FI::EquationEnum::Mass)]) 
    //         return conservatives[EIndex(FI::EquationEnum::MomentumX)];
    //     else if (equation == FI::Equations[EIndex(FI::EquationEnum::MomentumX)]) {
    //         return (conservatives[EIndex(FI::EquationEnum::MomentumX)] * conservatives[EIndex(FI::EquationEnum::MomentumX)]) / conservatives[EIndex(FI::EquationEnum::Mass)] + pressure;
    //     }
    //     else if (equation == FI::Equations[EIndex(FI::EquationEnum::Energy)]) {
    //         return conservatives[EIndex(FI::EquationEnum::MomentumX)] / conservatives[EIndex(FI::EquationEnum::Mass)] * (conservatives[EIndex(FI::EquationEnum::Energy)] + pressure);
    //     }
    //     else
    //         throw std::invalid_argument("Invalid conservative state!");
    // }
}

namespace StencilApplyUtils {
    template<typename StencilType>
    constexpr double StencilApply(const double (&array)[Stencil<StencilType>::StencilSize()]) {
        constexpr StencilType stencil = StencilType();
        return stencil.template Apply<Stencil<StencilType>>(array);
    }
}
#endif