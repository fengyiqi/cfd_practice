#include "utils.h"

double Utils::ConvertConservativesToFlux(const double (&conservatives)[FI::EN()], const double pressure, const std::underlying_type<FI::EquationEnum>::type equation) {
    if (equation == FI::Equations[EIndex(FI::EquationEnum::Mass)]) 
        return conservatives[EIndex(FI::EquationEnum::MomentumX)];
    else if (equation == FI::Equations[EIndex(FI::EquationEnum::MomentumX)]) {
        return (conservatives[EIndex(FI::EquationEnum::MomentumX)] * conservatives[EIndex(FI::EquationEnum::MomentumX)]) / conservatives[EIndex(FI::EquationEnum::Mass)] + pressure;
    }
    else if (equation == FI::Equations[EIndex(FI::EquationEnum::Energy)]) {
        return conservatives[EIndex(FI::EquationEnum::MomentumX)] / conservatives[EIndex(FI::EquationEnum::Mass)] * (conservatives[EIndex(FI::EquationEnum::Energy)] + pressure);
    }
    else
        throw std::invalid_argument("Invalid conservative state!");
}