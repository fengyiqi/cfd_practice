#ifndef UTILS_H
#define UTILS_H

#include "user_specification.h"
#include <stdexcept>

double ConvertConservativesToFlux(const double (&conservatives)[FI::EN()], const double pressure, const std::underlying_type<FI::EquationEum>::type equation) {
    if (equation == FI::Equations[EIndex(FI::EquationEum::Mass)]) 
        return conservatives[EIndex(FI::EquationEum::MomentumX)];
    else if (equation == FI::Equations[EIndex(FI::EquationEum::MomentumX)]) {
        return (conservatives[EIndex(FI::EquationEum::MomentumX)] * conservatives[EIndex(FI::EquationEum::MomentumX)]) / conservatives[EIndex(FI::EquationEum::Mass)] + pressure;
    }
    else if (equation == FI::Equations[EIndex(FI::EquationEum::Energy)]) {
        return conservatives[EIndex(FI::EquationEum::MomentumX)] / conservatives[EIndex(FI::EquationEum::Mass)] * (conservatives[EIndex(FI::EquationEum::Energy)] + pressure);
    }
    else
        throw std::invalid_argument("Invalid conservative state!");
}

#endif