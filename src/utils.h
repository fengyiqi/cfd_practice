#ifndef UTILS_H
#define UTILS_H

#include "user_specification.h"
#include <stdexcept>

namespace Utils {
    double ConvertConservativesToFlux(const double (&conservatives)[FI::EN()], const double pressure, const std::underlying_type<FI::EquationEnum>::type equation);
}
#endif