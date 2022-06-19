#ifndef STIFFENED_GAS_H
#define STIFFENED_GAS_H

#include "equation_of_state.h"
#include <cmath>

class StiffendGas : public EquationOfState {
public:
    StiffendGas() = delete;
    explicit StiffendGas(const double gamma) : EquationOfState(gamma) {}
    explicit StiffendGas(const double gamma, const double background_pressure) : EquationOfState(gamma, background_pressure) {}
    inline double ComputePressure(const double rho, const double e) {
        return (gamma_ - 1) * rho * e - gamma_ * background_pressure_;
    }
    inline double ComputeSpeedOfSound(const double rho, const double pressure) {
        return std::sqrt(gamma_ * (pressure + background_pressure_) / rho);
    }
};

#endif