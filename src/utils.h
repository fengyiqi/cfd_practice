#ifndef UTILS_H
#define UTILS_H

#include "user_specification.h"
#include <stdexcept>

double ConvertConservativesToFlux(const double (&conservatives)[FI::CN()], const double pressure, const unsigned int& state) {
    if (state == ConservativePool::Mass) 
        return conservatives[ConservativePool::MomentumX];
    else if (state == ConservativePool::MomentumX) {
        return (conservatives[ConservativePool::MomentumX] * conservatives[ConservativePool::MomentumX]) / conservatives[ConservativePool::Mass] + pressure;
    }
    else if (state == ConservativePool::Energy) {
        return conservatives[ConservativePool::MomentumX] / conservatives[ConservativePool::Mass] * (conservatives[ConservativePool::Energy] + pressure);
    }
    else
        throw std::invalid_argument("Invalid conservative state!");
}

#endif