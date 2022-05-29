#include "fixed_value.h"
#include <iostream>

FixedValueBoundaryCondition::FixedValueBoundaryCondition(const double (&left_side)[FI::PN()], const double (&right_side)[FI::PN()]) {
    for (unsigned int pr = 0; pr < FI::PN(); pr++) {
        left_primitives_[pr]  = left_side[pr];
        right_primitives_[pr] = right_side[pr];
    }
    // temporary solution
    double gamma = 1.4;
    // mass
    left_conservatives_[ConservativePool::Mass] = left_side[PrimeStatePool::Density];
    right_conservatives_[ConservativePool::Mass] = right_side[PrimeStatePool::Density];
    // momentum
    left_conservatives_[ConservativePool::MomentumX] = left_side[PrimeStatePool::Density] * left_side[PrimeStatePool::VelocityX];
    right_conservatives_[ConservativePool::MomentumX] = right_side[PrimeStatePool::Density] * right_side[PrimeStatePool::VelocityX];
    // energy
    left_conservatives_[ConservativePool::Energy] = left_side[PrimeStatePool::Density] * (left_side[PrimeStatePool::Pressure] / (left_side[PrimeStatePool::Density] * (gamma - 1.0)) +
                                                        0.5 * left_side[PrimeStatePool::VelocityX] * left_side[PrimeStatePool::VelocityX]);
    right_conservatives_[ConservativePool::Energy] = right_side[PrimeStatePool::Density] * (right_side[PrimeStatePool::Pressure] / (right_side[PrimeStatePool::Density] * (gamma - 1.0)) +
                                                        0.5 * right_side[PrimeStatePool::VelocityX] * right_side[PrimeStatePool::VelocityX]);


}

void FixedValueBoundaryCondition::Apply(double (&buffer)[FI::PN()][GI::TCX()], unsigned int state) {
    if (state == States::Primitives) {
        for (unsigned int i = 0; i < FI::PN(); i++) {
            // update left side
            for (unsigned int j = 0; j < GI::GCX(); j++)
                buffer[i][j] = left_primitives_[i];
            // update right side
            for (unsigned int j = GI::FHHX(); j < GI::TCX(); j++)
                buffer[i][j] = right_primitives_[i];
        }
    }
    else {  // conservative states
        for (unsigned int i = 0; i < FI::CN(); i++) {
            // update left side
            for (unsigned int j = 0; j < GI::GCX(); j++)
                buffer[i][j] = left_conservatives_[i];
            // update right side
            for (unsigned int j = GI::FHHX(); j < GI::TCX(); j++)
                buffer[i][j] = right_conservatives_[i];
        }
    }
}