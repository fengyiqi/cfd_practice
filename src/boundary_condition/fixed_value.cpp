#include "fixed_value.h"

FixedValueBoundaryCondition::FixedValueBoundaryCondition(const double (&left_side)[FI::PN()], const double (&right_side)[FI::PN()]) {
    for (auto& p : FI::PrimeStates) {
        left_primitives_[p]  = left_side[p];
        right_primitives_[p] = right_side[p];
    }
    // temporary solution
    double gamma = 1.4;
    // mass
    left_conservatives_[EIndex(FI::EquationEnum::Mass)] = left_side[PIndex(FI::PrimeStateEnum::Density)];
    right_conservatives_[EIndex(FI::EquationEnum::Mass)] = right_side[PIndex(FI::PrimeStateEnum::Density)];
    // momentum
    left_conservatives_[EIndex(FI::EquationEnum::MomentumX)] = left_side[PIndex(FI::PrimeStateEnum::Density)] * left_side[PIndex(FI::PrimeStateEnum::VelocityX)];
    right_conservatives_[EIndex(FI::EquationEnum::MomentumX)] = right_side[PIndex(FI::PrimeStateEnum::Density)] * right_side[PIndex(FI::PrimeStateEnum::VelocityX)];
    // energy
    left_conservatives_[EIndex(FI::EquationEnum::Energy)] = left_side[PIndex(FI::PrimeStateEnum::Pressure)] / (gamma - 1.0) // rho*e
                                                           + 0.5 * left_side[PIndex(FI::PrimeStateEnum::Density)] * left_side[PIndex(FI::PrimeStateEnum::VelocityX)] * left_side[PIndex(FI::PrimeStateEnum::VelocityX)]; // rho*uu/2
    right_conservatives_[EIndex(FI::EquationEnum::Energy)] = right_side[PIndex(FI::PrimeStateEnum::Pressure)] / (gamma - 1.0) // rho*e
                                                           + 0.5 * right_side[PIndex(FI::PrimeStateEnum::Density)] * right_side[PIndex(FI::PrimeStateEnum::VelocityX)] * right_side[PIndex(FI::PrimeStateEnum::VelocityX)]; // rho*uu/2
}

void FixedValueBoundaryCondition::Apply(double (&buffer)[FI::PN()][GI::TCX()], FI::States state) {
    if (state == FI::States::Primitives) {
        for (auto& p : FI::PrimeStates) {
            // update left side
            for (unsigned int i = 0; i < GI::GCX(); i++)
                buffer[p][i] = left_primitives_[p];
            // update right side
            for (unsigned int i = GI::FRGX(); i < GI::TCX(); i++)
                buffer[p][i] = right_primitives_[p];
        }
    }
    else {  // conservative states
        for (auto& e : FI::Equations) {
            // update left side
            for (unsigned int i = 0; i < GI::GCX(); i++)
                buffer[e][i] = left_conservatives_[e];
            // update right side
            for (unsigned int i = GI::FRGX(); i < GI::TCX(); i++)
                buffer[e][i] = right_conservatives_[e];
        }
    }
}