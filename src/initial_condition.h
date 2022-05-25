#ifndef INITIAL_CONDITION_H
#define INITIAL_CONDITION_H

#include <array>
#include "user_specification.h"
#include "boundary_condition.h"

class InitialCondition {
    double start_;
    double end_;
    double length_;
    double cell_size_;
    double gamma_ = 1.4;
    BoundaryCondition boundary_condition_;
    // unsigned int cells = GI::TCX();
public:
    InitialCondition() = delete;
    InitialCondition(const double& x_start, const double& x_end);
    // InitialCondition()
    inline double GetGamma() const { return gamma_; }
    void DefineCoordinate(double (&x)[GI::TCX()]);
    void DefineInitialPrimitiveStates(const double (&x)[GI::TCX()], double (&primitives)[FI::PN()][GI::TCX()]);
    void ConvertToConservativeState(const double (&primitives)[FI::PN()][GI::TCX()], double (&conservatives)[FI::CN()][GI::TCX()]);
    void TestInitialCondition();
};

#endif