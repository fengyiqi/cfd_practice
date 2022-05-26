#ifndef BOUNDARY_CONDITION_H
#define BOUNDARY_CONDITION_H
#include "user_specification.h"

class BoundaryCondition {

public:
    BoundaryCondition() = default;
    ~BoundaryCondition() = default;
    void ApplyPeriodicBondaryCondition(double (&buffer)[FI::PN()][GI::TCX()]);
    void ApplySymmetricBondaryCondition(double (&buffer)[FI::PN()][GI::TCX()]);
    void TestBondaryCondition();
};

#endif