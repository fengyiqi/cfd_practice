#ifndef BOUNDARY_CONDITION_H
#define BOUNDARY_CONDITION_H
#include "user_specification.h"

class BoundaryCondition {
public:
    BoundaryCondition() = default;
    ~BoundaryCondition() = default;
    virtual void Apply(double (&buffer)[FI::PN()][GI::TCX()]) = 0;
};

#endif