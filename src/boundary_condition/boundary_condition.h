#ifndef BOUNDARY_CONDITION_H
#define BOUNDARY_CONDITION_H
#include "user_specification.h"

class BoundaryCondition {
public:
    BoundaryCondition() = default;
    virtual ~BoundaryCondition() = default;
    virtual void Apply(double (&buffer)[FI::PN()][GI::TCX()], unsigned int state = States::Primitives) = 0;
};

#endif