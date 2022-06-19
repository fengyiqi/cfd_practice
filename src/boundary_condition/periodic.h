//
// Created by yiqif on 2022/5/27.
//

#ifndef PERIODIC_H
#define PERIODIC_H

#include "boundary_condition.h"

class PeriodicBoundaryCondition : public BoundaryCondition {
public:
    PeriodicBoundaryCondition() = default;
    ~PeriodicBoundaryCondition() = default;
    void Apply(double (&buffer)[FI::PN()][GI::TCX()], FI::States state = FI::States::Primitives) override;
};

#endif
