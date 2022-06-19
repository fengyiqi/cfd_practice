//
// Created by yiqif on 2022/5/27.
//

#ifndef SYMMETRIC_H
#define SYMMETRIC_H

#include "boundary_condition.h"

class SymmetricBoundaryCondition : public BoundaryCondition {
public:
    SymmetricBoundaryCondition() = default;
    ~SymmetricBoundaryCondition() = default;
    void Apply(double (&buffer)[FI::PN()][GI::TCX()], FI::States state = FI::States::Primitives) override;
};

#endif
