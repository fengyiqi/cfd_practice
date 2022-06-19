#ifndef FIXED_VALUE_H
#define FIXED_VALUE_H

#include "boundary_condition.h"

class FixedValueBoundaryCondition : public BoundaryCondition {
    double left_primitives_[FI::PN()];
    double right_primitives_[FI::PN()];
    double left_conservatives_[FI::EN()];
    double right_conservatives_[FI::EN()];
public:
    FixedValueBoundaryCondition() = delete;
    explicit FixedValueBoundaryCondition(const double (&left_side)[FI::PN()], const double (&right_side)[FI::PN()]);
    ~FixedValueBoundaryCondition() = default;
    void Apply(double (&buffer)[FI::PN()][GI::TCX()], FI::States state = FI::States::Primitives) override;
};

#endif