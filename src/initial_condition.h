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
    double t_end_;
    double t_step_ = 0.0001;
    BoundaryCondition boundary_condition_;
public:
    InitialCondition() = delete;
    InitialCondition(const double& x_start, const double& x_end, const double& t_end);
    void DefineCoordinate(double (&x)[GI::TCX()]);
    void DefineInitialPrimitiveStates(const double (&x)[GI::TCX()], double (&primitives)[FI::PN()][GI::TCX()]);
    void TestInitialCondition();
    inline double GetCellSize() const { return cell_size_; }
    inline double GetXStart() const { return start_; }
    inline double GetXEnd() const { return end_; }
    inline double GetGamma() const { return gamma_; }
    inline double GetTEnd() const { return t_end_; }
    inline double GetTStep() const { return t_step_; }
    inline BoundaryCondition& GetBoundaryCondition() { return boundary_condition_;}
};

#endif