#ifndef RUSANOV_RIEMANN_SOLVER_H
#define RUSANOV_RIEMANN_SOLVER_H

#include "user_specification.h"
#include "stencils/weno5.h"
#include "stencils/first_order.h"
#include "block.h"

class RusanovRiemannSolver {
    WENO5 stencil_;
    Block& block_;
public:
    RusanovRiemannSolver() = delete;
    RusanovRiemannSolver(Block& block) : block_(block) {};
    ~RusanovRiemannSolver() = default;
    double ConvertConservativesToFlux(const double (&conservatives)[FI::EN()], const std::underlying_type<FI::EquationEum>::type equation);
    void ComputeCellFaceFlux(const double (&conservatives)[FI::EN()][GI::TCX()], double (&flux)[FI::EN()][GI::ICX() + 1]);
};

#endif