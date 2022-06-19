#ifndef HLLC_RIEMANN_SOLVER_H
#define HLLC_RIEMANN_SOLVER_H

#include "user_specification.h"
#include "stencils/weno5.h"
#include "stencils/first_order.h"
#include "block.h"


class HllcRiemannSolver {
    FirstOrder stencil_;
    Block& block_;
public:
    HllcRiemannSolver() = delete;
    HllcRiemannSolver(Block& block) : block_(block) {};
    ~HllcRiemannSolver() = default;
    void ComputeCellFaceFlux();
};

#endif