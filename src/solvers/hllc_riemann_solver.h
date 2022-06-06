#ifndef HLLC_RIEMANN_SOLVER_H
#define HLLC_RIEMANN_SOLVER_H

#include "user_specification.h"
#include "stencils/weno5.h"
#include "block.h"

class HllcRiemannSolver {
    static constexpr unsigned int reconstruction_start_ = GI::FICX() - 3;
    WENO5 weno5_;
    Block& block_;
public:
    HllcRiemannSolver() = delete;
    HllcRiemannSolver(Block& block) : block_(block) {};
    ~HllcRiemannSolver() = default;
    void ComputeCellFaceFlux();
};

#endif