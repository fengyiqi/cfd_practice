#ifndef HLLC_RIEMANN_SOLVER_H
#define HLLC_RIEMANN_SOLVER_H

#include "user_specification.h"
#include "stencils/weno5.h"

class HllcRiemannSolver {
    static constexpr unsigned int reconstruction_start_ = GI::FICX() - 3;
    static constexpr unsigned int reconstruction_end_   = reconstruction_start_ + GI::ICX() + 1;
    WENO5 weno5_;
    // temporary solution
    static constexpr double gamma_ = 1.4;
public:
    HllcRiemannSolver();
    ~HllcRiemannSolver() = default;
    double ConvertConservativesToFlux(const double (&conservatives)[FI::CN()], const unsigned int& state);
    void ComputeCellFaceFlux(const double (&conservatives)[FI::CN()][GI::TCX()], double (&flux)[FI::CN()][GI::ICX() + 1]);
};

#endif