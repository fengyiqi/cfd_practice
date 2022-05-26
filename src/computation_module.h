#ifndef COMPUTATION_MODULE_H
#define COMPUTATION_MODULE_H

#include "block.h"
#include "solvers/hllc_riemann_solver.h"
#include "time_integration/runge_kutta_3.h"

class ComputationModule{
    Block& block_;
    HllcRiemannSolver& hllc_;
    RungeKutta3& runge_kutta_3_;
    double start_time_;
    
public:
    ComputationModule() = delete;
    ComputationModule(Block& block, HllcRiemannSolver& hllc, RungeKutta3& runge_kutta_3);
    void UpdateRightHandSide();
    void TimeIntegration();
    void Solve();
};

#endif