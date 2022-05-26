#include "initial_condition.h"
#include "boundary_condition.h"
#include "stencils/weno5.h"
#include "block.h"
#include "computation_module.h"
#include "time_integration/runge_kutta_3.h"
#include "solvers/hllc_riemann_solver.h"

int main() {
    InitialCondition initial_condition = InitialCondition(0, 1, 0.2);
    Block block(initial_condition);
    HllcRiemannSolver hllc;
    RungeKutta3 runge_kutta_3;

    ComputationModule computation_module(block, hllc, runge_kutta_3);
    computation_module.Solve();


    // initial_condition.TestInitialCondition();
    // BoundaryCondition boundary_condition;
    // boundary_condition.TestBondaryCondition();
    // WENO5 weno5;
    // weno5.Test();

    return 0;
}