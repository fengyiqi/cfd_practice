#include "initial_condition.h"
#include "boundary_condition.h"
#include "stencils/weno5.h"
#include "block.h"
#include "computation_module.h"
#include "time_integration/runge_kutta_3.h"
#include "solvers/hllc_riemann_solver.h"
#include "solvers/rusanov_riemann_solver.h"
#include "output.h"
#include <iostream>
#include <string>

int main() {
    
    InitialCondition initial_condition = InitialCondition(0, 1, 0.2);
    Block block(initial_condition);
    HllcRiemannSolver riemann_solver;
    RungeKutta3 runge_kutta_3;

    ComputationModule computation_module(block, riemann_solver, runge_kutta_3);
    computation_module.Solve();

    OutputWriter csv_writer(block, "sod.csv");
    csv_writer.WriteCSV();

    return 0;
}