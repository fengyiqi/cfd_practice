#include "boundary_condition/boundary_condition.h"
#include "stencils/weno5.h"
#include "block.h"
#include "computation_module.h"
#include "time_integration/runge_kutta_3.h"
#include "solvers/hllc_riemann_solver.h"
#include "solvers/rusanov_riemann_solver.h"
#include "output.h"
#include <iostream>
#include <string>
#include <memory>

int main() {
    std::shared_ptr<CaseSpecification> sod_shock_tube = std::make_shared<SodShockTube>();
    std::shared_ptr<BoundaryCondition> boundary_condition = std::make_shared<SymmetricBoundaryCondition>();
    Block block(sod_shock_tube, boundary_condition);
    HllcRiemannSolver riemann_solver;
    RungeKutta3 runge_kutta_3;

    ComputationModule computation_module(block, riemann_solver, runge_kutta_3);
    computation_module.Solve();

    OutputWriter csv_writer(block, "sod.csv");
    csv_writer.WriteCSV();

    return 0;
}