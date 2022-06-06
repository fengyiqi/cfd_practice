#include "boundary_condition/boundary_condition.h"
#include "stencils/weno5.h"
#include "block.h"
#include "computation_module.h"
#include "time_integration/runge_kutta_3.h"
#include "solvers/hllc_riemann_solver.h"
#include "solvers/rusanov_riemann_solver.h"
#include "output.h"
#include "initial_condition/case_parameters.h"

#include <iostream>
#include <string>
#include <memory>

int main()
{   
    RungeKutta3 runge_kutta_3;
    {
        

        SodShockTube sod_shock_tube;
        Block sod_shock_tube_block(sod_shock_tube);
        HllcRiemannSolver riemann_solver(sod_shock_tube_block);
        ComputationModule computation_module(sod_shock_tube_block, riemann_solver, runge_kutta_3);
        computation_module.Solve();
        OutputWriter csv_writer(sod_shock_tube_block, "sod.csv");
        csv_writer.WriteCSV();
    }

    {
        ShuOsher shu_osher;
        Block shu_osher_block(shu_osher);
        HllcRiemannSolver riemann_solver(shu_osher_block);
        ComputationModule computation_module(shu_osher_block, riemann_solver, runge_kutta_3);
        computation_module.Solve();
        OutputWriter csv_writer(shu_osher_block, "shu.csv");
        csv_writer.WriteCSV();
    }
    return 0;
}