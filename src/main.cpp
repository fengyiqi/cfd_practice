#include "block.h"

#include "solvers/riemann_solvers.h"
#include "output.h"
#include "initial_condition/case_parameters.h"
#include "stencils/stencil_instance.h"
#include "computation_module.h"

int main()
{   
    {
        SodShockTube sod_shock_tube;
        Block sod_shock_tube_block(sod_shock_tube);
        ComputationModule<HllcRiemannSolver<WENO5>> computation_module(sod_shock_tube_block);
        computation_module.Solve();
        OutputWriter csv_writer(sod_shock_tube_block, "sod.csv");
        csv_writer.WriteCSV();
    }

    {
        ShuOsher shu_osher;
        Block shu_osher_block(shu_osher);
        ComputationModule<HllcRiemannSolver<WENO5>> computation_module(shu_osher_block);
        computation_module.Solve();
        OutputWriter csv_writer(shu_osher_block, "shu.csv");
        csv_writer.WriteCSV();
    }
    return 0;
}