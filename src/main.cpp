#include "initial_condition.h"
#include "boundary_condition.h"
#include "stencils/weno5.h"
#include "block.h"

int main() {
    InitialCondition initial_condition = InitialCondition(0, 1);
    Block block(initial_condition);

    // initial_condition.TestInitialCondition();
    // BoundaryCondition boundary_condition;
    // boundary_condition.TestBondaryCondition();
    // WENO5 weno5;
    // weno5.Test();

    return 0;
}