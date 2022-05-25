#include "initial_condition.h"
#include "boundary_condition.h"

int main() {
    InitialCondition initial_condition = InitialCondition(0, 1);
    initial_condition.TestInitialCondition();
    BoundaryCondition boundary_condition;
    boundary_condition.TestBondaryCondition();

    return 0;
}