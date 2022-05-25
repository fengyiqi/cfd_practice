#include "block.h"

Block::Block(InitialCondition& initial_condition_) : 
                                                    cell_size_(initial_condition_.GetCellSize()),
                                                    gamma_(initial_condition_.GetGamma()) 
{
    initial_condition_.DefineCoordinate(x_coordinate_buffer_);
    initial_condition_.DefineInitialPrimitiveStates(x_coordinate_buffer_, primitive_buffer_);
    initial_condition_.ConvertToConservativeState(primitive_buffer_, conservative_buffer_);
}