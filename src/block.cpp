#include "block.h"
#include <iostream>

Block::Block(InitialCondition& init_cond) : 
                                            initial_condition_(init_cond),
                                            cell_size_(init_cond.GetCellSize()),
                                            gamma_(init_cond.GetGamma()),
                                            t_end_(init_cond.GetTEnd()),
                                            t_step_(init_cond.GetTStep()),
                                            boundary_condition_(init_cond.GetBoundaryCondition())
                                            
{
    initial_condition_.DefineCoordinate(x_coordinate_buffer_);
    initial_condition_.DefineInitialPrimitiveStates(x_coordinate_buffer_, primitive_buffer_);
    initial_condition_.ConvertToConservativeState(primitive_buffer_, conservative_buffer_);

    for (unsigned int cs = 0; cs < FI::CN(); cs++)
        for (unsigned int i = 0; i < GI::ICX(); i++) {
            rhs_hand_side_[cs][i] = 0.0;
            conservative_buffer_next_[cs][i] = 0.0;
        }
    
    for (unsigned int cs = 0; cs < FI::CN(); cs++)
        for (unsigned int i = 0; i < GI::ICX() + 1; i++)
            cell_face_flux_temp_[cs][i] = 0.0;
    std::cout << t_end_ << std::endl;
}