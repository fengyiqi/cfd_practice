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

    ConvertPrimitiveToConservativeStates(primitive_buffer_, conservative_buffer_);
    for (unsigned int cs = 0; cs < FI::CN(); cs++)
        for (unsigned int i = 0; i < GI::ICX(); i++) {
            rhs_hand_side_[cs][i] = 0.0;
            // conservative_buffer_next_[cs][i] = 0.0;
        }
    
    for (unsigned int cs = 0; cs < FI::CN(); cs++)
        for (unsigned int i = 0; i < GI::ICX() + 1; i++)
            cell_face_flux_temp_[cs][i] = 0.0;
}

void Block::ConvertPrimitiveToConservativeStates(const double (&primitives)[FI::PN()][GI::TCX()], double (&conservatives)[FI::CN()][GI::TCX()]) {
    for (unsigned int i = 0; i < GI::TCX(); i++) {
        conservatives[ConservativePool::Mass][i]      = primitives[PrimeStatePool::Density][i];
        conservatives[ConservativePool::MomentumX][i] = primitives[PrimeStatePool::Density][i] * primitives[PrimeStatePool::VelocityX][i];
        conservatives[ConservativePool::Energy][i]    = primitives[PrimeStatePool::Density][i] * 
                                                        (primitives[PrimeStatePool::Pressure][i] / (primitives[PrimeStatePool::Density][i] * (gamma_ - 1.0)) + 
                                                        0.5 * primitives[PrimeStatePool::VelocityX][i] * primitives[PrimeStatePool::VelocityX][i]);  
    }
}

void Block::ConvertConservativeToPrimitiveStates(const double (&conservatives)[FI::CN()][GI::TCX()], double (&primitives)[FI::PN()][GI::TCX()]) {
    for (unsigned int i = 0; i < GI::TCX(); i++) {
        primitives[PrimeStatePool::Density][i] = conservatives[ConservativePool::Mass][i];
        primitives[PrimeStatePool::VelocityX][i] = conservatives[ConservativePool::MomentumX][i] / conservatives[ConservativePool::Mass][i];
        primitives[PrimeStatePool::Pressure][i] = (gamma_ - 1.0) * (conservatives[ConservativePool::Energy][i] - 
                                                  0.5 * conservatives[ConservativePool::MomentumX][i] * conservatives[ConservativePool::MomentumX][i] / conservatives[ConservativePool::Mass][i]);
    }
}