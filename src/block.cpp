#include "block.h"
#include <iostream>

Block::Block(CaseSpecification& case_spec) :
cell_size_((case_spec.x_end_ - case_spec.x_start_) / GI::ICX()),
gamma_(case_spec.GetGamma()),
t_end_(case_spec.GetTEnd()),
t_step_(case_spec.GetTStep()),
boundary_condition_(case_spec.GetBoundaryCondition()),
eos_(case_spec.GetEquationOfState())
{
    case_spec.DefineCoordinate(x_coordinate_buffer_);
    case_spec.DefineInitialPrimitiveStates(x_coordinate_buffer_, primitive_buffer_);
    boundary_condition_.Apply(primitive_buffer_, States::Primitives);
    ConvertPrimitiveToConservativeStates();

    for (unsigned int cs = 0; cs < FI::CN(); cs++)
        for (unsigned int i = 0; i < GI::ICX(); i++) {
            rhs_hand_side_[cs][i] = 0.0;
        }
    
    for (unsigned int cs = 0; cs < FI::CN(); cs++)
        for (unsigned int i = 0; i < GI::ICX() + 1; i++)
            cell_face_flux_temp_[cs][i] = 0.0;
}

void Block::ConvertPrimitiveToConservativeStates() {
    for (unsigned int i = 0; i < GI::TCX(); i++) {
        conservative_buffer_[ConservativePool::Mass][i]      = primitive_buffer_[PrimeStatePool::Density][i];
        conservative_buffer_[ConservativePool::MomentumX][i] = primitive_buffer_[PrimeStatePool::Density][i] * primitive_buffer_[PrimeStatePool::VelocityX][i];
        conservative_buffer_[ConservativePool::Energy][i]    = primitive_buffer_[PrimeStatePool::Density][i] *
                                                        (primitive_buffer_[PrimeStatePool::Pressure][i] / (primitive_buffer_[PrimeStatePool::Density][i] * (gamma_ - 1.0)) +
                                                        0.5 * primitive_buffer_[PrimeStatePool::VelocityX][i] * primitive_buffer_[PrimeStatePool::VelocityX][i]);
    }
}

void Block::ConvertConservativeToPrimitiveStates() {
    for (unsigned int i = 0; i < GI::TCX(); i++) {
        primitive_buffer_[PrimeStatePool::Density][i] = conservative_buffer_[ConservativePool::Mass][i];
        primitive_buffer_[PrimeStatePool::VelocityX][i] = conservative_buffer_[ConservativePool::MomentumX][i] / conservative_buffer_[ConservativePool::Mass][i];
        primitive_buffer_[PrimeStatePool::Pressure][i] = (gamma_ - 1.0) * (conservative_buffer_[ConservativePool::Energy][i] -
                                                  0.5 * conservative_buffer_[ConservativePool::MomentumX][i] * conservative_buffer_[ConservativePool::MomentumX][i] / conservative_buffer_[ConservativePool::Mass][i]);
    }
}