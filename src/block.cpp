#include "block.h"
#include <iostream>

Block::Block(std::shared_ptr<CaseSpecification>& case_spec, std::shared_ptr<BoundaryCondition>& boundary) :
cell_size_((case_spec->x_end - case_spec->x_start) / GI::ICX()),
gamma_(case_spec->gamma),
t_end_(case_spec->t_end),
t_step_(case_spec->t_step),
boundary_condition_(boundary)
{

    case_spec->DefineCoordinate(x_coordinate_buffer_);
    case_spec->DefineInitialPrimitiveStates(x_coordinate_buffer_, primitive_buffer_);
    boundary_condition_->Apply(primitive_buffer_);
    ConvertPrimitiveToConservativeStates();
//    boundary_condition_->Apply(conservative_buffer_);

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