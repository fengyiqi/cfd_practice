#include "block.h"

Block::Block(CaseSpecification& case_spec) :
cell_size_((case_spec.GetXEnd() - case_spec.GetXStart()) / GI::ICX()),
gamma_(case_spec.GetGamma()),
t_end_(case_spec.GetTEnd()),
t_step_(case_spec.GetTStep()),
boundary_condition_(case_spec.GetBoundaryCondition()),
eos_(case_spec.GetEquationOfState())
{
    case_spec.DefineCoordinate(x_coordinate_buffer_);
    case_spec.DefineInitialPrimitiveStates(x_coordinate_buffer_, primitive_buffer_);
    boundary_condition_.Apply(primitive_buffer_, FI::States::Primitives);
    ConvertPrimitiveToConservativeStates();

    for (auto& e : FI::Equations)
        for (unsigned int i = 0; i < GI::ICX(); i++) {
            rhs_hand_side_[e][i] = 0.0;
        }
    
    for (auto& e : FI::Equations)
        for (unsigned int i = 0; i < GI::ICX() + 1; i++)
            cell_face_flux_temp_[e][i] = 0.0;
}

void Block::ConvertPrimitiveToConservativeStates() {
    for (unsigned int i = 0; i < GI::TCX(); i++) {
        conservative_buffer_[EIndex(FI::EquationEnum::Mass)][i]      = primitive_buffer_[PIndex(FI::PrimeStateEnum::Density)][i];
        conservative_buffer_[EIndex(FI::EquationEnum::MomentumX)][i] = primitive_buffer_[PIndex(FI::PrimeStateEnum::Density)][i] * primitive_buffer_[PIndex(FI::PrimeStateEnum::VelocityX)][i];
        conservative_buffer_[EIndex(FI::EquationEnum::Energy)][i]    = primitive_buffer_[PIndex(FI::PrimeStateEnum::Density)][i] *
                                                        (primitive_buffer_[PIndex(FI::PrimeStateEnum::Pressure)][i] / (primitive_buffer_[PIndex(FI::PrimeStateEnum::Density)][i] * (gamma_ - 1.0)) +
                                                        0.5 * primitive_buffer_[PIndex(FI::PrimeStateEnum::VelocityX)][i] * primitive_buffer_[PIndex(FI::PrimeStateEnum::VelocityX)][i]);
    }
}

void Block::ConvertConservativeToPrimitiveStates() {
    for (unsigned int i = 0; i < GI::TCX(); i++) {
        primitive_buffer_[PIndex(FI::PrimeStateEnum::Density)][i] = conservative_buffer_[EIndex(FI::EquationEnum::Mass)][i];
        primitive_buffer_[PIndex(FI::PrimeStateEnum::VelocityX)][i] = conservative_buffer_[EIndex(FI::EquationEnum::MomentumX)][i] / conservative_buffer_[EIndex(FI::EquationEnum::Mass)][i];
        primitive_buffer_[PIndex(FI::PrimeStateEnum::Pressure)][i] = (gamma_ - 1.0) * (conservative_buffer_[EIndex(FI::EquationEnum::Energy)][i] -
                                                  0.5 * conservative_buffer_[EIndex(FI::EquationEnum::MomentumX)][i] * conservative_buffer_[EIndex(FI::EquationEnum::MomentumX)][i] / conservative_buffer_[EIndex(FI::EquationEnum::Mass)][i]);
    }
}