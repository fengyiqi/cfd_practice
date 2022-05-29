#ifndef BLOCK_H
#define BLOCK_H

#include "initial_condition/case_parameters.h"
#include "user_specification.h"
#include "boundary_condition/boundary_condition.h"
#include <memory>

class Block {
public:
    double cell_size_;
    double gamma_;
    double t_end_;
    double t_step_;
    double primitive_buffer_[FI::PN()][GI::TCX()];
    double conservative_buffer_[FI::CN()][GI::TCX()];
    double conservative_buffer_next_[FI::CN()][GI::TCX()];
    double x_coordinate_buffer_[GI::TCX()];

    double rhs_hand_side_[FI::CN()][GI::ICX()];
    double cell_face_flux_temp_[FI::CN()][GI::ICX()+1];
    BoundaryCondition& boundary_condition_;

public:
    Block() = delete;
    explicit Block(CaseSpecification& case_spec);
    ~Block() = default;
    void ConvertPrimitiveToConservativeStates();
    void ConvertConservativeToPrimitiveStates();
};

#endif