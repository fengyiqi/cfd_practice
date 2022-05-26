#ifndef BLOCK_H
#define BLOCK_H

#include "initial_condition.h"

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

    InitialCondition& initial_condition_;
    BoundaryCondition& boundary_condition_;

public:
    Block() = delete;
    Block(InitialCondition& init_cond);
    ~Block() = default;
    void ConvertPrimitiveToConservativeStates(const double (&primitives)[FI::PN()][GI::TCX()], double (&conservatives)[FI::CN()][GI::TCX()]);
    void ConvertConservativeToPrimitiveStates(const double (&conservatives)[FI::CN()][GI::TCX()], double (&primitives)[FI::PN()][GI::TCX()]);
    void ConvertPrimitiveToConservativeStates() { ConvertPrimitiveToConservativeStates(primitive_buffer_, conservative_buffer_); }
    void ConvertConservativeToPrimitiveStates() { ConvertConservativeToPrimitiveStates(conservative_buffer_, primitive_buffer_); }
};

#endif