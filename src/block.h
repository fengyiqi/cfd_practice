#ifndef BLOCK_H
#define BLOCK_H

#include "initial_condition.h"
// #include "computation_module.h"

class Block {
    // friend RungeKutta3;
    // friend ComputationModule;
public:
    double cell_size_;
    double gamma_;
    double t_end_;
    double t_step_;
    double primitive_buffer_[FI::PN()][GI::TCX()];
    double conservative_buffer_[FI::PN()][GI::TCX()];
    double conservative_buffer_next_[FI::PN()][GI::TCX()];
    double x_coordinate_buffer_[GI::TCX()];
    double rhs_hand_side_[FI::CN()][GI::ICX()];
    double cell_face_flux_temp_[FI::CN()][GI::ICX()+1];

    InitialCondition& initial_condition_;
    BoundaryCondition& boundary_condition_;

public:
    Block() = delete;
    Block(InitialCondition& init_cond);
    ~Block() = default;
};

#endif