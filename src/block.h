#ifndef BLOCK_H
#define BLOCK_H

#include "initial_condition.h"

class Block {
    double cell_size_;
    double gamma_;
    double primitive_buffer_[FI::PN()][GI::TCX()];
    double conservative_buffer_[FI::PN()][GI::TCX()];
    double x_coordinate_buffer_[GI::TCX()];

public:
    Block() = delete;
    Block(InitialCondition& initial_condition_);
    ~Block() = default;
};

#endif