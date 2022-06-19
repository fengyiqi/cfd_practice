#include "runge_kutta_3.h"


void RungeKutta3::Advance(const unsigned int& stage, Block& block) {
        for (auto& e : FI::Equations){
            for (unsigned int i = GI::FICX(); i < GI::FRGX(); i++){
                block.conservative_buffer_next_[e][i] = coef_[stage][0] * block.conservative_buffer_[e][i] + 
                                                         coef_[stage][1] * block.conservative_buffer_next_[e][i] + 
                                                         coef_[stage][2] * block.rhs_hand_side_[e][i - GI::GCX()] * block.t_step_;
            }
        }

    }